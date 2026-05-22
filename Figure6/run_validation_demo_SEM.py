#!/usr/bin/env python3
"""
run_validation_demo_SEM.py
==========================

回应审稿人 "feature selection leakage" concern 的简化管道，SEM-only。

只跑两件事 + 一个对比 + 一个汇总 Excel：
  [1] Repeated stratified 5-fold CV (panel-locked)
      → paper-aligned robust AUC ± 95% CI
  [2] Permutation test (fold-internal RFE, n=100)
      → 严格条件下的 empirical p value
  [3] Performance comparison summary
      → 把 [1] 的真实 AUC 分布和 [2] 的 null 分布并排展示
  [4] Compile all results into a multi-sheet Excel file
      → 直接可用作 Supplementary Table

输出：
  validation_results_SEM/
    fixed_panel_20features.csv
    repeated_kfold_cv_SEM.csv
    repeated_kfold_cv_SEM_folds.csv
    permutation_test_SEM.csv
    permutation_test_SEM_null_dist.csv
    comparison_summary.csv
    Supplementary_Table_ProMeta-BC_validation.xlsx   ← 多 sheet 汇总
"""
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier, StackingClassifier
from sklearn.neural_network import MLPClassifier

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.utils.dataframe import dataframe_to_rows

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validation_methods import (
    repeated_kfold_cv,
    permutation_test_cv,
)

OUTDIR = "./validation_results_SEM"
os.makedirs(OUTDIR, exist_ok=True)


def build_sem_only(random_state=40):
    glm = LogisticRegression(penalty="elasticnet", l1_ratio=0.5, solver="saga",
                             max_iter=3000, random_state=random_state)
    lr = LogisticRegression(penalty="l2", solver="liblinear",
                            max_iter=3000, random_state=random_state)
    svm = SVC(kernel="linear", probability=True, random_state=random_state)
    gbm = GradientBoostingClassifier(learning_rate=0.01, n_estimators=200,
                                     max_depth=1, subsample=0.8, random_state=random_state)
    mlp = MLPClassifier(learning_rate_init=0.001, hidden_layer_sizes=(64, 32),
                        activation="relu", alpha=50, solver="adam",
                        max_iter=2000, random_state=random_state)
    sem = StackingClassifier(
        estimators=[("glm", glm), ("lr", lr), ("svm", svm), ("mlp", mlp), ("gbm", gbm)],
        final_estimator=LogisticRegression(penalty="l2", solver="liblinear",
                                           max_iter=3000, random_state=random_state),
        passthrough=False, n_jobs=1,
    )
    return {"SEM": sem}


def load_cohort1():
    pro = pd.read_csv("pro_coh1_filtered.csv").set_index("Gene").T
    meta = pd.read_csv("Meta_coh1_filtered.csv").set_index("Metabolite").T
    clinical = pd.read_csv("Clinical_coh1.csv").set_index("Sample")

    sp = StandardScaler(); sm = StandardScaler()
    pro_s = pd.DataFrame(sp.fit_transform(pro), index=pro.index, columns=pro.columns)
    meta_s = pd.DataFrame(sm.fit_transform(meta), index=meta.index, columns=meta.columns)

    df = clinical.join(pro_s, how="inner").join(meta_s, how="inner")
    feature_cols = list(pro_s.columns) + list(meta_s.columns)
    X = df[feature_cols].fillna(0)
    y = df["Label"]
    return X, y


def select_fixed_panel(X_train, y_train, n_features=20, random_state=40):
    """在整个 training set 上 RFE 选出固定 panel（= paper 工作流）。"""
    y_bin = (y_train == "Malignant").values.astype(int)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_train)
    rfe = RFE(estimator=SVC(kernel="linear", random_state=random_state),
              n_features_to_select=n_features, step=10)
    rfe.fit(X_scaled, y_bin)
    return X_train.columns[rfe.support_].tolist()


# ============================================================
# Excel export helpers
# ============================================================
HEADER_FILL = PatternFill("solid", start_color="2E86AB")
HEADER_FONT = Font(name="Arial", size=11, bold=True, color="FFFFFF")
TITLE_FONT  = Font(name="Arial", size=13, bold=True, color="1F1F1F")
NOTE_FONT   = Font(name="Arial", size=10, italic=True, color="555555")
CELL_FONT   = Font(name="Arial", size=10)
THIN        = Side(style="thin", color="CCCCCC")
BORDER      = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)


def _write_df(ws, df, start_row=1):
    """Write DataFrame to ws starting at start_row with styled header."""
    for j, col in enumerate(df.columns, start=1):
        c = ws.cell(row=start_row, column=j, value=col)
        c.font = HEADER_FONT; c.fill = HEADER_FILL
        c.alignment = Alignment(horizontal="center", vertical="center")
        c.border = BORDER
    for i, row in enumerate(dataframe_to_rows(df, index=False, header=False),
                            start=start_row + 1):
        for j, val in enumerate(row, start=1):
            c = ws.cell(row=i, column=j, value=val)
            c.font = CELL_FONT
            c.border = BORDER
            c.alignment = Alignment(horizontal="left" if isinstance(val, str) else "right",
                                    vertical="center")
    for j, col in enumerate(df.columns, start=1):
        max_len = max([len(str(col))] + [len(str(v)) for v in df.iloc[:, j - 1].astype(str)])
        ws.column_dimensions[get_column_letter(j)].width = min(max(max_len + 2, 12), 50)
    return start_row + len(df) + 1


def export_supplementary_xlsx(outdir, final_panel, rep_summary, rep_folds,
                              perm_summary, perm_dist, comparison):
    """把所有结果汇总到一个多 sheet 的 xlsx。"""
    out_xlsx = os.path.join(outdir, "Supplementary_Table_ProMeta-BC_validation.xlsx")
    wb = Workbook(); wb.remove(wb.active)

    # ----- Sheet 1: README -----
    ws = wb.create_sheet("README")
    ws["A1"] = "Supplementary Table — ProMeta-BC validation results"
    ws["A1"].font = TITLE_FONT
    ws.merge_cells("A1:B1")
    readme_rows = [
        ("", ""),
        ("Sheet", "Content"),
        ("Final_Panel",
         "20 features (15 proteins + 5 metabolites) selected by SVM-RFE on the full Cohort 2 training set; used as the panel-locked input for ProMeta-BC."),
        ("Repeated_CV_Summary",
         "Summary statistics of repeated stratified 5-fold cross-validation (5 folds × 10 random seeds = 50 folds) of the panel-locked SEM model on the Cohort 2 training set."),
        ("Repeated_CV_Folds",
         "All 50 fold-level AUC values from repeated cross-validation."),
        ("Permutation_Summary",
         "Summary statistics of the label-permutation test (n = 100). For each permutation, malignant/benign labels were shuffled and the full SVM-RFE + standardization + SEM training pipeline was re-executed."),
        ("Permutation_Null",
         "All 100 permutation AUC values forming the null distribution."),
        ("Comparison",
         "Side-by-side numerical comparison of the two robustness tests."),
        ("", ""),
        ("Methodology notes:", ""),
        ("• Repeated CV uses the panel-locked 20-feature panel from the full training set (paper-aligned).", ""),
        ("• Permutation test uses fold-internal RFE for each permutation, providing a stricter null reference.", ""),
        ("• Both analyses use the same SEM (stacking ensemble) implementation as the main ProMeta-BC model.", ""),
    ]
    for r, (a, b) in enumerate(readme_rows, start=2):
        ca = ws.cell(row=r, column=1, value=a); cb = ws.cell(row=r, column=2, value=b)
        if r == 3:
            ca.font = HEADER_FONT; ca.fill = HEADER_FILL
            cb.font = HEADER_FONT; cb.fill = HEADER_FILL
            ca.alignment = Alignment(horizontal="center")
            cb.alignment = Alignment(horizontal="center")
        elif a.startswith("•") or a.startswith("Methodology"):
            ca.font = NOTE_FONT
        else:
            ca.font = CELL_FONT; cb.font = CELL_FONT
            ca.alignment = Alignment(vertical="top")
            cb.alignment = Alignment(wrap_text=True, vertical="top")
    ws.column_dimensions["A"].width = 26
    ws.column_dimensions["B"].width = 90

    # ----- Sheet 2: Final panel -----
    ws = wb.create_sheet("Final_Panel")
    ws["A1"] = "Final ProMeta-BC feature panel (20 features)"
    ws["A1"].font = TITLE_FONT
    fp_df = pd.DataFrame({"rank": np.arange(1, len(final_panel) + 1),
                          "feature": final_panel})
    _write_df(ws, fp_df, start_row=3)

    # ----- Sheet 3-7: results -----
    for sheet_name, title, df in [
        ("Repeated_CV_Summary", "Repeated 5-fold CV (panel-locked SEM) — summary", rep_summary),
        ("Repeated_CV_Folds",   "Repeated 5-fold CV — fold-level AUC values",       rep_folds),
        ("Permutation_Summary", "Permutation test (fold-internal RFE) — summary",   perm_summary),
        ("Permutation_Null",    "Permutation test — null distribution of CV AUC",   perm_dist),
        ("Comparison",          "Performance comparison — Repeated CV vs Permutation test", comparison),
    ]:
        ws = wb.create_sheet(sheet_name)
        ws["A1"] = title
        ws["A1"].font = TITLE_FONT
        _write_df(ws, df, start_row=3)

    wb.save(out_xlsx)
    return out_xlsx


def main():
    print("=" * 78)
    print("SEM-only validation pipeline (simplified)")
    print("Repeated CV (panel-locked) + Permutation test (fold-internal RFE)")
    print("=" * 78)

    X_all, y_all = load_cohort1()
    print(f"Cohort: {X_all.shape[0]} samples, {X_all.shape[1]} candidate features")
    print(f"Label: {y_all.value_counts().to_dict()}")

    X_train, X_test, y_train, y_test = train_test_split(
        X_all, y_all, test_size=0.3, random_state=46, stratify=y_all)
    print(f"\nTraining set: {X_train.shape[0]} samples "
          f"({(y_train == 'Malignant').sum()} Mal / {(y_train == 'Benign').sum()} Ben)")

    N_FEATURES = 20
    sem_models = build_sem_only()

    # ============================================================
    # Step 0: Paper 的 final panel (在 full training set 上 RFE)
    # ============================================================
    print(f"\n>>> Step 0: select FINAL panel from full training set ({N_FEATURES} features)")
    final_panel = select_fixed_panel(X_train, y_train, n_features=N_FEATURES)
    print(f"Final panel ({len(final_panel)} features):")
    for i, f in enumerate(final_panel, 1):
        print(f"  {i:2d}. {f}")
    pd.DataFrame({"feature": final_panel}).to_csv(
        os.path.join(OUTDIR, "fixed_panel_20features.csv"), index=False)

    # ============================================================
    # [1] Repeated CV with FIXED panel
    # ============================================================
    print("\n" + "=" * 78)
    print("[1] Repeated CV with FIXED panel (paper-aligned robust AUC)")
    print("=" * 78)
    rep_summary, rep_folds = repeated_kfold_cv(
        X_train, y_train,
        models=sem_models,
        fixed_features=final_panel,           # ✅ panel-locked
        n_splits=5, n_repeats=10,
        random_state=42,
        output_csv=os.path.join(OUTDIR, "repeated_kfold_cv_SEM.csv"),
    )
    print("\n--- Repeated CV summary ---")
    print(rep_summary.to_string(index=False))

    # ============================================================
    # [2] Permutation test (fold-internal RFE)
    # ============================================================
    N_PERM = 100   # 正式投稿可改成 1000
    print("\n" + "=" * 78)
    print(f"[2] Permutation test (fold-internal RFE, n={N_PERM})")
    print("=" * 78)
    perm_summary, perm_dist = permutation_test_cv(
        X_train, y_train,
        models=sem_models,
        fixed_features=None,                  # ✅ permutation 用 fold-internal
        n_features=N_FEATURES, n_splits=5,
        n_permutations=N_PERM,
        random_state=42,
        output_csv=os.path.join(OUTDIR, "permutation_test_SEM.csv"),
    )
    print("\n--- Permutation test summary ---")
    print(perm_summary.to_string(index=False))

    # ============================================================
    # [3] Performance comparison
    # ============================================================
    print("\n" + "=" * 78)
    print("[3] Performance comparison: Repeated CV vs Permutation test")
    print("=" * 78)

    rep_mean = rep_summary.iloc[0]["Mean_AUC"]
    rep_std  = rep_summary.iloc[0]["Std_AUC"]
    rep_lo, rep_hi = rep_summary.iloc[0]["CI95_low"], rep_summary.iloc[0]["CI95_high"]

    perm_obs       = perm_summary.iloc[0]["True_CV_AUC"]
    perm_null_mean = perm_summary.iloc[0]["Null_Mean"]
    perm_null_std  = perm_summary.iloc[0]["Null_Std"]
    perm_p         = perm_summary.iloc[0]["p_value"]
    n_sds          = (perm_obs - perm_null_mean) / perm_null_std if perm_null_std > 0 else np.nan

    comparison = pd.DataFrame({
        "Test": (["Repeated 5-fold CV (panel-locked)"] * 4
                 + ["Permutation test (fold-internal RFE)"] * 7),
        "Metric": [
            "Mean AUC (true labels)",
            "SD",
            "95% empirical CI",
            "N folds",
            "Observed CV AUC (true labels)",
            "Null mean",
            "Null SD",
            "Delta (observed - null mean)",
            "Delta in SDs of null",
            "Empirical p",
            "N permutations",
        ],
        "Value": [
            f"{rep_mean:.3f}",
            f"{rep_std:.3f}",
            f"[{rep_lo:.3f}, {rep_hi:.3f}]",
            f"{len(rep_folds)}",
            f"{perm_obs:.3f}",
            f"{perm_null_mean:.3f}",
            f"{perm_null_std:.3f}",
            f"{perm_obs - perm_null_mean:+.3f}",
            f"{n_sds:.1f}",
            f"{perm_p:.4f}",
            f"{N_PERM}",
        ],
    })
    comparison.to_csv(os.path.join(OUTDIR, "comparison_summary.csv"), index=False)
    print(comparison.to_string(index=False))

    # ============================================================
    # [4] Export consolidated Excel for Supplementary Table
    # ============================================================
    print("\n" + "=" * 78)
    print("[4] Exporting multi-sheet Excel for Supplementary Table")
    print("=" * 78)
    xlsx_path = export_supplementary_xlsx(
        OUTDIR, final_panel, rep_summary, rep_folds,
        perm_summary, perm_dist, comparison,
    )
    print(f"✅ {xlsx_path}")

    # ============================================================
    # 总览：response letter 用的数字
    # ============================================================
    print("\n" + "=" * 78)
    print("SUMMARY for response letter")
    print("=" * 78)
    print(f"  Final panel:               {len(final_panel)} features")
    print(f"  [1] Repeated CV:           mean AUC = {rep_mean:.3f} ± {rep_std:.3f} "
          f"(95% CI [{rep_lo:.3f}, {rep_hi:.3f}], N = {len(rep_folds)} folds)")
    print(f"  [2] Permutation test:      observed AUC = {perm_obs:.3f}, "
          f"null mean = {perm_null_mean:.3f} (Δ = {n_sds:.1f} SDs), p = {perm_p:.4f}")

    print(f"\nAll outputs in: {OUTDIR}/")


if __name__ == "__main__":
    main()
