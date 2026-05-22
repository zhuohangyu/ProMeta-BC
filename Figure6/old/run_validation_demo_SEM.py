#!/usr/bin/env python3
"""
run_validation_demo_SEM.py
==========================

回应审稿人 "feature selection leakage" concern 的完整管道，SEM-only。

设计思路（双轨制）：
  Track A (panel-locked, 与 paper 一致):
    - Step 1: 在整个 training set 上 RFE 选出固定 20-feature panel（这是 paper 报告的 ProMeta-BC）
    - Repeated CV with fixed panel → 给出 paper 风格的 robust AUC ± 95% CI

  Track B (fold-internal, 回应 leakage):
    - Nested CV with fold-internal RFE → 无偏估计，证明无 leakage
    - Permutation test with fold-internal RFE → 严格条件下的 p 值
    - Feature stability matrix → 证明 final panel 的特征是稳定的

输出：
  validation_results_SEM/
    fixed_panel_20features.csv         # paper 的 final panel
    repeated_kfold_cv_SEM.csv          # Track A: panel-locked AUC ± CI
    nested_cv_SEM.csv                  # Track B: 无偏 AUC ± CI
    nested_cv_SEM_folds.csv            # 含每个 outer fold 选了哪些特征
    permutation_test_SEM.csv           # Track B: p value
    feature_stability_matrix.csv       # Track B: 特征 × fold 矩阵
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validation_methods import (
    repeated_kfold_cv,
    nested_cv_evaluation,
    permutation_test_cv,
    compute_feature_stability,
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


def main():
    print("=" * 78)
    print("SEM-only validation pipeline (dual-track)")
    print("=" * 78)

    X_all, y_all = load_cohort1()
    print(f"Cohort 1: {X_all.shape[0]} samples, {X_all.shape[1]} candidate features")
    print(f"Label: {y_all.value_counts().to_dict()}")

    X_train, X_test, y_train, y_test = train_test_split(
        X_all, y_all, test_size=0.3, random_state=46, stratify=y_all)
    print(f"\nTraining set: {X_train.shape[0]} samples "
          f"({(y_train == 'Malignant').sum()} Mal / {(y_train == 'Benign').sum()} Ben)")

    N_FEATURES = 20
    candidate_features = list(X_train.columns)
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
    # TRACK A: Panel-locked Repeated CV (与 paper AUC 对齐)
    # ============================================================
    print("\n" + "=" * 78)
    print("[Track A] Repeated CV with FIXED panel — paper-aligned robust AUC")
    print("=" * 78)
    rep_summary, rep_folds = repeated_kfold_cv(
        X_train, y_train,
        models=sem_models,
        fixed_features=final_panel,           # ✅ panel-locked
        n_splits=5, n_repeats=10,
        random_state=42,
        output_csv=os.path.join(OUTDIR, "repeated_kfold_cv_SEM.csv"),
    )
    print("\n--- Track A summary (paper-style) ---")
    print(rep_summary.to_string(index=False))

    # ============================================================
    # TRACK B: Fold-internal nested CV (回应 leakage concern)
    # ============================================================
    print("\n" + "=" * 78)
    print("[Track B-1] Nested CV with FOLD-INTERNAL RFE — unbiased, addresses leakage")
    print("=" * 78)
    nest_summary, nest_folds = nested_cv_evaluation(
        X_train, y_train,
        n_features=N_FEATURES,
        fixed_features=None,                  # ✅ 关键：每个 fold 内重做 RFE
        outer_splits=5, inner_splits=3,
        random_state=42,
        output_csv=os.path.join(OUTDIR, "nested_cv_SEM.csv"),
        models_subset=["SEM"],
    )
    print("\n--- Track B-1: Nested CV summary (unbiased) ---")
    print(nest_summary.to_string(index=False))

    # ============================================================
    # TRACK B 增强: Feature stability matrix
    # ============================================================
    print("\n" + "=" * 78)
    print("[Track B-2] Feature stability matrix across outer folds")
    print("=" * 78)
    stab = compute_feature_stability(
        nest_folds,
        all_candidate_features=candidate_features,
        final_panel=final_panel,
        model_name="SEM",
        output_csv=os.path.join(OUTDIR, "feature_stability_matrix.csv"),
    )
    print("\n--- Feature stability matrix (top 25 by selection_count) ---")
    print(stab.head(25).to_string(index=False))

    # 量化：final panel 里多少 features 是 "highly stable" (≥4/5 folds)
    fp_stab = stab.loc[stab["in_final_panel"] == "Yes"]
    n_stable = (fp_stab["selection_count"] >= 4).sum()
    print(f"\n>>> Final panel stability: "
          f"{n_stable}/{len(final_panel)} features were selected in ≥4/5 outer folds")
    print(f">>> Median selection_count of final panel features: "
          f"{fp_stab['selection_count'].median():.1f} / 5")

    # ============================================================
    # TRACK B-3: Permutation test (fold-internal RFE)
    # ============================================================
    N_PERM = 100   # 正式投稿改成 1000
    print("\n" + "=" * 78)
    print(f"[Track B-3] Permutation test (fold-internal RFE, n={N_PERM})")
    print("=" * 78)
    perm_summary, perm_dist = permutation_test_cv(
        X_train, y_train,
        models=sem_models,
        fixed_features=None,                  # ✅ permutation 也用 fold-internal
        n_features=N_FEATURES, n_splits=5,
        n_permutations=N_PERM,
        random_state=42,
        output_csv=os.path.join(OUTDIR, "permutation_test_SEM.csv"),
    )
    print("\n--- Track B-3: Permutation test summary ---")
    print(perm_summary.to_string(index=False))

    # ============================================================
    # 总览：response Letter 用的数字
    # ============================================================
    print("\n" + "=" * 78)
    print("SUMMARY for response letter")
    print("=" * 78)
    print(f"  Final panel (paper):       {len(final_panel)} features (saved to fixed_panel_20features.csv)")
    print(f"  Track A — Repeated CV:     AUC = {rep_summary.iloc[0]['Mean_AUC']:.3f} "
          f"[{rep_summary.iloc[0]['CI95_low']:.3f}, {rep_summary.iloc[0]['CI95_high']:.3f}]  "
          f"(panel-locked, paper-aligned)")
    print(f"  Track B — Nested CV:       AUC = {nest_summary.iloc[0]['Nested_Mean_AUC']:.3f} "
          f"[{nest_summary.iloc[0]['Nested_CI95_low']:.3f}, "
          f"{nest_summary.iloc[0]['Nested_CI95_high']:.3f}]  "
          f"(fold-internal RFE, leakage-free)")
    print(f"  Track B — Permutation:     p = {perm_summary.iloc[0]['p_value']:.4f} "
          f"(observed AUC = {perm_summary.iloc[0]['True_CV_AUC']:.3f}, "
          f"null mean = {perm_summary.iloc[0]['Null_Mean']:.3f})")
    print(f"  Track B — Panel stability: {n_stable}/{len(final_panel)} features stable across ≥4/5 folds")

    print(f"\nAll outputs in: {OUTDIR}/")


if __name__ == "__main__":
    main()
