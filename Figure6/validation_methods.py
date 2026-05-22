#!/usr/bin/env python3
"""
validation_methods.py
=====================

为 BreastCancerAnalysis (process_benign.py) 补充三种统计学严格的验证方法：

1) repeated_kfold_cv     : Repeated Stratified K-Fold CV -> 减少 fold 分割方差，给均值±95%CI
2) nested_cv_evaluation  : Nested CV -> 内层调超参 + 外层评估，给无偏泛化估计
3) permutation_test_cv   : Permutation test -> 通过打乱 y 构造 null 分布，给 p 值

设计原则
--------
- RFE / Scaler 全部放在 fold 内部，避免特征选择 leakage（与你的 select_model_by_cv_on_training 一致）。
- 输入: X (DataFrame) + y (Series 或 array)，可以无缝放到现有 pipeline 的 self.X_train / self.y_train。
- 输出: long-format CSV，方便后续 ggplot / seaborn。

依赖：numpy, pandas, scikit-learn, (可选) xgboost
"""
import os
import time
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier, StackingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score
from sklearn.base import clone

try:
    from xgboost import XGBClassifier
    HAS_XGB = True
except ImportError:
    HAS_XGB = False


# ============================================================
# 候选模型集合（与你 evaluate_models 的超参保持一致）
# ============================================================
def build_default_models(random_state=40, include_sem=True):
    """返回 dict[str -> estimator]。"""
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

    models = {"GLM": glm, "LR": lr, "SVM": svm, "GBM": gbm, "MLP": mlp}

    if HAS_XGB:
        models["XGB"] = XGBClassifier(
            n_estimators=300, max_depth=3, learning_rate=0.05,
            subsample=0.8, colsample_bytree=0.8,
            objective="binary:logistic", eval_metric="logloss",
            random_state=random_state, n_jobs=1, verbosity=0,  # n_jobs=1 for Windows compat
        )

    if include_sem:
        base = [("glm", clone(glm)), ("lr", clone(lr)), ("svm", clone(svm)),
                ("mlp", clone(mlp)), ("gbm", clone(gbm))]
        sem = StackingClassifier(
            estimators=base,
            final_estimator=LogisticRegression(
                penalty="l2", solver="liblinear", max_iter=3000, random_state=random_state),
            passthrough=False, n_jobs=1,
        )
        models["SEM"] = sem
    return models


# ============================================================
# fold 内部 pipeline: 一切 fit_on_train_fold_only
# ============================================================
def _fit_and_score_one_fold(X_tr, X_va, y_tr, y_va, models,
                            n_features=20, fixed_features=None,
                            rfe_step=10, rfe_random_state=40):
    """
    单 fold 内：(可选 SVM-RFE) -> StandardScaler -> 训练所有模型 -> AUC on val
    返回 dict[model_name -> auc]

    fixed_features : list[str] | None
        - None: fold-internal RFE 选 n_features 个特征（避免 selection leakage）
        - list: 直接用给定的固定 panel，跳过 RFE
                  → 与 paper 中"先选 panel 再 evaluate"的工作流一致
    """
    if fixed_features is not None:
        # ===== 用固定 panel（paper 工作流） =====
        if not isinstance(X_tr, pd.DataFrame):
            raise ValueError("使用 fixed_features 时 X 必须是 DataFrame（要按列名取子集）")
        missing = [f for f in fixed_features if f not in X_tr.columns]
        if missing:
            raise ValueError(f"fixed_features 中以下列不在 X 里: {missing[:5]}")
        X_tr_sel = X_tr[fixed_features].values
        X_va_sel = X_va[fixed_features].values
    else:
        # ===== fold-internal RFE =====
        scaler_rfe = StandardScaler()
        X_tr_for_rfe = scaler_rfe.fit_transform(X_tr)
        n_sel = min(n_features, X_tr.shape[1])
        rfe = RFE(estimator=SVC(kernel="linear", random_state=rfe_random_state),
                  n_features_to_select=n_sel, step=rfe_step)
        rfe.fit(X_tr_for_rfe, y_tr)
        mask = rfe.support_
        if isinstance(X_tr, pd.DataFrame):
            X_tr_sel = X_tr.iloc[:, mask].values
            X_va_sel = X_va.iloc[:, mask].values
        else:
            X_tr_sel = X_tr[:, mask]
            X_va_sel = X_va[:, mask]

    # Fold 内重新 standardize（fit on train fold only），这一步两种模式都必须做
    scaler = StandardScaler()
    X_tr2 = scaler.fit_transform(X_tr_sel)
    X_va2 = scaler.transform(X_va_sel)

    out = {}
    for name, clf in models.items():
        try:
            clf_c = clone(clf)
            clf_c.fit(X_tr2, y_tr)
            proba = clf_c.predict_proba(X_va2)[:, 1]
            out[name] = float(roc_auc_score(y_va, proba))
        except Exception:
            out[name] = np.nan
    return out


def _binarize_y(y):
    """统一 y 为 0/1 数组。"""
    if isinstance(y, pd.Series):
        if y.dtype == object or pd.api.types.is_string_dtype(y):
            return (y == "Malignant").values.astype(int)
        return y.values.astype(int)
    arr = np.asarray(y)
    if arr.dtype.kind in ("U", "O", "S"):
        return (arr == "Malignant").astype(int)
    return arr.astype(int)


# ============================================================
# 1) Repeated K-Fold CV
# ============================================================
def repeated_kfold_cv(X, y, models=None, n_features=20, fixed_features=None,
                      n_splits=5, n_repeats=10, random_state=42,
                      output_csv=None, verbose=True):
    """
    Repeated Stratified K-Fold CV
    -----------------------------
    n_splits × n_repeats 个 fold AUC，每次 repeat 用不同 random seed。

    fixed_features : list[str] | None
        给定就用这个固定 panel（跳过 fold-internal RFE），与 paper 工作流一致。

    Returns
    -------
    summary_df, folds_df
    """
    if models is None:
        models = build_default_models()
    y_bin = _binarize_y(y)

    rkf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats,
                                  random_state=random_state)

    t0 = time.time()
    fold_rows = []
    for k, (tr_idx, va_idx) in enumerate(rkf.split(X, y_bin), start=1):
        repeat_idx = (k - 1) // n_splits + 1
        fold_idx = (k - 1) % n_splits + 1

        X_tr = X.iloc[tr_idx] if isinstance(X, pd.DataFrame) else X[tr_idx]
        X_va = X.iloc[va_idx] if isinstance(X, pd.DataFrame) else X[va_idx]
        y_tr, y_va = y_bin[tr_idx], y_bin[va_idx]

        if verbose and fold_idx == 1:
            print(f"  [Repeated CV] repeat {repeat_idx}/{n_repeats} "
                  f"(elapsed {time.time() - t0:.1f}s)")

        fold_aucs = _fit_and_score_one_fold(X_tr, X_va, y_tr, y_va, models,
                                            n_features=n_features,
                                            fixed_features=fixed_features)
        for name, auc in fold_aucs.items():
            fold_rows.append({"Model": name, "Repeat": repeat_idx,
                              "Fold": fold_idx, "AUC": auc})

    folds_df = pd.DataFrame(fold_rows)

    rows = []
    for name in models.keys():
        s = folds_df.loc[folds_df["Model"] == name, "AUC"].dropna().values
        if len(s) == 0:
            continue
        rows.append({
            "Model": name,
            "Mean_AUC": float(np.mean(s)),
            "Std_AUC": float(np.std(s, ddof=1)),
            "CI95_low": float(np.percentile(s, 2.5)),
            "CI95_high": float(np.percentile(s, 97.5)),
            "n_folds": int(len(s)),
        })
    summary_df = pd.DataFrame(rows).sort_values("Mean_AUC", ascending=False).reset_index(drop=True)

    if output_csv:
        os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
        summary_df.to_csv(output_csv, index=False)
        folds_df.to_csv(output_csv.replace(".csv", "_folds.csv"), index=False)
        if verbose:
            print(f"✅ Repeated CV saved: {output_csv}")
    if verbose:
        print(f"✅ Repeated CV done in {time.time() - t0:.1f}s")
    return summary_df, folds_df


# ============================================================
# 2) Nested CV (with per-model hyperparameter search)
# ============================================================
def nested_cv_evaluation(X, y, n_features=20, fixed_features=None,
                         outer_splits=5, inner_splits=3,
                         random_state=42, output_csv=None, verbose=True,
                         models_subset=None):
    """
    Nested CV:
      外层 fold 评估，内层 fold 选超参（每个模型有自己的小网格）

    fixed_features : list[str] | None
        给定就用固定 panel（跳过 fold-internal RFE），与 paper 工作流一致。
    models_subset : list[str] | None
        如 ["SEM"] 只跑 SEM；None = 跑所有。
    """
    y_bin = _binarize_y(y)

    # 轻量级超参网格（保持 demo 可跑；正式可加大）
    grids = {
        "LR": (LogisticRegression(penalty="l2", solver="liblinear",
                                  max_iter=3000, random_state=random_state),
               {"C": [0.1, 1.0, 10.0]}),
        "SVM": (SVC(kernel="linear", probability=True, random_state=random_state),
                {"C": [0.1, 1.0, 10.0]}),
        "GBM": (GradientBoostingClassifier(subsample=0.8, random_state=random_state),
                {"n_estimators": [100, 200],
                 "learning_rate": [0.01, 0.05],
                 "max_depth": [1, 2]}),
        "MLP": (MLPClassifier(activation="relu", solver="adam", max_iter=2000,
                              random_state=random_state),
                {"alpha": [1.0, 10.0, 50.0],
                 "hidden_layer_sizes": [(32,), (64, 32)]}),
    }
    if HAS_XGB:
        grids["XGB"] = (XGBClassifier(
                            objective="binary:logistic", eval_metric="logloss",
                            random_state=random_state, n_jobs=1, verbosity=0),
                        {"n_estimators": [100, 300],
                         "max_depth": [2, 3],
                         "learning_rate": [0.05, 0.1]})

    # === SEM (Stacking) ===
    # SEM 本身没有"标准超参"，唯一有意义的可调点是 final_estimator 的 C
    glm_b = LogisticRegression(penalty="elasticnet", l1_ratio=0.5, solver="saga",
                               max_iter=3000, random_state=random_state)
    lr_b = LogisticRegression(penalty="l2", solver="liblinear",
                              max_iter=3000, random_state=random_state)
    svm_b = SVC(kernel="linear", probability=True, random_state=random_state)
    gbm_b = GradientBoostingClassifier(learning_rate=0.01, n_estimators=200,
                                       max_depth=1, subsample=0.8, random_state=random_state)
    mlp_b = MLPClassifier(learning_rate_init=0.001, hidden_layer_sizes=(64, 32),
                          activation="relu", alpha=50, solver="adam",
                          max_iter=2000, random_state=random_state)
    sem_est = StackingClassifier(
        estimators=[("glm", glm_b), ("lr", lr_b), ("svm", svm_b),
                    ("mlp", mlp_b), ("gbm", gbm_b)],
        final_estimator=LogisticRegression(penalty="l2", solver="liblinear",
                                           max_iter=3000, random_state=random_state),
        passthrough=False, n_jobs=1,
    )
    grids["SEM"] = (sem_est, {"final_estimator__C": [0.1, 1.0, 10.0]})

    if models_subset is not None:
        grids = {k: v for k, v in grids.items() if k in models_subset}
        if not grids:
            raise ValueError(f"models_subset={models_subset} 不在可选列表里")

    outer = StratifiedKFold(n_splits=outer_splits, shuffle=True, random_state=random_state)
    inner = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=random_state + 1)

    t0 = time.time()
    rows = []
    for model_name, (est, grid) in grids.items():
        if verbose:
            print(f"  [Nested CV] model={model_name} (elapsed {time.time() - t0:.1f}s)")
        for k, (outer_tr_idx, outer_te_idx) in enumerate(outer.split(X, y_bin), start=1):
            X_otr = X.iloc[outer_tr_idx] if isinstance(X, pd.DataFrame) else X[outer_tr_idx]
            X_ote = X.iloc[outer_te_idx] if isinstance(X, pd.DataFrame) else X[outer_te_idx]
            y_otr, y_ote = y_bin[outer_tr_idx], y_bin[outer_te_idx]

            # === 特征选择：固定 panel 或 fold-internal RFE ===
            if fixed_features is not None:
                if not isinstance(X_otr, pd.DataFrame):
                    raise ValueError("fixed_features 模式下 X 必须是 DataFrame")
                X_otr_sel = X_otr[fixed_features].values
                X_ote_sel = X_ote[fixed_features].values
                fold_selected = list(fixed_features)  # 每个 fold 都一样
            else:
                scaler_rfe = StandardScaler()
                X_otr_for_rfe = scaler_rfe.fit_transform(X_otr)
                n_sel = min(n_features, X_otr.shape[1])
                rfe = RFE(estimator=SVC(kernel="linear", random_state=random_state),
                          n_features_to_select=n_sel, step=10)
                rfe.fit(X_otr_for_rfe, y_otr)
                mask = rfe.support_
                if isinstance(X_otr, pd.DataFrame):
                    X_otr_sel = X_otr.iloc[:, mask].values
                    X_ote_sel = X_ote.iloc[:, mask].values
                    fold_selected = X_otr.columns[mask].tolist()
                else:
                    X_otr_sel = X_otr[:, mask]
                    X_ote_sel = X_ote[:, mask]
                    fold_selected = [f"feature_{i}" for i, m in enumerate(mask) if m]

            scaler = StandardScaler()
            X_otr2 = scaler.fit_transform(X_otr_sel)
            X_ote2 = scaler.transform(X_ote_sel)

            # 内层 GridSearch (n_jobs=1: 避免 Windows joblib loky 兼容性问题)
            gs = GridSearchCV(clone(est), grid, scoring="roc_auc",
                              cv=inner, n_jobs=1, refit=True)
            gs.fit(X_otr2, y_otr)

            # 外层评估
            proba = gs.predict_proba(X_ote2)[:, 1]
            outer_auc = float(roc_auc_score(y_ote, proba))
            rows.append({
                "Model": model_name,
                "Outer_Fold": k,
                "Outer_AUC": outer_auc,
                "Inner_Best_AUC": float(gs.best_score_),
                "Best_Params": str(gs.best_params_),
                "Selected_Features": ";".join(fold_selected),
                "n_selected": len(fold_selected),
            })

    folds_df = pd.DataFrame(rows)
    summary_rows = []
    for name in grids.keys():
        s = folds_df.loc[folds_df["Model"] == name, "Outer_AUC"].values
        summary_rows.append({
            "Model": name,
            "Nested_Mean_AUC": float(np.mean(s)),
            "Nested_Std_AUC": float(np.std(s, ddof=1)),
            "Nested_CI95_low": float(np.percentile(s, 2.5)),
            "Nested_CI95_high": float(np.percentile(s, 97.5)),
            "n_outer_folds": int(len(s)),
        })
    summary_df = pd.DataFrame(summary_rows).sort_values(
        "Nested_Mean_AUC", ascending=False).reset_index(drop=True)

    if output_csv:
        os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
        summary_df.to_csv(output_csv, index=False)
        folds_df.to_csv(output_csv.replace(".csv", "_folds.csv"), index=False)
        if verbose:
            print(f"✅ Nested CV saved: {output_csv}")
    if verbose:
        print(f"✅ Nested CV done in {time.time() - t0:.1f}s")
    return summary_df, folds_df


# ============================================================
# 3) Permutation test
# ============================================================
def permutation_test_cv(X, y, models=None, n_features=20, fixed_features=None,
                        n_splits=5, n_permutations=200,
                        random_state=42, output_csv=None, verbose=True):
    """
    Permutation test
    ----------------
    对每个候选模型：
      true_score = CV(real y)
      perm_scores[i] = CV(shuffled y), i=1..N
      p = (#{perm >= true} + 1) / (N + 1)
    所有模型共用同一组 permutation（paired），结果可比。

    fixed_features : list[str] | None
        给定就用固定 panel（跳过 fold-internal RFE），与 paper 工作流一致。
    """
    if models is None:
        models = build_default_models()
    y_bin = _binarize_y(y)
    rng = np.random.RandomState(random_state)

    t0 = time.time()

    def _cv_score_for_y(y_use, seed):
        """用给定 y 跑一次 stratified k-fold，返回 dict[model -> mean fold AUC]"""
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        per_model = {m: [] for m in models.keys()}
        for tr_idx, va_idx in skf.split(X, y_use):
            X_tr = X.iloc[tr_idx] if isinstance(X, pd.DataFrame) else X[tr_idx]
            X_va = X.iloc[va_idx] if isinstance(X, pd.DataFrame) else X[va_idx]
            y_tr, y_va = y_use[tr_idx], y_use[va_idx]
            fold_aucs = _fit_and_score_one_fold(X_tr, X_va, y_tr, y_va, models,
                                                n_features=n_features,
                                                fixed_features=fixed_features)
            for m, a in fold_aucs.items():
                per_model[m].append(a)
        return {m: float(np.nanmean(v)) for m, v in per_model.items()}

    # (1) 真实分数
    if verbose:
        print(f"  [Permutation] computing true CV scores ...")
    true_scores = _cv_score_for_y(y_bin, seed=random_state)

    # (2) Null
    perm_rows = []
    for p in range(n_permutations):
        y_perm = rng.permutation(y_bin)
        perm_scores = _cv_score_for_y(y_perm, seed=random_state)
        for m, s in perm_scores.items():
            perm_rows.append({"Model": m, "Perm_idx": p + 1, "Perm_AUC": s})
        if verbose and (p + 1) % max(1, n_permutations // 10) == 0:
            print(f"    perm {p + 1}/{n_permutations} "
                  f"(elapsed {time.time() - t0:.1f}s)")

    perm_df = pd.DataFrame(perm_rows)

    # (3) p 值
    rows = []
    for m in models.keys():
        null = perm_df.loc[perm_df["Model"] == m, "Perm_AUC"].values
        true_s = true_scores[m]
        pval = (np.sum(null >= true_s) + 1) / (n_permutations + 1)
        rows.append({
            "Model": m,
            "True_CV_AUC": float(true_s),
            "Null_Mean": float(np.mean(null)),
            "Null_Std": float(np.std(null, ddof=1)),
            "p_value": float(pval),
            "n_permutations": int(n_permutations),
        })
    summary_df = pd.DataFrame(rows).sort_values("True_CV_AUC", ascending=False).reset_index(drop=True)

    if output_csv:
        os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
        summary_df.to_csv(output_csv, index=False)
        perm_df.to_csv(output_csv.replace(".csv", "_null_dist.csv"), index=False)
        if verbose:
            print(f"✅ Permutation test saved: {output_csv}")
    if verbose:
        print(f"✅ Permutation test done in {time.time() - t0:.1f}s")
    return summary_df, perm_df


# ============================================================
# Feature stability matrix from nested CV
# ============================================================
def compute_feature_stability(folds_df, all_candidate_features,
                              final_panel=None, model_name=None,
                              output_csv=None):
    """
    根据 nested_cv_evaluation 的 folds_df 计算 feature stability matrix。

    Parameters
    ----------
    folds_df : pd.DataFrame
        nested_cv_evaluation 返回的第二个 DataFrame，需要包含 'Selected_Features' 列。
    all_candidate_features : list[str]
        所有候选特征（你的 27+7 = 34 个），用来确保未被任何 fold 选中的特征也出现在表里。
    final_panel : list[str] | None
        你 paper 报告的 final panel（在整个 training set 上 RFE 选出的）。
        给定后会在输出表加一列 'in_final_panel'。
    model_name : str | None
        如果 folds_df 含多个 model（一般 nested CV 是按 model 跑的），
        只统计这个 model；None = 用第一个 model。
    output_csv : str | None
        保存路径。

    Returns
    -------
    stability_df : pd.DataFrame
        每行一个候选特征，列含每个 outer fold 是否被选中 (0/1)、selection_count、in_final_panel
    """
    if "Selected_Features" not in folds_df.columns:
        raise ValueError("folds_df 缺少 'Selected_Features' 列，"
                         "说明 nested_cv_evaluation 是在 fixed_features 模式下跑的，无 stability 可算")

    if model_name is None:
        model_name = folds_df["Model"].iloc[0]
    sub = folds_df.loc[folds_df["Model"] == model_name].sort_values("Outer_Fold")

    n_folds = len(sub)
    # 矩阵 (n_features × n_folds)
    matrix = {}
    for _, row in sub.iterrows():
        fold_idx = int(row["Outer_Fold"])
        selected = set(row["Selected_Features"].split(";")) if row["Selected_Features"] else set()
        for f in all_candidate_features:
            matrix.setdefault(f, {})[f"fold_{fold_idx}"] = int(f in selected)

    stab = pd.DataFrame.from_dict(matrix, orient="index").fillna(0).astype(int)
    stab.index.name = "feature"
    fold_cols = [c for c in stab.columns if c.startswith("fold_")]
    stab["selection_count"] = stab[fold_cols].sum(axis=1)
    stab["selection_freq"] = stab["selection_count"] / n_folds

    if final_panel is not None:
        stab["in_final_panel"] = np.where(
            stab.index.isin(set(final_panel)), "Yes", "No"
        )

    stab = stab.sort_values(["selection_count", "feature"], ascending=[False, True])
    stab = stab.reset_index()

    if output_csv:
        os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
        stab.to_csv(output_csv, index=False)
        print(f"✅ Feature stability matrix saved: {output_csv}")

    return stab
