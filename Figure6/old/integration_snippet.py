"""
integration_snippet.py
======================

把 validation_methods.py 的三个方法接到 BreastCancerAnalysis 里。
你只需要：
  1) 把 validation_methods.py 放到 process_benign.py 同一目录下
  2) 把下面这三个方法 paste 到 BreastCancerAnalysis class 里
  3) 在 run_complete_analysis() 里调用（见末尾示例）
"""

# ===== Step 1. 文件顶部 import =====
from validation_methods import (
    repeated_kfold_cv,
    nested_cv_evaluation,
    permutation_test_cv,
    build_default_models,
)


# ===== Step 2. 加进 BreastCancerAnalysis class =====
# 全部放在 select_model_by_cv_on_training 后面就行

def repeated_cv_on_training(self, n_splits=5, n_repeats=10, n_features=20, random_state=42):
    """
    在 training set 上做 Repeated Stratified K-Fold CV。
    输出 mean ± std ± 95% CI（基于 50 fold 的经验分布）。
    """
    out_csv = f"{self.output_dir}/RepeatedCV_{self.feature_mode}.csv"
    summary, folds = repeated_kfold_cv(
        self.X_train, self.y_train,
        n_features=n_features,
        n_splits=n_splits, n_repeats=n_repeats,
        random_state=random_state,
        output_csv=out_csv,
    )
    print(f"\n=== Repeated CV ({n_splits}×{n_repeats}) — {self.feature_mode} ===")
    print(summary.to_string(index=False))
    self.repeated_cv_summary = summary
    self.repeated_cv_folds = folds
    return summary, folds


def nested_cv_on_training(self, outer_splits=5, inner_splits=3, n_features=20, random_state=42):
    """
    在 training set 上做 Nested CV。
    每个模型有自己的小型超参 grid，内层 grid search + 外层评估。
    """
    out_csv = f"{self.output_dir}/NestedCV_{self.feature_mode}.csv"
    summary, folds = nested_cv_evaluation(
        self.X_train, self.y_train,
        n_features=n_features,
        outer_splits=outer_splits, inner_splits=inner_splits,
        random_state=random_state,
        output_csv=out_csv,
    )
    print(f"\n=== Nested CV (outer={outer_splits}, inner={inner_splits}) — {self.feature_mode} ===")
    print(summary.to_string(index=False))
    self.nested_cv_summary = summary
    self.nested_cv_folds = folds
    return summary, folds


def permutation_test_on_training(self, n_permutations=1000, n_splits=5,
                                  n_features=20, random_state=42,
                                  models_subset=None):
    """
    在 training set 上做 Permutation test，得到每个模型的 p value。
    正式 paper 用 n_permutations=1000；快速 demo 用 100。

    models_subset: list[str], 子集模型名（去掉重模型加速）。
                   None = 全模型。例如 ["LR", "SVM", "GBM", "XGB"] 比较快。
    """
    out_csv = f"{self.output_dir}/PermutationTest_{self.feature_mode}.csv"

    models = build_default_models(include_sem=False)  # SEM permutation 太重
    if models_subset is not None:
        models = {k: v for k, v in models.items() if k in models_subset}

    summary, perm_dist = permutation_test_cv(
        self.X_train, self.y_train,
        models=models,
        n_features=n_features, n_splits=n_splits,
        n_permutations=n_permutations,
        random_state=random_state,
        output_csv=out_csv,
    )
    print(f"\n=== Permutation test (n={n_permutations}) — {self.feature_mode} ===")
    print(summary.to_string(index=False))
    self.permutation_summary = summary
    self.permutation_null = perm_dist
    return summary, perm_dist


# ===== Step 3. 在 run_complete_analysis 里调用 =====
# 把下面这段插在 self.feature_selection(n_features=20) 之后、self.evaluate_models() 之前

"""
# --- 严格的统计学验证（training set only） ---
self.repeated_cv_on_training(n_splits=5, n_repeats=10, n_features=20)
self.nested_cv_on_training(outer_splits=5, inner_splits=3, n_features=20)

# permutation test 比较重，正式 paper 用 1000；调试可用 100
self.permutation_test_on_training(
    n_permutations=1000,
    models_subset=["LR", "SVM", "GBM", "XGB"]   # 跳过 SEM/MLP 太重
)
"""
