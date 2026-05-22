#!/usr/bin/env python3
"""
plot_validation_results.py (simplified, multi-model)
====================================================

对 validation_results/ 下的输出做可视化（多模型版本）。

只画两张图：
  Fig 1. Repeated CV boxplot (50 fold AUC distribution per model)
  Fig 2. Repeated CV vs Permutation null — overlay per model

（与 plot_validation_results_SEM.py 的区别：这个支持多模型；
 如果只用 SEM，直接用 SEM 版即可。）
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats as sstats

OUTDIR = "./validation_results"
os.makedirs(OUTDIR, exist_ok=True)

plt.rcParams["font.size"] = 11
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"
sns.set_style("whitegrid")

COL_NULL = "#cccccc"
COL_OBS  = "#c0392b"


# ============================================================
# Fig 1. Repeated CV boxplot (per-model)
# ============================================================
folds = pd.read_csv(os.path.join(OUTDIR, "repeated_kfold_cv_folds.csv"))
order = (folds.groupby("Model")["AUC"].mean().sort_values(ascending=False).index.tolist())

fig, ax = plt.subplots(figsize=(8, 5))
sns.boxplot(data=folds, x="Model", y="AUC", order=order, palette="Set2",
            width=0.55, ax=ax, fliersize=0)
sns.stripplot(data=folds, x="Model", y="AUC", order=order,
              color="black", alpha=0.35, size=3, jitter=0.18, ax=ax)
ax.set_ylim(0.6, 1.02)
ax.set_xlabel("Model")
ax.set_ylabel("AUC (50 folds: 5-fold × 10 repeats)")
ax.set_title("Repeated Stratified K-Fold CV on training set")
ax.axhline(0.5, color="gray", linestyle="--", linewidth=1, alpha=0.7)
plt.tight_layout()
out1 = os.path.join(OUTDIR, "fig1_repeated_cv_boxplot.png")
plt.savefig(out1, dpi=300, bbox_inches="tight")
plt.close()
print(f"✅ {out1}")


# ============================================================
# Fig 2. Repeated CV vs Permutation null (per-model overlay)
# ============================================================
perm         = pd.read_csv(os.path.join(OUTDIR, "permutation_test_null_dist.csv"))
perm_summary = pd.read_csv(os.path.join(OUTDIR, "permutation_test.csv"))
rep          = pd.read_csv(os.path.join(OUTDIR, "repeated_kfold_cv.csv"))

models = perm_summary.sort_values("True_CV_AUC", ascending=False)["Model"].tolist()
common = [m for m in models if m in set(rep["Model"])]
n_models = len(common)

fig, axes = plt.subplots(1, n_models, figsize=(5.5 * n_models, 4.6), sharey=False)
if n_models == 1:
    axes = [axes]

for ax, m in zip(axes, common):
    null = perm.loc[perm["Model"] == m, "Perm_AUC"].values
    folds_m = folds.loc[folds["Model"] == m, "AUC"].values
    rep_mean = rep.loc[rep["Model"] == m, "Mean_AUC"].values[0]
    null_mean = perm_summary.loc[perm_summary["Model"] == m, "Null_Mean"].values[0]
    null_std  = perm_summary.loc[perm_summary["Model"] == m, "Null_Std"].values[0]
    pval      = perm_summary.loc[perm_summary["Model"] == m, "p_value"].values[0]

    # KDE null
    kn = sstats.gaussian_kde(null)
    xs1 = np.linspace(null.min() - 0.02, null.max() + 0.02, 400)
    ax.fill_between(xs1, 0, kn(xs1), color=COL_NULL, alpha=0.7,
                    label=f"Permutation null (n={len(null)})")
    # KDE repeated CV
    kr = sstats.gaussian_kde(folds_m)
    xs2 = np.linspace(folds_m.min() - 0.02, folds_m.max() + 0.02, 400)
    ax.fill_between(xs2, 0, kr(xs2), color="#2E86AB", alpha=0.55,
                    label=f"Repeated CV (n={len(folds_m)} folds)")

    ax.axvline(null_mean, color="#34495e", linestyle="--", linewidth=1.3, alpha=0.85)
    ax.axvline(rep_mean, color=COL_OBS, linewidth=2.5)

    gap_sd = (rep_mean - null_mean) / null_std if null_std > 0 else np.nan
    p_text = "p < 0.001" if pval < 0.001 else f"p = {pval:.4f}"
    ax.set_title(f"{m}  —  Δ = {gap_sd:.1f} SDs   |   {p_text}", fontsize=12)
    ax.set_xlabel("AUC")
    ax.set_xlim(0.35, 1.02)
    ax.legend(loc="upper center", fontsize=9, frameon=True)

axes[0].set_ylabel("Density")
fig.suptitle("Repeated CV performance vs Permutation null distribution",
             fontweight="bold", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
out2 = os.path.join(OUTDIR, "fig2_comparison_repCV_vs_permNull.png")
plt.savefig(out2, dpi=300, bbox_inches="tight")
plt.close()
print(f"✅ {out2}")

print(f"\nAll figures saved to {OUTDIR}/")
