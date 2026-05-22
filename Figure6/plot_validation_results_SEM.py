#!/usr/bin/env python3
"""
plot_validation_results_SEM.py (simplified)
===========================================
对 validation_results_SEM/ 下的输出做可视化。

只画两张图：
  Fig 1. Repeated CV violin  (panel-locked AUC distribution, 50 folds)
  Fig 2. Repeated CV vs Permutation null — overlay
         （核心图：真实 AUC 分布 vs null 分布的对比，直观看到 gap）
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats as sstats

OUTDIR = "./validation_results_SEM"
os.makedirs(OUTDIR, exist_ok=True)

# ===== 统一样式 =====
sns.set_style("white")
plt.rcParams.update({
    "font.size": 11,
    "font.family": "DejaVu Sans",
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})

# 配色
COL_A    = "#2E86AB"     # Repeated CV (蓝)
COL_OBS  = "#C0392B"     # Observed / true AUC
COL_NULL = "#BDC3C7"     # Permutation null (灰)


# ============================================================
# 读数据
# ============================================================
folds      = pd.read_csv(os.path.join(OUTDIR, "repeated_kfold_cv_SEM_folds.csv"))
summary_A  = pd.read_csv(os.path.join(OUTDIR, "repeated_kfold_cv_SEM.csv"))
perm       = pd.read_csv(os.path.join(OUTDIR, "permutation_test_SEM_null_dist.csv"))
perm_summ  = pd.read_csv(os.path.join(OUTDIR, "permutation_test_SEM.csv"))

mean_a = summary_A.iloc[0]["Mean_AUC"]
std_a  = summary_A.iloc[0]["Std_AUC"]
lo_a, hi_a = summary_A.iloc[0]["CI95_low"], summary_A.iloc[0]["CI95_high"]

null      = perm["Perm_AUC"].values
true_auc  = perm_summ.iloc[0]["True_CV_AUC"]
null_mean = perm_summ.iloc[0]["Null_Mean"]
null_std  = perm_summ.iloc[0]["Null_Std"]
pval      = perm_summ.iloc[0]["p_value"]
n_perm    = int(perm_summ.iloc[0]["n_permutations"])


# ============================================================
# Fig 1. Repeated CV (panel-locked) violin
# ============================================================
fig, ax = plt.subplots(figsize=(5, 6))

parts = ax.violinplot([folds["AUC"].values], positions=[0], widths=0.7,
                      showmeans=False, showmedians=False, showextrema=False)
for body in parts["bodies"]:
    body.set_facecolor(COL_A); body.set_edgecolor(COL_A); body.set_alpha(0.25)

ax.boxplot([folds["AUC"].values], positions=[0], widths=0.18,
           patch_artist=True, showfliers=False,
           medianprops=dict(color="black", linewidth=2),
           boxprops=dict(facecolor=COL_A, edgecolor=COL_A, linewidth=1.5, alpha=0.85),
           whiskerprops=dict(color=COL_A, linewidth=1.5),
           capprops=dict(color=COL_A, linewidth=1.5))

ax.scatter(np.random.RandomState(0).normal(0, 0.04, len(folds)),
           folds["AUC"].values, s=14, color="black", alpha=0.45, zorder=3)
ax.scatter(0, mean_a, marker="D", s=80, color=COL_OBS, zorder=4,
           edgecolor="white", linewidth=1.5)

ax.set_xticks([0])
ax.set_xticklabels(["SEM\n(panel-locked)"], fontsize=12)
ax.set_ylim(0.80, 1.01)
ax.set_ylabel("AUC", fontsize=12)

stats_txt = (f"Mean = {mean_a:.3f}\n"
             f"SD = {std_a:.3f}\n"
             f"95% CI = [{lo_a:.3f}, {hi_a:.3f}]\n"
             f"N = {len(folds)} folds")
ax.text(0.55, 0.05, stats_txt, transform=ax.transAxes, fontsize=10,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                  edgecolor=COL_A, alpha=0.95))
ax.grid(axis="y", alpha=0.3); ax.grid(axis="x", visible=False)
plt.tight_layout()
out1 = os.path.join(OUTDIR, "fig1_repeated_cv.png")
plt.savefig(out1, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ {out1}")


# ============================================================
# Fig 2. Repeated CV distribution vs Permutation null — overlay
# ============================================================
folds_auc = folds["AUC"].values

fig, ax = plt.subplots(figsize=(10, 5))

# KDE for null distribution
kde_null = sstats.gaussian_kde(null)
xs_null = np.linspace(null.min() - 0.02, null.max() + 0.02, 400)
ax.fill_between(xs_null, 0, kde_null(xs_null), color=COL_NULL, alpha=0.7,
                label=f"Permutation null (n = {n_perm})")
ax.plot(xs_null, kde_null(xs_null), color="#7f8c8d", linewidth=1.6, alpha=0.9)

# KDE for repeated CV distribution
kde_rep = sstats.gaussian_kde(folds_auc)
xs_rep = np.linspace(folds_auc.min() - 0.02, folds_auc.max() + 0.02, 400)
ax.fill_between(xs_rep, 0, kde_rep(xs_rep), color=COL_A, alpha=0.55,
                label=f"Repeated CV (n = {len(folds_auc)} folds)")
ax.plot(xs_rep, kde_rep(xs_rep), color=COL_A, linewidth=1.6, alpha=0.95)

# Mean lines
ax.axvline(null_mean, color="#34495e", linestyle="--", linewidth=1.5, alpha=0.85,
           label=f"Null mean = {null_mean:.3f}")
ax.axvline(mean_a, color=COL_OBS, linewidth=2.6,
           label=f"Repeated CV mean = {mean_a:.3f}")
ax.axvline(0.5, color="gray", linestyle=":", linewidth=1, alpha=0.5)

# Gap annotation
gap = mean_a - null_mean
n_sds_a = gap / null_std
y_max = max(kde_null(xs_null).max(), kde_rep(xs_rep).max())
y_arrow = y_max * 0.55
ax.annotate("", xy=(mean_a - 0.005, y_arrow), xytext=(null_mean + 0.005, y_arrow),
            arrowprops=dict(arrowstyle="->", color=COL_OBS, lw=1.5))
ax.text((null_mean + mean_a) / 2, y_arrow * 1.12,
        f"Δ = {gap:.3f}  ({n_sds_a:.1f} SDs of null)",
        ha="center", color=COL_OBS, fontsize=11, fontweight="bold")

ax.set_xlim(0.35, 1.02)
ax.set_xlabel("AUC", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
title = (f"Repeated CV vs Permutation null  —  "
         f"p {'< 0.001' if pval < 0.001 else f'= {pval:.4f}'}")
ax.set_title(title, fontsize=12.5)
ax.legend(loc="upper center", fontsize=10, frameon=True, framealpha=0.95, ncol=2)
ax.grid(axis="y", alpha=0.3); ax.grid(axis="x", alpha=0.2)
plt.tight_layout()
out2 = os.path.join(OUTDIR, "fig2_comparison_repCV_vs_permNull.png")
plt.savefig(out2, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ {out2}")

print(f"\nAll figures saved to {OUTDIR}/")
