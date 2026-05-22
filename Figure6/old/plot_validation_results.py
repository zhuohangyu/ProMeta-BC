#!/usr/bin/env python3
"""
plot_validation_results.py
==========================
对 validation_results/ 下三个 CSV 做可视化：
1) Repeated CV: per-model boxplot
2) Nested CV vs Repeated CV: 并排 forest plot 比对
3) Permutation null distribution + 真实 AUC 标注
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

OUTDIR = "./validation_results"
os.makedirs(OUTDIR, exist_ok=True)

plt.rcParams["font.size"] = 11
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"
sns.set_style("whitegrid")


# =====================================================
# Fig 1. Repeated CV boxplot (50 fold AUC distribution per model)
# =====================================================
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
ax.set_title("Repeated Stratified K-Fold CV on Cohort-1 training set")
ax.axhline(0.5, color="gray", linestyle="--", linewidth=1, alpha=0.7)
plt.tight_layout()
out1 = os.path.join(OUTDIR, "fig1_repeated_cv_boxplot.png")
plt.savefig(out1, dpi=300, bbox_inches="tight")
plt.close()
print(f"✅ {out1}")


# =====================================================
# Fig 2. Repeated CV vs Nested CV forest plot
# =====================================================
rep = pd.read_csv(os.path.join(OUTDIR, "repeated_kfold_cv.csv"))
nest = pd.read_csv(os.path.join(OUTDIR, "nested_cv.csv"))

common = sorted(set(rep["Model"]) & set(nest["Model"]),
                key=lambda m: -rep.loc[rep["Model"] == m, "Mean_AUC"].values[0])

fig, ax = plt.subplots(figsize=(8, 5))
y_pos = np.arange(len(common))
offset = 0.18

for i, m in enumerate(common):
    r = rep.loc[rep["Model"] == m].iloc[0]
    n = nest.loc[nest["Model"] == m].iloc[0]

    # Repeated CV (蓝)
    ax.errorbar(r["Mean_AUC"], i - offset,
                xerr=[[r["Mean_AUC"] - r["CI95_low"]],
                      [r["CI95_high"] - r["Mean_AUC"]]],
                fmt="o", color="#2b93db", capsize=4, markersize=8,
                label="Repeated CV (5x10)" if i == 0 else None, linewidth=2)
    # Nested CV (橙)
    ax.errorbar(n["Nested_Mean_AUC"], i + offset,
                xerr=[[n["Nested_Mean_AUC"] - n["Nested_CI95_low"]],
                      [n["Nested_CI95_high"] - n["Nested_Mean_AUC"]]],
                fmt="s", color="#FF6B35", capsize=4, markersize=8,
                label="Nested CV (5 outer × 3 inner)" if i == 0 else None, linewidth=2)

ax.set_yticks(y_pos)
ax.set_yticklabels(common)
ax.set_xlabel("AUC (95% CI)")
ax.set_title("Repeated CV vs Nested CV — point estimate and 95% CI")
ax.axvline(0.5, color="gray", linestyle="--", linewidth=1, alpha=0.6)
ax.set_xlim(0.6, 1.0)
ax.invert_yaxis()
ax.legend(loc="lower left", frameon=True)
plt.tight_layout()
out2 = os.path.join(OUTDIR, "fig2_repeated_vs_nested_forest.png")
plt.savefig(out2, dpi=300, bbox_inches="tight")
plt.close()
print(f"✅ {out2}")


# =====================================================
# Fig 3. Permutation null distribution + true AUC marker
# =====================================================
perm = pd.read_csv(os.path.join(OUTDIR, "permutation_test_null_dist.csv"))
perm_summary = pd.read_csv(os.path.join(OUTDIR, "permutation_test.csv"))

models = perm_summary.sort_values("True_CV_AUC", ascending=False)["Model"].tolist()
n_models = len(models)

fig, axes = plt.subplots(1, n_models, figsize=(4.5 * n_models, 4.2), sharey=True)
if n_models == 1:
    axes = [axes]

for ax, m in zip(axes, models):
    null = perm.loc[perm["Model"] == m, "Perm_AUC"].values
    true_auc = perm_summary.loc[perm_summary["Model"] == m, "True_CV_AUC"].values[0]
    pval = perm_summary.loc[perm_summary["Model"] == m, "p_value"].values[0]

    ax.hist(null, bins=22, color="#cccccc", edgecolor="white", alpha=0.85,
            label=f"Null (n={len(null)})")
    ax.axvline(true_auc, color="#c0392b", linewidth=2.5,
               label=f"Observed AUC = {true_auc:.3f}")
    ax.axvline(0.5, color="gray", linestyle="--", linewidth=1, alpha=0.6)

    ax.set_title(f"{m} — p = {pval:.4f}", fontsize=13)
    ax.set_xlabel("CV AUC")
    ax.set_xlim(0.35, 1.0)
    ax.legend(loc="upper center", fontsize=9, frameon=True)

axes[0].set_ylabel("Permutation count")
fig.suptitle("Permutation test — null distribution of CV AUC under shuffled labels",
             fontweight="bold", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
out3 = os.path.join(OUTDIR, "fig3_permutation_null_distribution.png")
plt.savefig(out3, dpi=300, bbox_inches="tight")
plt.close()
print(f"✅ {out3}")

print("\nAll figures saved to", OUTDIR)
