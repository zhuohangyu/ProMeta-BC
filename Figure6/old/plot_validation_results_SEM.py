#!/usr/bin/env python3
"""
plot_validation_results_SEM.py (v2 - polished)
==============================================
对 validation_results_SEM/ 下的输出做可视化。

改进点：
- Fig 4 stability heatmap: 高稳定性特征在顶部，去除冗余 0/1 文本，
  bar 按 stability tier 着色，增加 panel 摘要标注
- Fig 1: 加 violin 叠加 + mean/median 标注
- Fig 2: 加 delta AUC 标注
- Fig 3: 加 KDE 曲线 + null 分布统计 + gap 标注
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
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
COL_A      = "#2E86AB"
COL_B      = "#E07A5F"
COL_OBS    = "#C0392B"
COL_NULL   = "#BDC3C7"
COL_PANEL  = "#C0392B"


# ============================================================
# Fig 1. Repeated CV (Track A)
# ============================================================
folds = pd.read_csv(os.path.join(OUTDIR, "repeated_kfold_cv_SEM_folds.csv"))
summary_A = pd.read_csv(os.path.join(OUTDIR, "repeated_kfold_cv_SEM.csv"))

mean_a = summary_A.iloc[0]["Mean_AUC"]
std_a  = summary_A.iloc[0]["Std_AUC"]
lo_a, hi_a = summary_A.iloc[0]["CI95_low"], summary_A.iloc[0]["CI95_high"]

fig, ax = plt.subplots(figsize=(5, 6))
parts = ax.violinplot([folds["AUC"].values], positions=[0], widths=0.7,
                       showmeans=False, showmedians=False, showextrema=False)
for body in parts["bodies"]:
    body.set_facecolor(COL_A)
    body.set_edgecolor(COL_A)
    body.set_alpha(0.25)

ax.boxplot([folds["AUC"].values], positions=[0], widths=0.18,
            patch_artist=True, showfliers=False,
            medianprops=dict(color="black", linewidth=2),
            boxprops=dict(facecolor=COL_A, edgecolor=COL_A, linewidth=1.5, alpha=0.85),
            whiskerprops=dict(color=COL_A, linewidth=1.5),
            capprops=dict(color=COL_A, linewidth=1.5))

ax.scatter(np.random.RandomState(0).normal(0, 0.04, len(folds)),
           folds["AUC"].values, s=14, color="black", alpha=0.45, zorder=3)
ax.scatter(0, mean_a, marker="D", s=80, color=COL_OBS, zorder=4,
           edgecolor="white", linewidth=1.5, label=f"Mean = {mean_a:.3f}")

# ax.axhspan(lo_a, hi_a, color=COL_OBS, alpha=0.08, zorder=0)
ax.set_xticks([0])
ax.set_xticklabels(["SEM\n(panel-locked)"], fontsize=12)
ax.set_ylim(0.80, 1.01)
ax.set_ylabel("AUC", fontsize=12)
# ax.set_title("Track A — Repeated CV\n(5-fold × 10 repeats = 50 folds)",
#              fontsize=13)
stats_txt = (f"Mean = {mean_a:.3f}\n"
             f"SD = {std_a:.3f}\n"
             f"95% CI = [{lo_a:.3f}, {hi_a:.3f}]\n"
             f"N = {len(folds)} folds")
ax.text(0.55, 0.05, stats_txt, transform=ax.transAxes, fontsize=10,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                  edgecolor=COL_A, alpha=0.95))
ax.grid(axis="y", alpha=0.3)
ax.grid(axis="x", visible=False)
plt.tight_layout()
out1 = os.path.join(OUTDIR, "fig1_repeated_cv_TrackA.png")
plt.savefig(out1, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ {out1}")


# ============================================================
# Fig 2. Track A vs Track B forest plot
# ============================================================
summary_B = pd.read_csv(os.path.join(OUTDIR, "nested_cv_SEM.csv"))
mean_b = summary_B.iloc[0]["Nested_Mean_AUC"]
lo_b, hi_b = summary_B.iloc[0]["Nested_CI95_low"], summary_B.iloc[0]["Nested_CI95_high"]
delta = mean_a - mean_b

fig, ax = plt.subplots(figsize=(9, 3.8))
ax.errorbar(mean_a, 1, xerr=[[mean_a - lo_a], [hi_a - mean_a]],
            fmt="o", color=COL_A, capsize=7, markersize=13, linewidth=2.8,
            markeredgecolor="white", markeredgewidth=1.5,
            label=f"Track A · panel-locked Repeated CV: {mean_a:.3f}  [{lo_a:.3f}, {hi_a:.3f}]")
ax.errorbar(mean_b, 0, xerr=[[mean_b - lo_b], [hi_b - mean_b]],
            fmt="s", color=COL_B, capsize=7, markersize=13, linewidth=2.8,
            markeredgecolor="white", markeredgewidth=1.5,
            label=f"Track B · fold-internal RFE Nested CV: {mean_b:.3f}  [{lo_b:.3f}, {hi_b:.3f}]")

# ax.annotate("", xy=(mean_b, 0.55), xytext=(mean_a, 0.45),
#             arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
ax.text((mean_a + mean_b) / 2, 0.5, f"  Δ = {delta:+.3f}",
        ha="left", va="center", fontsize=11, fontweight="bold")

ax.axvline(0.5, color="gray", linestyle=":", linewidth=1.2, alpha=0.7)
ax.set_yticks([0, 1])
ax.set_yticklabels(["Nested 5-fold\ncross-validation", "Repeated 5-fold\ncross-validation"],
                    fontsize=11)
ax.set_xlabel("AUC (with 95% CI)", fontsize=12)
ax.set_xlim(min(0.7, lo_b - 0.05), 1.0)
ax.set_ylim(-0.6, 1.7)
# ax.set_title("Track A vs Track B — performance is concordant\n"
#              "→ feature selection leakage does not meaningfully inflate AUC",
#              fontsize=12.5)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.28), frameon=True,
          fontsize=10, framealpha=0.95, ncol=1)
ax.grid(axis="x", alpha=0.3)
ax.grid(axis="y", visible=False)
plt.tight_layout()
out2 = os.path.join(OUTDIR, "fig2_TrackA_vs_TrackB.png")
plt.savefig(out2, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ {out2}")


# ============================================================
# Fig 3. Permutation null distribution + KDE
# ============================================================
perm = pd.read_csv(os.path.join(OUTDIR, "permutation_test_SEM_null_dist.csv"))
perm_summary = pd.read_csv(os.path.join(OUTDIR, "permutation_test_SEM.csv"))

null = perm["Perm_AUC"].values
true_auc  = perm_summary.iloc[0]["True_CV_AUC"]
null_mean = perm_summary.iloc[0]["Null_Mean"]
null_std  = perm_summary.iloc[0]["Null_Std"]
pval      = perm_summary.iloc[0]["p_value"]
n_perm    = int(perm_summary.iloc[0]["n_permutations"])

fig, ax = plt.subplots(figsize=(9, 4.8))
counts_h, bins, _ = ax.hist(null, bins=30, color=COL_NULL,
                              edgecolor="white", alpha=0.85,
                              label=f"Null distribution (n = {n_perm} permutations)")
# KDE
kde = sstats.gaussian_kde(null)
xs = np.linspace(null.min() - 0.02, null.max() + 0.02, 200)
ys = kde(xs) * len(null) * (bins[1] - bins[0])
ax.plot(xs, ys, color="#7f8c8d", linewidth=2, alpha=0.9)

# Observed
ax.axvline(true_auc, color=COL_OBS, linewidth=3,
           label=f"Observed AUC = {true_auc:.3f}")
# Null mean (放在更靠下不挡 legend)
ax.axvline(null_mean, color="#34495e", linestyle="--", linewidth=1.2, alpha=0.7,
           label=f"Null mean = {null_mean:.3f}")
ax.axvline(0.5, color="gray", linestyle=":", linewidth=1, alpha=0.5)

# 注释 gap (放在中下区域)
y_arrow = ax.get_ylim()[1] * 0.35
ax.annotate("", xy=(true_auc - 0.005, y_arrow), xytext=(null_mean + 0.005, y_arrow),
            arrowprops=dict(arrowstyle="->", color=COL_OBS, lw=1.4))
ax.text((null_mean + true_auc) / 2, y_arrow * 1.18,
        f"Δ = {true_auc - null_mean:.3f}  ({(true_auc - null_mean)/null_std:.1f} SDs)",
        ha="center", color=COL_OBS, fontsize=10.5, fontweight="bold")

ax.set_xlim(0.35, 1.0)
ax.set_xlabel("CV AUC under shuffled labels", fontsize=12)
ax.set_ylabel("Permutation count", fontsize=12)
title = (f"Permutation test with fold-internal RFE  —  "
         f"p {'< 0.001' if pval < 0.001 else f'= {pval:.4f}'}")
ax.set_title(title, fontsize=13)
ax.legend(loc="upper center", fontsize=10, frameon=True, framealpha=0.95)
ax.grid(axis="y", alpha=0.3)
ax.grid(axis="x", alpha=0.2)
plt.tight_layout()
out3 = os.path.join(OUTDIR, "fig3_permutation_null.png")
plt.savefig(out3, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ {out3}")


# ============================================================
# Fig 4. Feature stability heatmap — 改良版
# ============================================================
stab = pd.read_csv(os.path.join(OUTDIR, "feature_stability_matrix.csv"))

fold_cols = [c for c in stab.columns if c.startswith("fold_")]
n_folds = len(fold_cols)

# 排序：highest 在顶部，同 count 内 panel 特征优先
stab_sorted = stab.sort_values(
    ["selection_count", "in_final_panel", "feature"],
    ascending=[False, False, True]
).reset_index(drop=True)

mat = stab_sorted[fold_cols].values
feature_names = stab_sorted["feature"].values
in_panel = stab_sorted["in_final_panel"].values
counts = stab_sorted["selection_count"].values
n_features = len(feature_names)

# Stability tier 颜色
tier_color = {
    5: "#3F007D",
    4: "#6A51A3",
    3: "#3182BD",
    2: "#FDAE6B",
    1: "#F16913",
    0: "#BDBDBD",
}
bar_colors = [tier_color[c] for c in counts]

fig_height = max(7, 0.32 * n_features + 1.5)
fig = plt.figure(figsize=(11, fig_height))
gs = fig.add_gridspec(1, 2, width_ratios=[3.4, 1.4], wspace=0.05)
ax_hm = fig.add_subplot(gs[0, 0])
ax_bar = fig.add_subplot(gs[0, 1], sharey=ax_hm)

# Heatmap - 不要 invert y axis，highest count 自动在 top
ax_hm.imshow(mat, aspect="auto", cmap="Blues", vmin=0, vmax=1,
              interpolation="nearest")
for i in range(n_features + 1):
    ax_hm.axhline(i - 0.5, color="white", linewidth=1.2)
for j in range(n_folds + 1):
    ax_hm.axvline(j - 0.5, color="white", linewidth=1.2)

ax_hm.set_xticks(np.arange(n_folds))
ax_hm.set_xticklabels([c.replace("fold_", "Fold ") for c in fold_cols], fontsize=11)
ax_hm.set_yticks(np.arange(n_features))
text_objs = ax_hm.set_yticklabels(feature_names, fontsize=10)
for txt, p in zip(text_objs, in_panel):
    if p == "Yes":
        txt.set_color(COL_PANEL)
        txt.set_fontweight("bold")
ax_hm.tick_params(axis="y", length=0)
ax_hm.tick_params(axis="x", length=4)
for spine in ax_hm.spines.values():
    spine.set_visible(False)
ax_hm.set_title("Selection across outer folds\n(red bold = in final ProMeta-BC panel)",
                 fontsize=12.5, pad=12)

# Bar
ax_bar.barh(np.arange(n_features), counts, color=bar_colors,
             edgecolor="white", linewidth=1.2, height=0.78)
for i, p in enumerate(in_panel):
    if p == "Yes":
        ax_bar.text(counts[i] + 0.15, i, "★", color=COL_PANEL,
                    fontsize=12, va="center", ha="left", fontweight="bold")

ax_bar.set_xlim(0, n_folds + 0.6)
ax_bar.set_xticks(range(n_folds + 1))
ax_bar.set_yticks(np.arange(n_features))
ax_bar.tick_params(axis="y", labelleft=False, length=0)
ax_bar.set_xlabel("Selection count (out of 5 folds)", fontsize=11)
ax_bar.set_title("Stability", fontsize=12.5, pad=12)
ax_bar.axvline(4, color="#7f8c8d", linestyle="--", linewidth=1.3, alpha=0.8)
ax_bar.grid(axis="x", alpha=0.3)
ax_bar.spines["left"].set_visible(False)

# 图例
legend_elements = [
    Patch(facecolor=tier_color[5], label="5/5"),
    Patch(facecolor=tier_color[4], label="4/5"),
    Patch(facecolor=tier_color[3], label="3/5"),
    Patch(facecolor=tier_color[2], label="2/5"),
    Patch(facecolor=tier_color[1], label="1/5"),
    Patch(facecolor=tier_color[0], label="0/5"),
    Line2D([0], [0], marker="*", color="w", markerfacecolor=COL_PANEL,
           markersize=14, label="in final panel"),
]
ax_bar.legend(handles=legend_elements, loc="upper center",
               bbox_to_anchor=(0.5, -0.06), ncol=4, frameon=True,
               fontsize=9, framealpha=0.95)

# 顶部 summary (放在 x label 下方避免和 subplot title 冲突)
fp = stab[stab["in_final_panel"] == "Yes"]
n_stable = int((fp["selection_count"] >= 4).sum())
median_count = int(fp["selection_count"].median())
# summary_txt = (f"Final panel stability:  {n_stable}/{len(fp)} features selected in "
#                f"≥4/5 outer folds   |   median selection count = {median_count}/5")

plt.tight_layout()
# 在 tight_layout 之后用 figtext 放摘要，避开 subplot title
# fig.text(0.5, 0.99, summary_txt, fontsize=12.5, fontweight="bold",
#          color="#2c3e50", ha="center", va="top",
#          bbox=dict(boxstyle="round,pad=0.6", facecolor="#ecf0f1",
#                    edgecolor="#bdc3c7", alpha=0.9))
plt.subplots_adjust(top=0.93)
out4 = os.path.join(OUTDIR, "fig4_feature_stability_heatmap.png")
plt.savefig(out4, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ {out4}")

print(f"\n>>> Final panel summary:")
print(f"    {n_stable}/{len(fp)} features in ≥4/5 outer folds")
print(f"    median selection count = {median_count}/5")
