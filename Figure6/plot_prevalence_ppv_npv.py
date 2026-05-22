#!/usr/bin/env python3
"""
plot_prevalence_ppv_npv.py
==========================
基于 Cohort 3 (external validation) 的 sensitivity/specificity (含 95% CI)，
依据 Bayes' theorem 绘制 PPV / NPV 随 prevalence 变化的曲线。

设计：
- 曲线 + 95%CI band
- 编号 marker (1, 2, 3, 4A, 4B, 4C, 5) 在曲线上，每个 marker 独立配色
- Scenario legend 在右侧第三个 panel 里, 跟曲线视觉对齐, 不重叠

BI-RADS 4 按 ACR 习惯拆为 4A / 4B / 4C 三档:
  4A (low suspicion):    10%
  4B (moderate):         20%
  4C / pre-op ALN:       40%

输出:
  SupplementaryFig_prevalence_PPV_NPV.png
  SupplementaryTable_prevalence_PPV_NPV.csv
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


# ============================================================
# 输入参数 — 改这里就行
# ============================================================
SENS       = 0.953271
SENS_LO    = 0.908144
SENS_HI    = 0.990478
SPEC       = 0.870370
SPEC_LO    = 0.769231
SPEC_HI    = 0.945455

COHORT3_PREV = 107 / 161

OUT_DIR  = "."
FIG_NAME = "SupplementaryFig_prevalence_PPV_NPV.png"
TBL_NAME = "SupplementaryTable_prevalence_PPV_NPV.csv"


# ============================================================
# Bayes
# ============================================================
def ppv(prev, sens, spec):
    return (sens * prev) / (sens * prev + (1 - spec) * (1 - prev))


def npv(prev, sens, spec):
    return (spec * (1 - prev)) / ((1 - sens) * prev + spec * (1 - prev))


# ============================================================
# 网格
# ============================================================
prev = np.concatenate([
    np.linspace(0.001, 0.01, 60),
    np.linspace(0.01, 0.70, 250),
])
ppv_mean = ppv(prev, SENS, SPEC)
npv_mean = npv(prev, SENS, SPEC)
ppv_lo = ppv(prev, SENS_LO, SPEC_LO)
ppv_hi = ppv(prev, SENS_HI, SPEC_HI)
npv_lo = npv(prev, SENS_LO, SPEC_LO)
npv_hi = npv(prev, SENS_HI, SPEC_HI)


# ============================================================
# 关键场景 (编号 marker + 配色)
# BI-RADS 4 拆成 4A / 4B / 4C，颜色用渐变（黄 → 橙 → 深橙）
# ============================================================
scenarios = [
    {"id": "1",  "prev": 0.001,        "label": "Population screening",        "color": "#7F8C8D"},
    {"id": "2",  "prev": 0.005,        "label": "High-risk asymptomatic",      "color": "#3498DB"},
    {"id": "3",  "prev": 0.05,         "label": "BI-RADS 3 / mild clinical",   "color": "#27AE60"},
    {"id": "4A", "prev": 0.10,         "label": "BI-RADS 4A (low suspicion)",  "color": "#F4D03F"},
    {"id": "4B", "prev": 0.20,         "label": "BI-RADS 4B (moderate)",       "color": "#F39C12"},
    {"id": "4C", "prev": 0.40,         "label": "BI-RADS 4C / pre-op ALN",     "color": "#E67E22"},
    {"id": "5",  "prev": COHORT3_PREV, "label": "Cohort 3 (observed)",         "color": "#C0392B"},
]


# ============================================================
# 绘图: 用 GridSpec 做 3 列布局 (2 个曲线 + 1 个 legend panel)
# ============================================================
fig = plt.figure(figsize=(16.5, 5.8))
gs  = GridSpec(1, 3, width_ratios=[1.0, 1.0, 1.05], wspace=0.32,
               left=0.05, right=0.99, top=0.86, bottom=0.13)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax_leg = fig.add_subplot(gs[0, 2])


def draw_marker(ax, x, y, color, num, size=320):
    """画一个带白边 + 编号 + 配色 marker。
    对于多字符 id (4A/4B/4C)，字号略缩小以保持在圆内。"""
    ax.scatter([x], [y], color=color, s=size, zorder=5,
               edgecolor="white", linewidth=2.2)
    fs = 10.5 if len(num) > 1 else 12
    ax.text(x, y, num, color="white", fontsize=fs, fontweight="bold",
            ha="center", va="center", zorder=6)


# --- (A) PPV ---
ax1.fill_between(prev * 100, ppv_lo * 100, ppv_hi * 100,
                  color="#3498DB", alpha=0.15)
ax1.plot(prev * 100, ppv_mean * 100, color="#1A5490", linewidth=2.8, zorder=3)

ax1.text(0.03, 0.95,
          f"Sens = {SENS:.3f}\nSpec = {SPEC:.3f}\n(Cohort 3)",
          transform=ax1.transAxes, fontsize=9.5, va="top", ha="left",
          bbox=dict(boxstyle="round,pad=0.45", facecolor="#EBF5FB",
                    edgecolor="#3498DB", linewidth=1.2),
          family="monospace")

for s in scenarios:
    y = ppv(s["prev"], SENS, SPEC) * 100
    draw_marker(ax1, s["prev"] * 100, y, s["color"], s["id"])

ax1.set_xscale("log")
ax1.set_xlim(0.08, 80)
ax1.set_ylim(0, 100)
ax1.set_xlabel("Disease prevalence (%, log scale)", fontsize=11, fontweight="bold")
ax1.set_ylabel("Positive Predictive Value, PPV (%)", fontsize=11, fontweight="bold")
ax1.set_title("(A) PPV vs prevalence", fontsize=12.5, fontweight="bold", loc="left")
ax1.grid(True, alpha=0.3, which="both")
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# --- (B) NPV ---
ax2.fill_between(prev * 100, npv_lo * 100, npv_hi * 100,
                  color="#27AE60", alpha=0.15)
ax2.plot(prev * 100, npv_mean * 100, color="#196F3D", linewidth=2.8, zorder=3)

ax2.text(0.03, 0.30,
          f"Sens = {SENS:.3f}\nSpec = {SPEC:.3f}\n(Cohort 3)",
          transform=ax2.transAxes, fontsize=9.5, va="top", ha="left",
          bbox=dict(boxstyle="round,pad=0.45", facecolor="#E8F8F0",
                    edgecolor="#27AE60", linewidth=1.2),
          family="monospace")

for s in scenarios:
    y = npv(s["prev"], SENS, SPEC) * 100
    draw_marker(ax2, s["prev"] * 100, y, s["color"], s["id"])

ax2.set_xscale("log")
ax2.set_xlim(0.08, 80)
ax2.set_ylim(85, 100.4)
ax2.set_xlabel("Disease prevalence (%, log scale)", fontsize=11, fontweight="bold")
ax2.set_ylabel("Negative Predictive Value, NPV (%)", fontsize=11, fontweight="bold")
ax2.set_title("(B) NPV vs prevalence", fontsize=12.5, fontweight="bold", loc="left")
ax2.grid(True, alpha=0.3, which="both")
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)


# --- (C) Legend Panel ---
ax_leg.axis("off")
ax_leg.set_xlim(0, 1)
ax_leg.set_ylim(0, 1)
ax_leg.set_title("Scenario lookup", fontsize=12.5, fontweight="bold", loc="left")

# 表头
y_header = 0.94
ax_leg.text(0.06, y_header, "#",         fontsize=10.5, fontweight="bold", ha="center")
ax_leg.text(0.14, y_header, "Scenario",  fontsize=10.5, fontweight="bold", ha="left")
ax_leg.text(0.66, y_header, "Prev.",     fontsize=10.5, fontweight="bold", ha="center")
ax_leg.text(0.80, y_header, "PPV",       fontsize=10.5, fontweight="bold", ha="center")
ax_leg.text(0.93, y_header, "NPV",       fontsize=10.5, fontweight="bold", ha="center")

# 表头分割线
ax_leg.plot([0.02, 0.99], [y_header - 0.035, y_header - 0.035],
             color="#34495E", linewidth=1.2, transform=ax_leg.transAxes)

# 每个 scenario 一行 — 7 行，行距略收紧
row_height = 0.094
y_start = 0.83
for i, s in enumerate(scenarios):
    y = y_start - i * row_height
    ppv_val = ppv(s["prev"], SENS, SPEC) * 100
    npv_val = npv(s["prev"], SENS, SPEC) * 100
    prev_pct = s["prev"] * 100

    # 配色圆形 + 编号
    ax_leg.scatter([0.06], [y], color=s["color"], s=260, zorder=5,
                    edgecolor="white", linewidth=2.0,
                    transform=ax_leg.transAxes)
    fs_marker = 9.5 if len(s["id"]) > 1 else 10.5
    ax_leg.text(0.06, y, s["id"], color="white", fontsize=fs_marker, fontweight="bold",
                 ha="center", va="center", zorder=6,
                 transform=ax_leg.transAxes)

    # Scenario label
    ax_leg.text(0.14, y, s["label"], fontsize=9.8, va="center", ha="left",
                 transform=ax_leg.transAxes)

    # Prev / PPV / NPV 数值
    prev_str = f"{prev_pct:.1f}%" if prev_pct >= 1 else f"{prev_pct:.2f}%"
    ax_leg.text(0.66, y, prev_str, fontsize=9.8, va="center", ha="center",
                 transform=ax_leg.transAxes, family="monospace")
    ax_leg.text(0.80, y, f"{ppv_val:.1f}%", fontsize=9.8, va="center", ha="center",
                 transform=ax_leg.transAxes, family="monospace",
                 color="#1A5490", fontweight="bold")
    ax_leg.text(0.93, y, f"{npv_val:.2f}%", fontsize=9.8, va="center", ha="center",
                 transform=ax_leg.transAxes, family="monospace",
                 color="#196F3D", fontweight="bold")

# 底部注释 — 关于外推的诚实声明
note = ("Markers 1–2 represent prevalences below the validation cohort's empirical "
        "range; PPV and NPV at these prevalences are theoretical projections under "
        "the assumption of population-invariant sensitivity and specificity. "
        "Shaded bands depict 95% CI propagated from sens/spec CIs.")
import textwrap
wrapped = "\n".join(textwrap.wrap(note, width=58))
ax_leg.text(0.0, 0.13, wrapped, fontsize=8.2, va="top", ha="left",
             transform=ax_leg.transAxes, style="italic", color="#566573")


# 大标题
fig.suptitle(
    "Supplementary Figure SXX — Prevalence-dependent PPV and NPV of ProMeta-BC "
    "(Bayes' theorem, Cohort 3 sens/spec)",
    fontsize=12.5, fontweight="bold", y=0.97)

fig_path = os.path.join(OUT_DIR, FIG_NAME)
fig.savefig(fig_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()


# ============================================================
# 数值表 — CSV 里也按 BI-RADS 4A/4B/4C 明确标注
# ============================================================
prev_table = [0.001, 0.005, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, COHORT3_PREV]
labels_table = [
    "Population screening (0.1%)",
    "Asymptomatic high-risk (0.5%)",
    "General population subset (1%)",
    "BI-RADS 3 / mild clinical (5%)",
    "BI-RADS 4A — low suspicion (10%)",
    "BI-RADS 4B — moderate (20%)",
    "BI-RADS 4B / 4C boundary (30%)",
    "BI-RADS 4C / pre-op ALN (40%)",
    "Clinically suspicious (50%)",
    f"Cohort 3 observed ({COHORT3_PREV*100:.1f}%)",
]
rows = []
for p, lbl in zip(prev_table, labels_table):
    rows.append({
        "Scenario":         lbl,
        "Prevalence_%":     round(p * 100, 2),
        "PPV_%":            round(ppv(p, SENS, SPEC) * 100, 2),
        "PPV_CI_low_%":     round(ppv(p, SENS_LO, SPEC_LO) * 100, 2),
        "PPV_CI_high_%":    round(ppv(p, SENS_HI, SPEC_HI) * 100, 2),
        "NPV_%":            round(npv(p, SENS, SPEC) * 100, 2),
        "NPV_CI_low_%":     round(npv(p, SENS_LO, SPEC_LO) * 100, 2),
        "NPV_CI_high_%":    round(npv(p, SENS_HI, SPEC_HI) * 100, 2),
    })
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT_DIR, TBL_NAME), index=False)

print(f"✅ {fig_path}")
print(f"✅ {os.path.join(OUT_DIR, TBL_NAME)}")
