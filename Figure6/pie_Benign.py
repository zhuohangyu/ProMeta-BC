import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
import matplotlib.patches as mpatches

# =========================
# 全局字体设置
# =========================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

# =========================
# Data
# =========================
sizes = [1, 8]
labels = ['Unpredicted', 'Predicted']   # 只用于 legend
colors = ['#FF6B6B', '#4ECDC4']
explode = (0.10, 0.02)

# =========================
# Figure
# =========================
fig, ax = plt.subplots(figsize=(18, 10), facecolor='#F8F9FA')

wedges, texts = ax.pie(
    sizes,
    labels=None,
    startangle=90,
    colors=colors,
    explode=explode,
    shadow=False,
    wedgeprops={
        'linewidth': 5,
        'antialiased': True,
        'alpha': 0.92,
        'edgecolor': 'white'
    }
)

# 边框
wedges[0].set_edgecolor('#E85A5A')
wedges[1].set_edgecolor('#45B7AF')

# =========================
# 百分比
# =========================
total = sum(sizes)
percentages = [s / total * 100 for s in sizes]

# =========================
# 饼图内部：只显示 数值 + 百分比
# =========================
r_label_big = 0.62
r_label_small = 0.55

for i, (wedge, size, pct) in enumerate(zip(wedges, sizes, percentages)):
    ang = (wedge.theta2 + wedge.theta1) / 2.0
    r_label = r_label_big if size == max(sizes) else r_label_small

    x = np.cos(np.deg2rad(ang)) * r_label
    y = np.sin(np.deg2rad(ang)) * r_label

    ax.text(
        x, y,
        f"{size}\n({pct:.1f}%)",
        ha='center', va='center',
        fontsize=46 if size == max(sizes) else 46,
        fontweight='bold',
        color='white',
        linespacing=1.2,
        bbox=dict(
            boxstyle='round,pad=0.45',
            facecolor=colors[i],
            edgecolor='white',
            linewidth=3,
            alpha=0.98
        ),
        zorder=10
    )


# =========================
# Legend（只在这里出现 Unpredicted / Predicted）
# =========================
legend_elements = [
    mpatches.Patch(facecolor=colors[0], edgecolor='white',
                   linewidth=3, label='Unpredicted'),
    mpatches.Patch(facecolor=colors[1], edgecolor='white',
                   linewidth=3, label='Predicted')
]

ax.legend(
    handles=legend_elements,
    loc='lower left',
    fontsize=38,
    frameon=True,
    fancybox=True,
    framealpha=0.95,
    edgecolor='#BDC3C7'
)

# # =========================
# # Title
# # =========================
# fig.text(
#     0.5, 0.96,
#     'Imaging Benign Surgical Malignant',
#     ha='center',
#     fontsize=26,
#     weight='bold',
#     color='#2C3E50',
#     family='Arial'
# )

# =========================
# Layout & Save
# =========================
ax.axis('equal')
plt.tight_layout()
plt.subplots_adjust(top=0.88, bottom=0.08)

plt.savefig(
    'pie_Benign.png',
    dpi=300,
    bbox_inches='tight',
    facecolor='white',
    edgecolor='none'
)
plt.show()

print("✓ Pie chart saved (labels only in legend, values inside pie)!")
