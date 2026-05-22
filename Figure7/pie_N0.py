import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch
import matplotlib.patches as mpatches
import matplotlib.cm as cm

# 设置更好的字体
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

# Data
sizes = [27, 5]
labels = ['Non-LNM', 'LNM']
# 使用更现代的配色 - 珊瑚红和天蓝色

# 从 colormap 中取颜色（数值 0~1，越大越深）
n0_color   = cm.get_cmap("Purples")(0.65)   # Non-LNM
npos_color = cm.get_cmap("Oranges")(0.65)   # LNM

colors = [n0_color, npos_color]
explode = (0.1, 0.02)  # 更明显的分离效果

# 创建图形
fig, ax = plt.subplots(figsize=(15, 12))

# 绘制饼图，使用更精致的样式
wedges, texts, autotexts = ax.pie(
    sizes,
    labels=None,  # 我们将自定义标签
    autopct='',  # 我们将自定义百分比显示
    startangle=90,
    colors=colors,
    explode=explode,
    shadow=False,  # 取消默认阴影，我们将添加自定义阴影
    textprops={'fontsize': 40, 'weight': 'bold'},
    wedgeprops={
        'linewidth': 5,
        'antialiased': True,
        'alpha': 0.9,
        'edgecolor': 'white'  # 白色边框更现代
    }
)

# 边框更柔和
wedges[0].set_edgecolor(cm.get_cmap("Purples")(0.8))
wedges[1].set_edgecolor(cm.get_cmap("Oranges")(0.8))

# 百分比
total = sum(sizes)
percentages = [s / total * 100 for s in sizes]

# ===== 标签放在扇形内部（关键）=====
# 你可以调 r_label：越大越靠外，越小越靠中心
r_label = 0.62

for i, (wedge, size, pct) in enumerate(zip(wedges, sizes, percentages)):
    ang = (wedge.theta2 + wedge.theta1) / 2.0
    x = np.cos(np.deg2rad(ang)) * r_label
    y = np.sin(np.deg2rad(ang)) * r_label

    ax.text(
        x, y,
        f"{size}\n({pct:.1f}%)",
        ha='center', va='center',
        fontsize=40, fontweight='bold',
        color='white',
        bbox=dict(
            boxstyle='round,pad=0.45',
            facecolor=colors[i],
            edgecolor='white',
            linewidth=3,
            alpha=0.98
        ),
        zorder=10
    )

# 装饰圆环（可留可删）
for radius, alpha in [(1.10, 0.12), (1.20, 0.06)]:
    circle = Circle((0, 0), radius,
                    fill=False,
                    edgecolor='#BDC3C7',
                    linewidth=2,
                    linestyle='--',
                    alpha=alpha,
                    zorder=0)
    ax.add_patch(circle)


# 图例（精简版）
legend_elements = [
    mpatches.Patch(facecolor=colors[0], edgecolor='white',
                   linewidth=3, label='Non-LNM'),
    mpatches.Patch(facecolor=colors[1], edgecolor='white',
                   linewidth=3, label='LNM')
]

legend = ax.legend(handles=legend_elements,
                   loc='lower right',
                   fontsize=38,
                   frameon=True,
                   fancybox=True,
                   shadow=True,
                   framealpha=0.95,
                   edgecolor='#BDC3C7')
legend.get_title().set_weight('bold')
legend.get_title().set_color('#2C3E50')

ax.axis('equal')
plt.tight_layout()
plt.subplots_adjust(top=0.88, bottom=0.1)
fig.patch.set_alpha(0)
ax.set_facecolor('none')

# 保存
plt.savefig('pie_N0.png',pad_inches=0.6,
            dpi=300, bbox_inches='tight',
            transparent=True,   edgecolor='none')
plt.show()