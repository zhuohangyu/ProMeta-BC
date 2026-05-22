import pandas as pd
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator




auc_file = "all_models_auc_N1.csv"
all_results = pd.read_csv(auc_file)


# 统一模型顺序（按你实际模型来改）
model_order = ["GLM", "LR", "SVM", "GBM", "XGB", "MLP", "SEM"]
all_results = all_results[all_results["Model"].isin(model_order)].copy()

# 按模型收集 AUC 列表，和你原来 dir_method_cv 的结构类似
dir_model_auc = {}
for m in model_order:
    vals = all_results.loc[all_results["Model"] == m, "AUC"].values
    if len(vals) > 0:
        dir_model_auc[m] = np.array(vals, dtype=np.float64)

# 确保顺序
models = list(dir_model_auc.keys())

# =========================
# 2. 全局画图参数，按你给的风格
# =========================
plt.rcParams['font.sans-serif'] = ['Arial']     # 字体
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
plt.rcParams['figure.figsize'] = (5.0, 5.0)     # 这里我缩到 5x5 inch，更适合文章
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

font = {
    'family': 'Arial',
    'size': 8
}

# =========================
# 3. 小提琴 + 箱线图
# =========================
plt.figure()

# 和你之前一致的 8 色列表，用前 7 个
colors = ["#a032cb","#efbed6","#fc496a", "#52c465",
          "#3372e0", "#7f6cd7","#c4ceff", "#fc9d40"]

data_for_plot = [dir_model_auc[m] for m in models]

# 小提琴图
violin = plt.violinplot(
    data_for_plot,
    widths=0.7,
    showextrema=False
)

# 箱线图
bplot = plt.boxplot(
    data_for_plot,
    widths=0.15,
    patch_artist=True,
    whis=10,
    flierprops={
        'marker': '.',
        'markeredgecolor': '#353b48',
        'markersize': 1
    },
    medianprops={
        'linestyle': '-',
        'color': '#353b48',
        'linewidth': 0.4
    },
    whiskerprops={
        'linestyle': '-',
        'color': '#353b48',
        'linewidth': 0.4
    },
    capprops={
        'linestyle': '-',
        'color': '#353b48',
        'linewidth': 0.4
    }
)

# 给箱线图上色
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_edgecolor('#353b48')
    patch.set_linewidth(0.4)

# 给小提琴上色
for body, color in zip(violin['bodies'], colors):
    body.set_facecolor(color)
    body.set_edgecolor('black')
    body.set_linewidth(0.2)
    body.set_alpha(0.5)


for i, (m, color) in enumerate(zip(models, colors)):
    mean_val = np.mean(dir_model_auc[m])
    plt.scatter(
        i + 1,
        mean_val,
        s=32,                 # 稍微比箱线点大一点
        facecolor=color,     # ⭐ 与模型颜色一致
        edgecolor="#353b48", # 深色描边，保证清晰
        linewidth=0.8,
        zorder=10
    )


# X 轴
x = range(1, len(models) + 1)
plt.xticks(x, models, fontsize=10)


# Y 轴范围按你的要求：0.75–1.0
plt.ylim(0.7, 0.9)


# Y 轴刻度可以加个主刻度间隔（可选）
ax = plt.gca()
y_major_locator = MultipleLocator(0.05)  # 每 0.05 一个刻度
ax.yaxis.set_major_locator(y_major_locator)

plt.xlabel("Machine learning methods", fontsize=12)
plt.ylabel("AUC", fontsize=12)
plt.tick_params(which='major', length=3, labelsize=10, pad=3)

plt.tight_layout()
plt.savefig("FigureS12A.png", transparent=True)
plt.close('all')

