import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap


def get_soft_cmap(cmap_name, minval=0.15, maxval=0.80, n=256):
    """裁剪 colormap，避免最深色过暗"""
    base = cm.get_cmap(cmap_name)
    new_colors = base(np.linspace(minval, maxval, n))
    return LinearSegmentedColormap.from_list(f"{cmap_name}_soft", new_colors)


# =========================
# 1) 文件与输出
# =========================
CSV_PATH = "BMISubgroup_metrics_combined.csv"   # ✅ BMI 文件
OUT_DIR_SINGLE = "confusion_by_row_BMI"
OUT_PNG_GRID = "ConfusionMatrix_BMISubgroups_grid.png"
os.makedirs(OUT_DIR_SINGLE, exist_ok=True)


def _sanitize_filename(s: str) -> str:
    s = str(s)
    s = re.sub(r"[^\w\-\.]+", "_", s)
    return s.strip("_")


# =========================
# 2) BMI -> cmap 规则
# =========================
def get_cmap_by_bmi_subset(subset: str):
    """
    低BMI/偏瘦/正常 -> Oranges
    高BMI/超重/肥胖 -> Purples
    """
    s = str(subset).lower()

    # 低 BMI 相关关键词（你可以按你自己的 Subset 命名再补充）
    low_keys = [
        "low", "lean", "normal", "under", "lt25", "<25", "<=25",
        "bmi<25", "bmi_lt25", "bmi< 25", "bmi <=25"
    ]
    # 高 BMI 相关关键词
    high_keys = [
        "high", "over", "obese", "overweight", "ge25", ">=25",
        "bmi>=25", "bmi_ge25", "bmi>= 25", "bmi >=25"
    ]

    if any(k in s for k in low_keys):
        return get_soft_cmap("Oranges", 0.15, 0.80)

    if any(k in s for k in high_keys):
        return get_soft_cmap("Purples", 0.15, 0.80)

    # 兜底
    return get_soft_cmap("Greens", 0.15, 0.80)


# =========================
# 3) 单张混淆矩阵
# =========================
def plot_cm_from_counts(
    TN: int, FP: int, FN: int, TP: int,
    title: str,
    out_png: str,
    cmap="Purples",
    display=("Benign", "Maligant"),
    annot_size=22,
    label_size=18,
    title_size=20
):
    cmx = np.array([[TN, FP],
                    [FN, TP]], dtype=int)

    fig, ax = plt.subplots(figsize=(4.8, 4.3))
    im = ax.imshow(cmx, cmap=cmap)

    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(display, fontsize=label_size)
    ax.set_yticklabels(display, fontsize=label_size)

    ax.set_xlabel("Predicted", fontsize=label_size, fontweight="bold")
    ax.set_ylabel("Actual", fontsize=label_size, fontweight="bold")
    # ax.set_title(title, fontsize=title_size)

    # --- 自动字体颜色阈值 ---
    threshold = cmx.max() * 0.5  # 深色区域切换白字

    for (i, j), v in np.ndenumerate(cmx):
        color = "white" if v > threshold else "black"

        ax.text(
            j, i, f"{v:d}",
            ha="center", va="center",
            fontsize=annot_size,
            color=color,
            fontweight="bold"
        )

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=label_size)
    cbar.set_label("Freq", fontsize=label_size)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


# =========================
# 4) 网格混淆矩阵（Subset × Cohort）
# =========================
def plot_cm_grid_from_df(
    df: pd.DataFrame,
    out_png: str,
    model_col="Model",
    subset_col="Subset",
    cohort_col="Cohort",
    tn_col="TN", fp_col="FP", fn_col="FN", tp_col="TP",
    display=("Benign", "Maligant"),
    annot_size=16,
    tick_size=12,
    title_size=14
):
    required = {model_col, subset_col, cohort_col, tn_col, fp_col, fn_col, tp_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV缺少必要列: {sorted(missing)}")

    def _natkey(x):
        return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(x))]

    for model, dmodel in df.groupby(model_col, dropna=False):
        subsets = sorted(dmodel[subset_col].dropna().unique())
        cohorts = sorted(dmodel[cohort_col].dropna().unique(), key=_natkey)

        nrows, ncols = len(subsets), len(cohorts)
        fig, axes = plt.subplots(nrows, ncols, figsize=(4.8*ncols, 4.5*nrows), squeeze=False)

        vmax = int(np.max(dmodel[[tn_col, fp_col, fn_col, tp_col]].values))
        last_im = None

        for i, subset in enumerate(subsets):
            cmap_use = get_cmap_by_bmi_subset(subset)  # ✅ 每个 Subset 使用自己的颜色系

            for j, cohort in enumerate(cohorts):
                ax = axes[i, j]
                row = dmodel[(dmodel[subset_col] == subset) & (dmodel[cohort_col] == cohort)]
                if row.empty:
                    ax.axis("off")
                    continue
                row = row.iloc[0]

                cmx = np.array([[int(row[tn_col]), int(row[fp_col])],
                                [int(row[fn_col]), int(row[tp_col])]], dtype=int)

                last_im = ax.imshow(cmx, cmap=cmap_use, vmin=0, vmax=vmax)  # ✅ 用 BMI cmap

                ax.set_xticks([0, 1])
                ax.set_yticks([0, 1])
                ax.set_xticklabels(display, fontsize=tick_size)
                ax.set_yticklabels(display, fontsize=tick_size)

                if i == nrows - 1:
                    ax.set_xlabel("Predicted", fontsize=tick_size)
                if j == 0:
                    ax.set_ylabel("Actual", fontsize=tick_size)

                ax.set_title(f"{subset} | {cohort}", fontsize=title_size, fontweight="bold")

                for (ii, jj), v in np.ndenumerate(cmx):
                    ax.text(jj, ii, f"{v:d}", ha="center", va="center", fontsize=annot_size)

        if last_im is not None:
            cbar = fig.colorbar(last_im, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02)
            cbar.ax.tick_params(labelsize=tick_size)
            cbar.set_label("Freq", fontsize=tick_size)

        fig.suptitle(f"Confusion Matrices ({model})", fontsize=title_size+2, fontweight="bold", y=1.02)

        out_path = out_png
        if df[model_col].nunique() > 1:
            root, ext = os.path.splitext(out_png)
            out_path = f"{root}_{_sanitize_filename(model)}{ext}"

        fig.tight_layout()
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved grid: {out_path}")


# =========================
# 5) 主程序
# =========================
if __name__ == "__main__":
    df = pd.read_csv(CSV_PATH)

    # ---- 1) 每行输出单张混淆矩阵 ----
    for _, r in df.iterrows():
        model = r["Model"]
        subset = r["Subset"]
        cohort = r["Cohort"]

        fname = f"CM_{_sanitize_filename(model)}_{_sanitize_filename(subset)}_{_sanitize_filename(cohort)}.png"
        out_png = os.path.join(OUT_DIR_SINGLE, fname)

        cmap_use = get_cmap_by_bmi_subset(subset)

        plot_cm_from_counts(
            TN=int(r["TN"]), FP=int(r["FP"]),
            FN=int(r["FN"]), TP=int(r["TP"]),
            title=f"{subset} | {cohort}",
            out_png=out_png,
            cmap=cmap_use
        )

    print(f"Saved single CMs to: {OUT_DIR_SINGLE}")

    # ---- 2) 网格汇总图 ----
    plot_cm_grid_from_df(df, out_png=OUT_PNG_GRID)
