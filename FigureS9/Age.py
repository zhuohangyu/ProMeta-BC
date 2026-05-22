import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap


def get_soft_cmap(cmap_name, minval=0.15, maxval=0.85, n=256):
    """
    裁剪 colormap，避免最深色过暗
    """
    base = cm.get_cmap(cmap_name)
    new_colors = base(np.linspace(minval, maxval, n))
    return LinearSegmentedColormap.from_list(
        f"{cmap_name}_soft", new_colors
    )

CSV_PATH = "AgeSubgroup_metrics_combined.csv"
OUT_DIR_SINGLE = "confusion_by_row"
OUT_PNG_GRID = "ConfusionMatrix_AgeSubgroups_grid.png"

os.makedirs(OUT_DIR_SINGLE, exist_ok=True)


def _sanitize_filename(s: str) -> str:
    s = str(s)
    s = re.sub(r"[^\w\-\.]+", "_", s)
    return s.strip("_")


def get_cmap_by_age_subset(subset: str):
    s = str(subset).lower()

    if any(k in s for k in ["young", "lt40", "<40", "<=40"]):
        return get_soft_cmap("Oranges", 0.15, 0.80)

    if any(k in s for k in ["old", "ge55", ">=55", "elder"]):
        return get_soft_cmap("Purples", 0.15, 0.80)

    return get_soft_cmap("Greens", 0.15, 0.80)



def plot_cm_from_counts(
    TN: int, FP: int, FN: int, TP: int,
    title: str,
    out_png: str,
    cmap="Purples",  # 👈 新增
    display=("Benign", "Maligant"),
    annot_size=22,
    label_size=18,
    title_size=20
):
    """
    用计数(TN,FP,FN,TP)绘制2x2混淆矩阵.
    矩阵布局:
        [[TN, FP],
         [FN, TP]]
    行=Actual(Non-LN, LN), 列=Predicted(Non-LN, LN)
    """
    cm = np.array([[TN, FP],
                   [FN, TP]], dtype=int)

    fig, ax = plt.subplots(figsize=(4.8, 4.3))
    im = ax.imshow(cm, cmap=cmap)   # 👈 使用 cmap

    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(display, fontsize=label_size)
    ax.set_yticklabels(display, fontsize=label_size)

    ax.set_xlabel("Predicted", fontsize=label_size, fontweight="bold")
    ax.set_ylabel("Actual", fontsize=label_size, fontweight="bold")
    # ax.set_title(title, fontsize=title_size)

    # 自动字体颜色阈值（深色背景 → 白字）
    threshold = cm.max() * 0.5

    for (i, j), v in np.ndenumerate(cm):
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
    """
    网格图：行=Subset，列=Cohort。
    若有多个Model，会为每个Model分别输出一张图（文件名自动加model后缀）。
    """
    required = {model_col, subset_col, cohort_col, tn_col, fp_col, fn_col, tp_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV缺少必要列: {sorted(missing)}")

    def _natkey(x):
        # 让 Cohort 2 < Cohort 10 这种自然排序
        return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(x))]

    for model, dmodel in df.groupby(model_col, dropna=False):
        subsets = sorted(dmodel[subset_col].dropna().unique())
        cohorts = sorted(dmodel[cohort_col].dropna().unique(), key=_natkey)

        nrows, ncols = len(subsets), len(cohorts)
        fig, axes = plt.subplots(nrows, ncols, figsize=(4.8*ncols, 4.5*nrows), squeeze=False)

        vmax = int(np.max(dmodel[[tn_col, fp_col, fn_col, tp_col]].values))

        last_im = None
        for i, subset in enumerate(subsets):
            for j, cohort in enumerate(cohorts):
                ax = axes[i, j]
                row = dmodel[(dmodel[subset_col] == subset) & (dmodel[cohort_col] == cohort)]
                if row.empty:
                    ax.axis("off")
                    continue
                row = row.iloc[0]

                cm = np.array([[int(row[tn_col]), int(row[fp_col])],
                               [int(row[fn_col]), int(row[tp_col])]], dtype=int)

                last_im = ax.imshow(cm, vmin=0, vmax=vmax)

                ax.set_xticks([0, 1])
                ax.set_yticks([0, 1])
                ax.set_xticklabels(display, fontsize=tick_size)
                ax.set_yticklabels(display, fontsize=tick_size)

                if i == nrows - 1:
                    ax.set_xlabel("Predicted", fontsize=tick_size)
                if j == 0:
                    ax.set_ylabel("Actual", fontsize=tick_size)

                ax.set_title(f"{subset} | {cohort}", fontsize=title_size, fontweight="bold")

                for (ii, jj), v in np.ndenumerate(cm):
                    ax.text(jj, ii, f"{v:d}", ha="center", va="center", fontsize=annot_size)

        if last_im is not None:
            cbar = fig.colorbar(last_im, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02)
            cbar.ax.tick_params(labelsize=tick_size)
            cbar.set_label("Freq", fontsize=tick_size)

        fig.suptitle(f"Confusion Matrices ({model})", fontsize=title_size+2, fontweight="bold", y=1.02)

        # 多个model则自动加后缀
        out_path = out_png
        if df[model_col].nunique() > 1:
            root, ext = os.path.splitext(out_png)
            out_path = f"{root}_{_sanitize_filename(model)}{ext}"

        fig.tight_layout()
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved grid: {out_path}")


if __name__ == "__main__":
    df = pd.read_csv(CSV_PATH)

    # ---- 1) 每行输出单张混淆矩阵 ----
    for _, r in df.iterrows():
        model = r["Model"]
        subset = r["Subset"]
        cohort = r["Cohort"]
        title = f"{model} | {subset} | {cohort}"

        fname = f"CM_{_sanitize_filename(model)}_{_sanitize_filename(subset)}_{_sanitize_filename(cohort)}.png"
        out_png = os.path.join(OUT_DIR_SINGLE, fname)

        cmap = get_cmap_by_age_subset(subset)

        plot_cm_from_counts(
            TN=int(r["TN"]), FP=int(r["FP"]),
            FN=int(r["FN"]), TP=int(r["TP"]),
            title=f"{subset} | {cohort}",
            out_png=out_png,
            cmap=cmap
        )

    print(f"Saved single CMs to: {OUT_DIR_SINGLE}")

    # ---- 2) 网格汇总图 ----
    plot_cm_grid_from_df(df, out_png=OUT_PNG_GRID)
