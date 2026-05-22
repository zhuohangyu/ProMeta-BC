## ===== 蛋白组学热图可视化 =====
## 模仿代谢组学可视化方式

## ===== 1) 加载包 =====
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(tibble)
library(grid)

## ===== 2) 读入数据 =====
prot <- read_csv("Sample1.csv")
clin <- read_csv("../Figure2/Clinical.csv")

## ===== 3) 处理临床数据：将数字分型转为文字 =====
clin <- clin %>%
  mutate(`Molecular Classification` = case_when(
    Label == "Benign" ~ "Benign",
    `Molecular Classification` == 1 ~ "HR+",
    `Molecular Classification` == 2 ~ "HER2+",
    `Molecular Classification` == 3 ~ "TNBC",
    TRUE ~ NA_character_
  ))

## ===== 4) 筛选蛋白（四个比较都不是 unchanged）=====
prot_filt <- prot %>%
  filter(Regulation_ALL != "unchanged",
         `Regulation_HR+`  != "unchanged",
         `Regulation_HER+` != "unchanged",
         Regulation_TNBC != "unchanged")

cat("筛选后蛋白数量:", nrow(prot_filt), "\n")

## 如果筛选结果太少，可以放宽条件（至少有一个不是unchanged）
if(nrow(prot_filt) < 10) {
  cat("筛选结果太少，放宽为至少3个比较显著变化...\n")
  prot_filt <- prot %>%
    mutate(n_changed = (Regulation_ALL != "unchanged") + 
             (`Regulation_HR+` != "unchanged") +
             (`Regulation_HER+` != "unchanged") +
             (Regulation_TNBC != "unchanged")) %>%
    filter(n_changed >= 3) %>%
    select(-n_changed)
  cat("放宽后蛋白数量:", nrow(prot_filt), "\n")
}

## ===== 5.1) 只保留 log2FC 前15的上、前15的下，用于可视化 =====
# 只要上调 / 下调，去掉 unchanged
prot_top <- prot_filt %>%
  filter(Regulation_ALL %in% c("up-regulated", "down-regulated"))

# 最上调的 15 个（log2FC 最大）
up_top <- prot_top %>%
  filter(Regulation_ALL == "up-regulated") %>%
  arrange(desc(log2FC_ALL)) %>%
  slice_head(n = 25)

# 最下调的 15 个（log2FC 最小，越负越下调）
down_top <- prot_top %>%
  filter(Regulation_ALL == "down-regulated") %>%
  arrange(log2FC_ALL) %>%   # 从小到大，最负在前
  slice_head(n = 25)

# 合并起来，保证先上调后下调（方便行 split）
prot_filt <- bind_rows(down_top, up_top) %>%
  mutate(Regulation_ALL = factor(
    Regulation_ALL,
    levels = c("down-regulated", "up-regulated")   # 下调在前
  )) %>%
  arrange(Regulation_ALL, log2FC_ALL)

cat("用于可视化的蛋白数量:", nrow(prot_filt), "\n")

## ===== 6) 构建表达矩阵（行=蛋白，列=样本）=====
## 找到样本列的起始位置（Regulation_TNBC 之后）
sample_start <- which(colnames(prot_filt) == "Regulation_TNBC") + 1

## 表达矩阵
expr_mat <- prot_filt %>%
  select(`Gene name`, all_of(sample_start:ncol(.))) %>%
  column_to_rownames("Gene name") %>%
  as.matrix()

## 将空字符串转为NA，再转为数值
expr_mat[expr_mat == ""] <- NA
expr_mat <- apply(expr_mat, 2, as.numeric)
rownames(expr_mat) <- prot_filt$`Gene name`

## ===== Z-score 标准化 (行方向) =====
expr_mat <- t(scale(t(expr_mat)))

## ===== 7) 列分组（根据 Clinical.csv）=====
## 确保只保留在表达矩阵中存在的样本
common_samples <- intersect(colnames(expr_mat), clin$Sample)
expr_mat <- expr_mat[, common_samples]

group <- clin$`Molecular Classification`
names(group) <- clin$Sample
group <- group[colnames(expr_mat)]
group <- factor(group, levels = c("Benign", "HR+", "HER2+", "TNBC"))

cat("样本分组统计:\n")
print(table(group))

## ===== 8) 行信息（log2FC, P-value）=====
row_df <- prot_filt %>%
  transmute(
    feature = `Gene name`,
    Regulation = Regulation_ALL,
    log2FC_0 = log2FC_ALL,   Pval_0 = P.Value_ALL,
    log2FC_1 = `log2FC_HR+`, Pval_1 = `P.Value_HR+`,
    log2FC_2 = `log2FC_HER+`, Pval_2 = `P.Value_HER+`,
    log2FC_3 = log2FC_TNBC,  Pval_3 = P.Value_TNBC
  ) %>% 
  column_to_rownames("feature")


## ===== 在构建 row_split_factor 之前清洗 unchanged 行 =====

# 只保留 up-regulated 和 down-regulated
keep_idx <- row_df$Regulation %in% c("up-regulated", "down-regulated") &
  !is.na(row_df$Regulation)

cat("原始行数:", nrow(row_df), " -> 清洗后行数:", sum(keep_idx), "\n")
print(table(row_df$Regulation, useNA = "ifany"))

# 同步子集 row_df 和 expr_mat
row_df   <- row_df[keep_idx, , drop = FALSE]
expr_mat <- expr_mat[rownames(row_df), , drop = FALSE]


## ===== 9) 颜色映射 =====
# 表达矩阵热图颜色
col_heat <- colorRamp2(c(-2, 0, 2), c("#47bdeb", "white", "#f9632c"))

# 调控方向颜色（用于行分组）
reg_cols <- c(
  "up-regulated" = "#e74c3c",
  "down-regulated" = "#3498db"
)

# P值颜色（-log10转换）
neglog10 <- function(x) -log10(x + .Machine$double.eps)
pval_all <- range(neglog10(c(row_df$Pval_0, row_df$Pval_1, 
                             row_df$Pval_2, row_df$Pval_3)), 
                  na.rm = TRUE, finite = TRUE)
col_pval <- colorRamp2(seq(pval_all[1], pval_all[2], length.out = 4),
                       c("#f1b6da", "#de77ae", "#c51b7d", "#7a0177"))

# # log2FC绝对值缩放（用于点大小）
# scale_fc_fun <- function(fc_vec, min_pt = 1.5, max_pt = 3.5) {
#   abs_fc <- abs(fc_vec)
#   vr <- range(abs_fc, na.rm = TRUE, finite = TRUE)
#   if(diff(vr) == 0) return(rep((min_pt + max_pt) / 2, length(fc_vec)))
#   (abs_fc - vr[1]) / (vr[2] - vr[1]) * (max_pt - min_pt) + min_pt
# }

## ===== 10) 左侧调控方向色块 + 行 split（务必用同一个 row_split_factor）=====

# 先定义 row_split_factor（基于 Regulation）
row_split_factor <- factor(
  row_df$Regulation,
  levels = c("down-regulated", "up-regulated")  # 下调在上
)
row_split_factor <- droplevels(row_split_factor)

split_levels <- levels(row_split_factor)
split_labels <- ifelse(split_levels == "down-regulated", "Down", "Up")


# 去掉没出现的水平，保证切片数真实
row_split_factor <- droplevels(row_split_factor)
split_levels <- levels(row_split_factor)  # 实际存在的分组

# 为每个实际存在的分组生成对应标签
split_labels <- ifelse(split_levels == "up-regulated", "Up", "Down")

# 行侧 block 注释，labels 长度=实际分组数
left_block <- rowAnnotation(
  RegBlock = anno_block(
    gp = gpar(
      fill = unname(reg_cols[split_levels]),  # 只取实际存在的 level 的颜色
      col  = NA
    ),
    labels = split_labels,  # 长度与 split_levels 一致
    labels_gp = gpar(col = "white", fontface = 2, cex = 0.9),
    labels_rot = 0
  ),
  width = unit(15, "mm")
)


# ===== 11) 右侧四个 lollipop 注释 =====
# pt0 <- scale_fc_fun(row_df$log2FC_0)
col0 <- col_pval(neglog10(row_df$Pval_0))
# pt1 <- scale_fc_fun(row_df$log2FC_1)
col1 <- col_pval(neglog10(row_df$Pval_1))
# pt2 <- scale_fc_fun(row_df$log2FC_2)
col2 <- col_pval(neglog10(row_df$Pval_2))
# pt3 <- scale_fc_fun(row_df$log2FC_3)
col3 <- col_pval(neglog10(row_df$Pval_3))

anno_right <- rowAnnotation(
  "Benign" = anno_points(row_df$log2FC_0, pch = 16,
                         gp = gpar(col = col0, fontsize = 6),
                         axis = TRUE, width = unit(18, "mm")),
  
  "HR+" = anno_points(row_df$log2FC_1, pch = 16,
                      gp = gpar(col = col1, fontsize = 6),
                      axis = TRUE, width = unit(18, "mm")),
  
  "HER2+" = anno_points(row_df$log2FC_2, pch = 16,
                        gp = gpar(col = col2, fontsize = 6),
                        axis = TRUE, width = unit(18, "mm")),
  
  "TNBC" = anno_points(row_df$log2FC_3, pch = 16,
                       gp = gpar(col = col3, fontsize = 6),
                       axis = TRUE, width = unit(18, "mm")),
  annotation_name_side = "top"
)

## ===== 12) 图例 =====
lgd_pval <- Legend(
  title = expression(-log[10](P)),
  col_fun = col_pval,
  at = floor(seq(pval_all[1], pval_all[2], length.out = 5)),
  direction = "horizontal"
)
# 
# lgd_fc <- Legend(
#   title = "|log2FC| (size)",
#   type = "points", pch = 16,
#   size = unit(c(1.5, 2.5, 3.5), "mm"),
#   labels = c("Low", "Med", "High"),
#   direction = "horizontal",
#   legend_gp = gpar(col = "grey20")
# )

## ===== 13) 主热图 =====
ht <- Heatmap(
  expr_mat,
  name = "Z-score",
  col = col_heat,
  na_col = "grey90",     # <--- 这里为 NA 专用颜色
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_split_factor,   # 用刚刚 droplevels 之后的这个
  column_split = group,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_max_width = unit(8, "cm"),
  row_names_gp = gpar(fontsize = 8),
  show_column_names = FALSE,
  left_annotation = left_block,   # 这里的 block 和 row_split 保持一致
  right_annotation = anno_right,
  row_title = NULL,
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(
    at = seq(-2, 2, by = 1),
    labels = seq(-2, 2, by = 1),
    direction = "horizontal",
    title_position = "topcenter"
  )
)


## ===== 14) 保存图片 =====
png("FigureS2B.png", width = 3000, height = max(2000, nrow(expr_mat) * 15 + 500), res = 300)

cat("row_split_factor levels:\n"); print(levels(row_split_factor))
cat("anno_block labels:\n");      print(split_labels)
table(row_df$Regulation, useNA = "ifany")


## ===== 15) 绘制热图 =====
draw(ht,
     # heatmap_legend_list = list(lgd_pval, lgd_fc),
     heatmap_legend_list = list(lgd_pval),
     merge_legends = FALSE,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom")

## ===== 16) 给四个 lollipop 每个 split 都加虚线 =====
anns <- c("Benign", "HR+", "HER2+", "TNBC")
n_slices <- length(levels(row_split_factor))  # 与上面的 row_split_factor 完全一致

for (ann in anns) {
  for (s in 1:n_slices) {
    decorate_annotation(ann, slice = s, {
      grid.lines(
        x = unit(c(0.5, 0.5), "npc"),
        y = unit(c(0, 1), "npc"),
        gp = gpar(lty = 3, col = "grey40")
      )
    })
  }
}


dev.off()