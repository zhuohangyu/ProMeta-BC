## ============================================================
## 多组学整合分析总脚本（美化版）：
## WGCNA蛋白模块 × mFuzz代谢簇
## ============================================================

options(stringsAsFactors = FALSE)

library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(purrr)
library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)
library(broom)
library(ggpubr)
library(scales)
library(RColorBrewer)

## ------------------------------------------------------------
## 自定义主题和配色
## -------------------------------------------------------------

# 自定义ggplot主题
theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      # 标题
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      plot.subtitle = element_text(size = rel(1), hjust = 0.5, color = "grey40"),
      
      # 坐标轴
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.text = element_text(size = rel(0.9)),
      axis.line = element_line(color = "black", size = 0.5),
      
      # 图例
      legend.title = element_text(face = "bold", size = rel(1)),
      legend.text = element_text(size = rel(0.9)),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "grey80"),
      
      # 网格线
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.8),
      
      # 背景
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
}

# 配色方案
color_stage <- c("Benign" = "#4DAF4A", "0-II" = "#377EB8", "III-IV" = "#E41A1C")
color_gradient <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                    "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")

## ------------------------------------------------------------
## 1. 读入数据 & 样本对齐
## ------------------------------------------------------------

# 临床
clinical <- read_csv("Clinical.csv")

# 蛋白表达：行 = 蛋白，列 = 样本
prot_expr <- read.csv("../Figure3/Combined_data_filtered.csv",
                      row.names = 1, check.names = FALSE)

# 代谢表达：第一列为 Metabolite，其余列为样本
metab_expr <- read_csv("Sample.csv")

# 代谢矩阵：行 = 代谢物，列 = 样本
metab_mat <- metab_expr %>%
  column_to_rownames("Metabolite") %>%
  as.matrix()

# 三方样本交集
common_samples <- Reduce(
  intersect,
  list(
    colnames(prot_expr),
    colnames(metab_mat),
    clinical$Sample
  )
)

cat("Common samples:", length(common_samples), "\n")

# 按交集子集化
prot_expr <- prot_expr[, common_samples]
metab_mat <- metab_mat[, common_samples]
clinical  <- clinical %>% filter(Sample %in% common_samples)

## ------------------------------------------------------------
## 2. WGCNA 模块 eigengenes
## ------------------------------------------------------------

ME <- read_csv("../Figure3/WGCNA/ModuleEigengenes_All_Samples.csv")
colnames(ME)[1] <- "SampleID"

ME <- ME %>%
  filter(SampleID %in% common_samples) %>%
  arrange(match(SampleID, common_samples))

ME_mat <- ME %>%
  column_to_rownames("SampleID") %>%
  as.matrix()

# 模块名（除去 SampleID）
module_names <- colnames(ME_mat)
cat("Number of modules:", length(module_names), "\n")

## ------------------------------------------------------------
## 3. mFuzz 结果 & 代谢簇 eigengene
## ------------------------------------------------------------

mfuzz_res <- read_csv("./mFuzz/mfuzz_differential_metabolites_results_with_benign.csv")

# 只保留高置信代谢物
mfuzz_res_filt <- mfuzz_res %>%
  filter(High_Confidence)

cat("High confidence metabolites:", nrow(mfuzz_res_filt), "\n")

# 从代谢矩阵中提取这些代谢物
metab_sub <- metab_mat[rownames(metab_mat) %in% mfuzz_res_filt$Metabolite,
                       common_samples, drop = FALSE]

# 为每个簇计算簇eigengene（簇内均值）
cluster_ids <- sort(unique(mfuzz_res_filt$Cluster))
cat("Number of metabolite clusters:", length(cluster_ids), "\n")

cluster_eig <- lapply(cluster_ids, function(cl) {
  mets <- mfuzz_res_filt$Metabolite[mfuzz_res_filt$Cluster == cl]
  mets <- intersect(mets, rownames(metab_sub))
  if (length(mets) == 0) {
    return(rep(NA, length(common_samples)))
  } else if (length(mets) == 1) {
    colMeans(t(metab_sub[mets, common_samples, drop = FALSE]))
  } else {
    colMeans(metab_sub[mets, common_samples, drop = FALSE])
  }
}) %>% do.call(cbind, .)

colnames(cluster_eig) <- paste0("Cluster", sprintf("%02d", cluster_ids))
rownames(cluster_eig) <- common_samples

cluster_names <- colnames(cluster_eig)

## ------------------------------------------------------------
## 4. 层级1：模块 ME × 单个代谢物 相关性
## ------------------------------------------------------------

ME_df  <- as.data.frame(ME_mat)                            # 样本 × 模块
met_df <- as.data.frame(t(metab_sub[, common_samples]))    # 样本 × 代谢物

cor_res <- corAndPvalue(ME_df, met_df, method = "spearman")

ME_met_cor <- cor_res$cor   # 行: 模块, 列: 代谢物
ME_met_p   <- cor_res$p

# FDR 校正
ME_met_p_adj_vec <- p.adjust(as.vector(ME_met_p), method = "BH")
ME_met_p_adj <- matrix(ME_met_p_adj_vec,
                       nrow = nrow(ME_met_p),
                       dimnames = dimnames(ME_met_p))

# 拉成长表并合并 mFuzz 信息
ME_met_long <- ME_met_cor %>%
  as.data.frame() %>%
  rownames_to_column("Module") %>%
  pivot_longer(-Module, names_to = "Metabolite", values_to = "cor") %>%
  mutate(p_adj = as.vector(ME_met_p_adj)) %>%
  left_join(mfuzz_res_filt %>% select(Metabolite, Cluster, Pattern),
            by = "Metabolite")

# 显著关联
sig_ME_met <- ME_met_long %>%
  filter(p_adj < 0.05 & abs(cor) > 0.3)

cat("Significant ME-metabolite correlations:", nrow(sig_ME_met), "\n")

# 导出
write.csv(sig_ME_met,
          "./relation/ME_metabolite_cor_sig.csv",
          row.names = FALSE)

## ------------------------------------------------------------
## 5. 层级2：模块 ME × 代谢簇 eigengene 相关性 + 美化热图
## ------------------------------------------------------------

ME_df2 <- ME_df[common_samples, , drop = FALSE]
MC_df  <- as.data.frame(cluster_eig[common_samples, , drop = FALSE])

cor_res2 <- corAndPvalue(ME_df2, MC_df, method = "spearman")

ME_MC_cor <- cor_res2$cor   # 行: 模块, 列: MC簇
ME_MC_p   <- cor_res2$p

ME_MC_p_adj_vec <- p.adjust(as.vector(ME_MC_p), method = "BH")
ME_MC_p_adj <- matrix(ME_MC_p_adj_vec,
                      nrow = nrow(ME_MC_p),
                      dimnames = dimnames(ME_MC_p))

cat("Significant ME-MC correlations (FDR<0.05):", 
    sum(ME_MC_p_adj < 0.05), "\n")

# 准备热图数据
cor_mat_plot <- ME_MC_cor
p_mat_plot   <- ME_MC_p_adj

# 改进的配色方案
col_fun <- colorRamp2(
  seq(-0.6, 0.6, length.out = 9),
  color_gradient
)

# 添加行列注释
# 模块颜色（13个模块的自定义颜色）
module_colors <- c("#ff7500", "#00ebff", "#2add00", "#007f45", "#563f10", "#bcc2df", "#ff70ff",
                   "#ff55ad", "#ffac00", "#950e1c", "#3b56d4", "#3eaed6", "#96cf6a")

row_ha <- rowAnnotation(
  Module = module_colors,
  col = list(Module = setNames(module_colors, module_colors)),
  show_legend = FALSE,
  simple_anno_size = unit(0.4, "cm")
)

# 列注释：代谢簇
cluster_anno <- data.frame(
  Cluster = colnames(cor_mat_plot)
)
cluster_colors <-  c("#6A3D9A", "#1F78B4", "#33A02C", "#FF7F00", "#E31A1C", "#B15928", "#FFD700", "#E91E63", "#00BCD4", "#607D8B")
names(cluster_colors) <- paste0("Cluster", sprintf("%02d", cluster_ids))

column_ha <- HeatmapAnnotation(
  Cluster = cluster_anno$Cluster,
  col = list(Cluster = cluster_colors),
  show_legend = FALSE,
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 14, fontface = "bold"),
  simple_anno_size     = unit(0.4, "cm")
)


# 创建美化的热图
ht_raw <- Heatmap(
  cor_mat_plot,
  name = "Spearman\nCorrelation",
  col = col_fun,
  
  # 聚类
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  
  # 显示设置
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_side = "bottom",
  row_names_gp      = gpar(fontsize = 16),         # 行名字体
  column_names_gp   = gpar(fontsize = 16),         # 列名字体
  column_names_rot = 45,
  
  # 注释
  left_annotation = row_ha,
  top_annotation = column_ha,
  
  # 单元格函数：添加显著性标记
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (p_mat_plot[i, j] < 0.001) {
      grid::grid.text("***", x, y,
                      gp = grid::gpar(fontsize = 12, col = "black", fontface = "bold"))
    } else if (p_mat_plot[i, j] < 0.01) {
      grid::grid.text("**", x, y,
                      gp = grid::gpar(fontsize = 12, col = "black", fontface = "bold"))
    } else if (p_mat_plot[i, j] < 0.05) {
      grid::grid.text("*", x, y,
                      gp = grid::gpar(fontsize = 12, col = "black", fontface = "bold"))
    }
  },
  
  # 热图尺寸
  width = unit(10, "cm"),
  height = unit(12, "cm"),
  
  # 图例
  heatmap_legend_param = list(
    title = "Spearman\nCorrelation",
    title_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(4, "cm"),
    legend_direction = "vertical",
    title_position = "topcenter"
  )
)

# 绘制热图
draw(ht_raw, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right")

# 保存高质量PDF
pdf("./relation/ME_MC_correlation_heatmap_enhanced.pdf", width = 10, height = 12)
draw(ht_raw, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

# 保存高分辨率PNG
png("./relation/ME_MC_correlation_heatmap_enhanced.png", 
    width = 10, height = 12, units = "in", res = 300)
draw(ht_raw, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

## ------------------------------------------------------------
## 10. 额外可视化：相关性气泡图
## ------------------------------------------------------------

# 将相关性矩阵转换为长格式
cor_long <- ME_MC_cor %>%
  as.data.frame() %>%
  rownames_to_column("Module") %>%
  pivot_longer(-Module, names_to = "Cluster", values_to = "correlation") %>%
  mutate(
    p_adj = as.vector(ME_MC_p_adj),
    abs_cor = abs(correlation),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    sig_level = case_when(
      p_adj < 0.001 ~ "p < 0.001",
      p_adj < 0.01 ~ "p < 0.01",
      p_adj < 0.05 ~ "p < 0.05",
      TRUE ~ "n.s."
    )
  )

# 创建完整的相关性气泡图（包含所有关联）
# 创建完整的相关性气泡图（包含所有关联）
p_bubble_all <- ggplot(cor_long, aes(x = Cluster, y = Module)) +
  geom_point(aes(size = -log10(p_adj), color = correlation), alpha = 0.85) +
  
  scale_color_gradient2(
    low = "#2166AC", 
    mid = "#F7F7F7", 
    high = "#B2182B",
    midpoint = 0,
    limits = c(-0.5, 0.5),
    name = "Spearman\nCorrelation"
  ) +
  
  scale_size_continuous(
    range = c(3, 14),
    breaks = c(1, 2, 3, 4),
    labels = c("1", "2", "3", "4"),  # ← 直接显示数字 *
    name = "-log10(P)"
  ) +
  
  labs(
    title = "Module-Cluster Correlation Network",
    subtitle = paste0("All associations (", nrow(cor_long), " pairs)"),
    x = "Metabolite Cluster",
    y = "Protein Module"
  ) +
  
  theme_publication() +
  theme(
    axis.text.x   = element_text(size = 18, angle = 45, hjust = 1, vjust = 1),
    axis.text.y   = element_text(size = 18),
    axis.title.x  = element_text(size = 20, face = "bold"),
    axis.title.y  = element_text(size = 20, face = "bold"),
    plot.title    = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 18, hjust = 0.5),
    
    panel.grid.major = element_line(color = "grey95", size = 0.5),
    panel.grid.minor = element_blank(),
    
    legend.position = "right",
    legend.title  = element_text(size = 20, face = "bold"),
    legend.text   = element_text(size = 18),
    
    ## 🔥 删除灰色框框
    legend.background   = element_blank(),
    legend.key          = element_blank(),
    legend.box.background = element_blank()
  )

print(p_bubble_all)

ggsave("./relation/Module_Cluster_correlation_bubble_Xcluster_YME.pdf", 
       p_bubble_all, width = 10, height = 10, dpi = 300)

ggsave("Figure5D.png", 
       p_bubble_all, width = 10, height = 10, dpi = 300)



## ------------------------------------------------------------
## 6. 构建回归用数据框：临床 + 模块 + 代谢簇
## ------------------------------------------------------------

data_for_lm <- clinical %>%
  select(Sample, Age, BMI, Cohort, Label, `Clinical stage`) %>%
  rename(SampleID = Sample) %>%
  left_join(
    ME %>% as.data.frame(),
    by = "SampleID"
  ) %>%
  left_join(
    cluster_eig %>%
      as.data.frame() %>%
      rownames_to_column("SampleID"),
    by = "SampleID"
  ) %>%
  mutate(
    Cohort = as.factor(Cohort),
    Label  = as.factor(Label)
  )


