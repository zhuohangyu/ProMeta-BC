# 清理工作空间
rm(list = ls())
gc()
# 蛋白质数据PCA分析 - 完整可视化版本
# 根据Clinical文件中的Label进行PCA分析并生成可视化图片

# 安装和加载必要的包
packages <- c("ggplot2", "dplyr", "ggrepel", "RColorBrewer", "gridExtra", "corrplot", "factoextra", "grid")
for(pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(corrplot)
library(factoextra)
library(grid)

# 设置图片保存参数
options(repr.plot.width=12, repr.plot.height=8)

# 读取数据 - 重要：设置check.names=FALSE保持原始列名
clinical_data <- read.csv("Clinical.csv", stringsAsFactors = FALSE)
protein_data <- read.csv("Combined_data_filtered.csv", stringsAsFactors = FALSE, check.names = FALSE)

# 查看数据结构
cat("临床数据维度:", dim(clinical_data), "\n")
cat("蛋白质数据维度:", dim(protein_data), "\n")
cat("Label分布:\n")
print(table(clinical_data$Label))

# 准备蛋白质表达矩阵
protein_matrix <- protein_data[, 4:ncol(protein_data)]
rownames(protein_matrix) <- protein_data[,1]  # 使用第一列作为行名

# 转置矩阵，使样本为行，蛋白质为列
protein_matrix_t <- t(protein_matrix)

# 确保样本顺序一致
sample_names <- rownames(protein_matrix_t)
clinical_data_ordered <- clinical_data[match(sample_names, clinical_data$Sample), ]

# 移除含有NA值的蛋白质
protein_matrix_clean <- protein_matrix_t[, !apply(is.na(protein_matrix_t), 2, any)]
cat("清理后的蛋白质数量:", ncol(protein_matrix_clean), "\n")

# 标准化数据
protein_matrix_scaled <- scale(protein_matrix_clean)

# 执行PCA分析
pca_result <- prcomp(protein_matrix_scaled, center = FALSE, scale. = FALSE)

# 获取PCA结果
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Sample <- rownames(pca_scores)
pca_scores$Label <- clinical_data_ordered$Label

# 计算方差解释比例
var_explained <- summary(pca_result)$importance[2, ] * 100

# 定义颜色 - 按照用户要求的配色方案
point_colors <- c("Benign" = "#1eb0ff", "Malignant" = "#ff1e3f")
fill_colors <- c("Benign" = "#87CEEB", "Malignant" = "#FFB6C1")

# === 图1: 主要的PCA散点图 (PC1 vs PC2) ===
p1 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Label)) +
  geom_point(size = 3, alpha = 1) +
  stat_ellipse(aes(fill = Label), alpha = 0.3, geom = "polygon") +
  scale_color_manual(values = point_colors) +
  scale_fill_manual(values = fill_colors) +
  labs(
    # title = "PCA Analysis of Protein Expression Data",
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 26),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24))

# 保存图1
ggsave("Figure4A.png", p1, width = 8, height = 7, dpi = 300)

# === 图2: PC1 vs PC3 ===
p2 <- ggplot(pca_scores, aes(x = PC1, y = PC3, color = Label)) +
  geom_point(size = 3, alpha = 1) +
  stat_ellipse(aes(fill = Label), alpha = 0.3, geom = "polygon") +
  scale_color_manual(values = point_colors) +
  scale_fill_manual(values = fill_colors) +
  labs(title = "PCA Analysis of Protein Expression Data (PC1 vs PC3)",
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC3 (", round(var_explained[3], 2), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 26),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24))

# 保存图2
ggsave("./PCA/PCA_PC1_vs_PC3.png", p2, width = 12, height = 8, dpi = 300)

# === 图3: PC2 vs PC3 ===
p3 <- ggplot(pca_scores, aes(x = PC2, y = PC3, color = Label)) +
  geom_point(size = 3, alpha = 1) +
  stat_ellipse(aes(fill = Label), alpha = 0.3, geom = "polygon") +
  scale_color_manual(values = point_colors) +
  scale_fill_manual(values = fill_colors) +
  labs(title = "PCA Analysis of Protein Expression Data (PC2 vs PC3)",
       x = paste0("PC2 (", round(var_explained[2], 2), "%)"),
       y = paste0("PC3 (", round(var_explained[3], 2), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 26),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24))

# 保存图3
ggsave("./PCA/PCA_PC2_vs_PC3.png", p3, width = 12, height = 8, dpi = 300)

# === 图4: 方差解释比例柱状图 ===
variance_df <- data.frame(
  PC = factor(paste0("PC", 1:15), levels = paste0("PC", 1:15)),
  Variance = var_explained[1:15],
  Cumulative = cumsum(var_explained[1:15])
)

p4 <- ggplot(variance_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#4CAF50", alpha = 0.8, color = "black") +
  geom_line(aes(y = Cumulative, group = 1), color = "#FF6B35", size = 1.5) +
  geom_point(aes(y = Cumulative), color = "#FF6B35", size = 3) +
  labs(title = "Principal Component Variance Explained",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 22),
        axis.text.y = element_text(size = 22),
        axis.title = element_text(size = 24)) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Cumulative Variance Explained (%)"))

# 保存图4
ggsave("./PCA/PCA_Variance_Explained.png", p4, width = 12, height = 8, dpi = 300)

# === 图5: 3D PCA图的2D投影组合 ===
p5_1 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Label)) +
  geom_point(size = 3, alpha = 1) +
  scale_color_manual(values = point_colors) +
  labs(title = "PC1 vs PC2", 
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)")) +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

p5_2 <- ggplot(pca_scores, aes(x = PC1, y = PC3, color = Label)) +
  geom_point(size = 3, alpha = 1) +
  scale_color_manual(values = point_colors) +
  labs(title = "PC1 vs PC3", 
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC3 (", round(var_explained[3], 2), "%)")) +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

p5_3 <- ggplot(pca_scores, aes(x = PC2, y = PC3, color = Label)) +
  geom_point(size = 3, alpha = 1) +
  scale_color_manual(values = point_colors) +
  labs(title = "PC2 vs PC3", 
       x = paste0("PC2 (", round(var_explained[2], 2), "%)"),
       y = paste0("PC3 (", round(var_explained[3], 2), "%)")) +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

# 创建图例
legend_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Label)) +
  geom_point(size = 4) +
  scale_color_manual(values = point_colors) +
  theme_void() +
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))

p5_combined <- grid.arrange(p5_1, p5_2, p5_3, legend_plot, 
                            layout_matrix = matrix(c(1,2,3,4), nrow=2),
                            top = textGrob("PCA Analysis - Multiple Views", 
                                           gp = gpar(fontsize = 24, fontface = "bold")))

# 保存图5
ggsave("./PCA/PCA_Combined_Views.png", p5_combined, width = 15, height = 10, dpi = 300)

# === 图6: 样本标签分布密度图 ===
p6_1 <- ggplot(pca_scores, aes(x = PC1, fill = Label)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = fill_colors) +
  labs(title = "PC1 Distribution", 
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"), 
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

p6_2 <- ggplot(pca_scores, aes(x = PC2, fill = Label)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = fill_colors) +
  labs(title = "PC2 Distribution", 
       x = paste0("PC2 (", round(var_explained[2], 2), "%)"), 
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

p6_combined <- grid.arrange(p6_1, p6_2, nrow = 2,
                            top = textGrob("Principal Component Distribution", 
                                           gp = gpar(fontsize = 24, fontface = "bold")))

# 保存图6
ggsave("./PCA/PCA_Density_Distribution.png", p6_combined, width = 10, height = 8, dpi = 300)

# === 图7: 载荷图 (Biplot) ===
# 使用factoextra包创建更好的biplot
p7 <- fviz_pca_biplot(pca_result, 
                      col.ind = pca_scores$Label,
                      palette = point_colors,
                      addEllipses = TRUE,
                      ellipse.alpha = 0.3,
                      label = "var",
                      col.var = "black",
                      repel = TRUE,
                      title = "PCA Biplot - Samples and Protein Loadings",
                      legend.title = "Sample Type",
                      pointsize = 3) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 26),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24))

# 保存图7
ggsave("./PCA/PCA_Biplot.png", p7, width = 14, height = 10, dpi = 300)

# === 输出结果摘要 ===
cat("\n=== PCA分析结果摘要 ===\n")
cat("总样本数:", nrow(pca_scores), "\n")
cat("Benign样本数:", sum(pca_scores$Label == "Benign"), "\n")
cat("Malignant样本数:", sum(pca_scores$Label == "Malignant"), "\n")
cat("使用的蛋白质数量:", ncol(protein_matrix_clean), "\n")
cat("前5个主成分的方差解释比例:\n")
for(i in 1:5) {
  cat(paste0("PC", i, ": ", round(var_explained[i], 2), "%\n"))
}
cat("前5个主成分累积方差解释比例:", round(sum(var_explained[1:5]), 2), "%\n")

# 保存数据
write.csv(pca_scores, "./PCA/pca_scores.csv", row.names = FALSE)

# 保存载荷矩阵
loadings <- pca_result$rotation[, 1:5]
loadings_df <- data.frame(
  Protein = rownames(loadings),
  loadings,
  stringsAsFactors = FALSE
)
write.csv(loadings_df, "./PCA/pca_loadings.csv", row.names = FALSE)

# 输出关键蛋白质
top_proteins_pc1 <- head(loadings_df[order(abs(loadings_df$PC1), decreasing = TRUE), ], 10)
top_proteins_pc2 <- head(loadings_df[order(abs(loadings_df$PC2), decreasing = TRUE), ], 10)

cat("\n对PC1贡献最大的前10个蛋白质:\n")
print(top_proteins_pc1[, c("Protein", "PC1")])

cat("\n对PC2贡献最大的前10个蛋白质:\n")
print(top_proteins_pc2[, c("Protein", "PC2")])

# 创建统计摘要表
summary_stats <- data.frame(
  Group = c("Benign", "Malignant", "Total"),
  Count = c(sum(pca_scores$Label == "Benign"), 
            sum(pca_scores$Label == "Malignant"), 
            nrow(pca_scores)),
  PC1_Mean = round(c(mean(pca_scores$PC1[pca_scores$Label == "Benign"]),
                     mean(pca_scores$PC1[pca_scores$Label == "Malignant"]),
                     mean(pca_scores$PC1)), 3),
  PC2_Mean = round(c(mean(pca_scores$PC2[pca_scores$Label == "Benign"]),
                     mean(pca_scores$PC2[pca_scores$Label == "Malignant"]),
                     mean(pca_scores$PC2)), 3),
  PC1_Variance_Explained = paste0(round(var_explained[1], 2), "%"),
  PC2_Variance_Explained = paste0(round(var_explained[2], 2), "%")
)

write.csv(summary_stats, "./PCA/pca_summary_stats.csv", row.names = FALSE)