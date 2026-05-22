# 清理工作空间
rm(list = ls())
gc()


# 加载必要的R包
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(pheatmap)
library(gridExtra)
library(cowplot)

# 设置工作目录（请根据实际情况修改）
# setwd("your_working_directory")

# 创建输出文件夹
if (!dir.exists("PCA_Results")) {
  dir.create("PCA_Results")
}
if (!dir.exists("PCA_Results/Plots")) {
  dir.create("PCA_Results/Plots")
}

# 1. 读取数据
print("正在读取数据...")
clinical_data <- read.csv("Clinical.csv", stringsAsFactors = FALSE)
sample_data <- read.csv("sample.csv", stringsAsFactors = FALSE)

# 2. 数据预处理
print("正在进行数据预处理...")

# 查看数据基本信息
print(paste("临床数据维度:", paste(dim(clinical_data), collapse = " x ")))
print(paste("样本数据维度:", paste(dim(sample_data), collapse = " x ")))

# 去除表头空白
names(clinical_data) <- trimws(names(clinical_data))
names(sample_data) <- trimws(names(sample_data))

# 转置sample_data，让样本作为行，代谢物作为列
sample_data_t <- t(sample_data[, -1])  # 去除第一列（代谢物名称）
colnames(sample_data_t) <- sample_data[, 1]  # 设置列名为代谢物名称
sample_data_t <- as.data.frame(sample_data_t)

# 添加样本ID作为行名
sample_data_t$Sample <- rownames(sample_data_t)

# 3. 合并数据
print("正在合并临床数据和样本数据...")

# 合并临床数据和转置后的样本数据
merged_data <- merge(clinical_data, sample_data_t, by = "Sample", all.x = TRUE)

# 检查合并后的数据
print(paste("合并后数据维度:", paste(dim(merged_data), collapse = " x ")))
print("Label分布:")
print(table(merged_data$Label))

# 4. 准备PCA数据
print("正在准备PCA分析数据...")

# 提取数值型代谢物数据（从第4列开始，前3列是Sample, Label, Molecular.Classification）
metabolite_data <- merged_data[, 6:ncol(merged_data)]

# 检查并处理缺失值
missing_count <- sum(is.na(metabolite_data))
if (missing_count > 0) {
  print(paste("发现", missing_count, "个缺失值，将使用均值填充"))
  for (i in 1:ncol(metabolite_data)) {
    metabolite_data[is.na(metabolite_data[, i]), i] <- mean(metabolite_data[, i], na.rm = TRUE)
  }
}

# 数据标准化
metabolite_data_scaled <- scale(metabolite_data)
str(metabolite_data)

metabolite_data$X
# 5. 进行PCA分析
print("正在进行PCA分析...")

# 使用prcomp进行PCA
pca_result <- prcomp(metabolite_data_scaled, center = FALSE, scale. = FALSE)

# 计算方差贡献率
variance_explained <- summary(pca_result)$importance
pc_variance <- variance_explained[2, 1:10]  # 前10个主成分的方差贡献率
var_explained <- pc_variance * 100  # 转换为百分比，与用户代码保持一致

# 创建PCA结果数据框
pca_scores <- data.frame(
  Sample = merged_data$Sample,
  Label = merged_data$Label,
  Molecular_Classification = merged_data$Molecular.Classification,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  PC4 = pca_result$x[, 4],
  PC5 = pca_result$x[, 5]
)

# 6. 使用用户指定样式创建可视化图表
print("正在生成用户指定样式的可视化图表...")

# 设置图片参数
plot_width <- 8
plot_height <- 10
plot_dpi <- 300

# 1. 基本PCA散点图（用户指定样式）
p1 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Label)) +
  geom_point(size = 3, alpha = 1) +
  stat_ellipse(aes(fill = Label), alpha = 0.3, geom = "polygon") +
  scale_color_manual(values = c("Benign" = "#1eb0ff", "Malignant" = "#ff1e3f")) +
  scale_fill_manual(values = c("Benign" = "#87CEEB", "Malignant" = "#FFB6C1")) +
  labs(title = "PCA Analysis of Metabolomics Data",
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 32, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 38),
        legend.text = element_text(size = 38),
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 36))

# 保存图1
ggsave("Figure4B.png", p1, 
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

# 2. PC1 vs PC3散点图（相同样式）
p2 <- ggplot(pca_scores, aes(x = PC1, y = PC3, color = Label)) +
  geom_point(size = 3, alpha = 1) +
  stat_ellipse(aes(fill = Label), alpha = 0.3, geom = "polygon") +
  scale_color_manual(values = c("Benign" = "#1eb0ff", "Malignant" = "#ff1e3f")) +
  scale_fill_manual(values = c("Benign" = "#87CEEB", "Malignant" = "#FFB6C1")) +
  labs(title = "PCA Analysis: PC1 vs PC3",
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
ggsave("PCA_Results/Plots/02_PCA_PC1_vs_PC3.png", p2, 
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

# 3. PC2 vs PC3散点图（相同样式）
p3 <- ggplot(pca_scores, aes(x = PC2, y = PC3, color = Label)) +
  geom_point(size = 3, alpha = 1) +
  stat_ellipse(aes(fill = Label), alpha = 0.3, geom = "polygon") +
  scale_color_manual(values = c("Benign" = "#1eb0ff", "Malignant" = "#ff1e3f")) +
  scale_fill_manual(values = c("Benign" = "#87CEEB", "Malignant" = "#FFB6C1")) +
  labs(title = "PCA Analysis: PC2 vs PC3",
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
ggsave("PCA_Results/Plots/03_PCA_PC2_vs_PC3.png", p3, 
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

# 4. 方差贡献率图（匹配用户样式）
variance_df <- data.frame(
  PC = paste0("PC", 1:10),
  Variance = var_explained[1:10]
)

p4 <- ggplot(variance_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#4A90E2", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = paste0(round(Variance, 1), "%")), 
            vjust = -0.5, size = 6, color = "black") +
  labs(title = "Variance Explained by Principal Components",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, max(variance_df$Variance) * 1.15)

# 保存方差图
ggsave("PCA_Results/Plots/04_Variance_Explained.png", p4, 
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

# 5. 累积方差贡献率图（匹配用户样式）
cumulative_variance <- cumsum(var_explained[1:10])
cum_var_df <- data.frame(
  PC = 1:10,
  Cumulative_Variance = cumulative_variance
)

p5 <- ggplot(cum_var_df, aes(x = PC, y = Cumulative_Variance)) +
  geom_line(color = "#4A90E2", size = 2) +
  geom_point(color = "#4A90E2", size = 4) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "#ff1e3f", 
             alpha = 0.8, size = 1) +
  geom_text(aes(label = paste0(round(Cumulative_Variance, 1), "%")), 
            vjust = -0.8, size = 6, color = "black") +
  annotate("text", x = 8.5, y = 82, label = "80% threshold", 
           color = "#ff1e3f", size = 6) +
  labs(title = "Cumulative Variance Explained",
       x = "Number of Principal Components",
       y = "Cumulative Variance Explained (%)") +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24))

# 保存累积方差图
ggsave("PCA_Results/Plots/05_Cumulative_Variance.png", p5, 
       width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

# 6. 载荷图（匹配用户样式）
loadings <- pca_result$rotation[, 1:2]
loadings_df <- data.frame(
  Metabolite = rownames(loadings),
  PC1_loading = loadings[, 1],
  PC2_loading = loadings[, 2]
)

# 计算载荷的重要性
loadings_df$importance <- sqrt(loadings_df$PC1_loading^2 + loadings_df$PC2_loading^2)

# 选择前20个最重要的代谢物
top_loadings <- loadings_df[order(-loadings_df$importance), ][1:20, ]

p6 <- ggplot(loadings_df, aes(x = PC1_loading, y = PC2_loading)) +
  geom_point(alpha = 0.5, color = "grey60", size = 2) +
  geom_point(data = top_loadings, aes(x = PC1_loading, y = PC2_loading), 
             color = "#ff1e3f", size = 3) +
  geom_text_repel(data = top_loadings, 
                  aes(x = PC1_loading, y = PC2_loading, label = Metabolite),
                  size = 4, max.overlaps = 20, 
                  box.padding = 0.3, point.padding = 0.3,
                  color = "#ff1e3f") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6, color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6, color = "grey40") +
  labs(title = "PCA Loadings Plot",
       x = paste0("PC1 Loading (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 Loading (", round(var_explained[2], 2), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24))

# 保存载荷图
ggsave("PCA_Results/Plots/06_Loadings_Plot.png", p6, 
       width = 14, height = 12, dpi = plot_dpi, bg = "white")

# 7. 双变量图（匹配用户样式）
n_arrows <- 15
top_loadings_biplot <- loadings_df[order(-loadings_df$importance), ][1:n_arrows, ]

# 缩放载荷向量
scale_factor <- max(abs(range(c(pca_scores$PC1, pca_scores$PC2)))) / 
  max(abs(range(c(loadings_df$PC1_loading, loadings_df$PC2_loading)))) * 0.7

top_loadings_biplot$PC1_scaled <- top_loadings_biplot$PC1_loading * scale_factor
top_loadings_biplot$PC2_scaled <- top_loadings_biplot$PC2_loading * scale_factor

p7 <- ggplot() +
  # 样本点和椭圆
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = Label), 
             size = 3, alpha = 1) +
  stat_ellipse(data = pca_scores, aes(x = PC1, y = PC2, color = Label, fill = Label), 
               alpha = 0.3, geom = "polygon") +
  scale_color_manual(values = c("Benign" = "#1eb0ff", "Malignant" = "#ff1e3f")) +
  scale_fill_manual(values = c("Benign" = "#87CEEB", "Malignant" = "#FFB6C1")) +
  # 载荷箭头
  geom_segment(data = top_loadings_biplot, 
               aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "#FF6B00", size = 1, alpha = 0.8) +
  # 载荷标签
  geom_text_repel(data = top_loadings_biplot, 
                  aes(x = PC1_scaled, y = PC2_scaled, label = Metabolite),
                  size = 4, color = "#FF6B00", fontface = "bold",
                  max.overlaps = n_arrows) +
  labs(title = "PCA Biplot",
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)"),
       color = "Label") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 26),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24))

# 保存双变量图
ggsave("PCA_Results/Plots/07_Biplot.png", p7, 
       width = 14, height = 12, dpi = plot_dpi, bg = "white")

# 8. 3D PCA图
library(scatterplot3d)
png("PCA_Results/Plots/08_3D_PCA.png", 
    width = plot_width*plot_dpi, height = plot_height*plot_dpi, res = plot_dpi,
    bg = "white")

color_vec <- c("Benign" = "#1eb0ff", "Malignant" = "#ff1e3f")[pca_scores$Label]

s3d <- scatterplot3d(pca_scores$PC1, pca_scores$PC2, pca_scores$PC3,
                     color = color_vec, pch = 19, cex = 1.5,
                     xlab = paste0("PC1 (", round(var_explained[1], 2), "%)"),
                     ylab = paste0("PC2 (", round(var_explained[2], 2), "%)"),
                     zlab = paste0("PC3 (", round(var_explained[3], 2), "%)"),
                     main = "3D PCA Plot",
                     grid = TRUE, box = FALSE,
                     bg = "white",
                     cex.main = 2, cex.lab = 1.5, cex.axis = 1.2)

legend("topright", legend = c("Benign", "Malignant"), 
       col = c("#1eb0ff", "#ff1e3f"), pch = 19, cex = 1.5,
       bg = "white")

dev.off()

# 9. 组合图（用户样式）
combined_user_style <- plot_grid(p1, p2, p4, p5,
                                 labels = c("A", "B", "C", "D"), 
                                 ncol = 2, nrow = 2,
                                 label_size = 20)

ggsave("PCA_Results/Plots/09_Combined_Plots.png", combined_user_style, 
       width = 20, height = 16, dpi = plot_dpi, bg = "white")

# 7. 输出结果摘要
print("=== PCA分析结果摘要 ===")
print(paste("总样本数:", nrow(pca_scores)))
print(paste("恶性样本数:", sum(pca_scores$Label == "Malignant")))
print(paste("良性样本数:", sum(pca_scores$Label == "Benign")))
print(paste("代谢物数量:", ncol(metabolite_data)))
print("")
print("前10个主成分方差贡献率:")
for (i in 1:10) {
  print(paste0("PC", i, ": ", round(var_explained[i], 2), "%"))
}
print(paste("前5个主成分累积贡献率:", round(sum(var_explained[1:5]), 2), "%"))
print(paste("前10个主成分累积贡献率:", round(sum(var_explained[1:10]), 2), "%"))

# 8. 保存结果数据
print("正在保存结果数据...")

# 保存PCA分数
write.csv(pca_scores, "PCA_Results/PCA_scores.csv", row.names = FALSE)

# 保存载荷矩阵
write.csv(loadings_df, "PCA_Results/PCA_loadings.csv", row.names = FALSE)

# 保存前20个最重要的代谢物
write.csv(top_loadings, "PCA_Results/Top20_important_metabolites.csv", row.names = FALSE)

# 保存方差解释信息
variance_summary <- data.frame(
  PC = paste0("PC", 1:10),
  Variance_Explained = round(var_explained[1:10], 2),
  Cumulative_Variance = round(cumsum(var_explained[1:10]), 2)
)
write.csv(variance_summary, "PCA_Results/Variance_Summary.csv", row.names = FALSE)

# 保存样本分组信息
sample_info <- data.frame(
  Sample = pca_scores$Sample,
  Label = pca_scores$Label,
  PC1_score = round(pca_scores$PC1, 4),
  PC2_score = round(pca_scores$PC2, 4),
  PC3_score = round(pca_scores$PC3, 4)
)
write.csv(sample_info, "PCA_Results/Sample_PCA_info.csv", row.names = FALSE)
