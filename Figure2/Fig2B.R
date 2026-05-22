# 清理工作空间
rm(list = ls())
gc()

# 加载必要的包
library(readr)
library(dplyr)
library(ggplot2)
library(factoextra)
library(corrplot)
library(plotly)
library(rgl)
library(car)  # 用于置信椭球计算
library(htmlwidgets)  # 用于保存交互式图表

# 读取数据文件
sample_data <- read_csv("Sample.csv")
clinical_data <- read_csv("Clinical.csv")
head(clinical_data)

# 提取基因表达数据（去除Gene name列）
gene_names <- sample_data[[1]]  # 第一列是基因名
expression_data <- sample_data[, -1]  # 去除第一列

# 检查并处理重复的基因名
cat("原始基因数量:", length(gene_names), "\n")
cat("唯一基因数量:", length(unique(gene_names)), "\n")

# 如果有重复基因名，保留第一个出现的
if(any(duplicated(gene_names))) {
  cat("发现重复基因名，正在处理...\n")
  # 找到重复的基因名
  duplicated_genes <- gene_names[duplicated(gene_names)]
  cat("重复的基因名:", paste(unique(duplicated_genes), collapse = ", "), "\n")
  
  # 保留第一次出现的基因，删除后续重复的
  keep_indices <- !duplicated(gene_names)
  gene_names <- gene_names[keep_indices]
  expression_data <- expression_data[keep_indices, ]
  
  cat("处理后基因数量:", length(gene_names), "\n")
}

# 处理列名不一致问题：将下划线替换为连字符
colnames(expression_data) <- gsub("_", "-", colnames(expression_data))

# 转置数据：样本作为行，基因作为列
expression_matrix <- t(expression_data)
colnames(expression_matrix) <- gene_names

# 将矩阵转换为数据框
expression_df <- as.data.frame(expression_matrix)
expression_df$Sample <- rownames(expression_df)

# 合并表达数据和临床数据
merged_data <- inner_join(clinical_data, expression_df, by = "Sample")

cat("合并后数据维度:", dim(merged_data), "\n")
cat("合并成功的样本数:", nrow(merged_data), "\n")

# 检查Label分布
table(merged_data$Label)

# 准备PCA分析的数据
# 去除Sample和Label列，只保留基因表达数据
pca_data <- merged_data[, -(1:2)]

# 强制转换所有列为数值型
pca_data_numeric <- data.frame(lapply(pca_data, function(x) {
  if(is.character(x) || is.factor(x)) {
    # 尝试将字符转换为数值，无法转换的设为NA
    as.numeric(as.character(x))
  } else {
    as.numeric(x)
  }
}))

# 恢复列名
colnames(pca_data_numeric) <- colnames(pca_data)
cat("转换后数值型列数:", sum(sapply(pca_data_numeric, is.numeric)), "\n")

# 移除含有NA值的基因
na_counts <- colSums(is.na(pca_data_numeric))
cat("含NA值的基因数:", sum(na_counts > 0), "\n")
if(sum(na_counts > 0) > 0) {
  cat("含NA值最多的前5个基因:\n")
  top_na <- sort(na_counts, decreasing = TRUE)[1:min(5, sum(na_counts > 0))]
  print(top_na)
}

pca_data_clean <- pca_data_numeric[, na_counts == 0]
cat("清理后的基因数量:", ncol(pca_data_clean), "\n")

# 移除方差为0的基因
var_genes <- apply(pca_data_clean, 2, var, na.rm = TRUE)
zero_var_count <- sum(var_genes == 0, na.rm = TRUE)
cat("零方差基因数:", zero_var_count, "\n")

pca_data_final <- pca_data_clean[, var_genes > 0 & !is.na(var_genes)]
cat("最终用于PCA的基因数量:", ncol(pca_data_final), "\n")

# 执行PCA分析
pca_result <- prcomp(pca_data_final, scale. = TRUE, center = TRUE)

# 查看PCA结果摘要
summary(pca_result)

# 创建PCA结果数据框
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Sample <- merged_data$Sample
pca_scores$Label <- merged_data$Label

# 计算方差解释比例
var_explained <- summary(pca_result)$importance[2, ] * 100

# 1. 基本PCA散点图
p1 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Label)) +
  geom_point(size = 3, alpha = 1) +
  stat_ellipse(aes(fill = Label), alpha = 0.3, geom = "polygon") +
  scale_color_manual(values = c("Benign" = "#1eb0ff", "Malignant" = "#ff1e3f")) +  # 反转颜色
  scale_fill_manual(values = c("Benign" = "#87CEEB", "Malignant" = "#FFB6C1")) +   # 反转颜色
  labs(#title = "PCA Analysis of Gene Expression Data",
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 34),
        legend.text = element_text(size = 34),
        axis.title = element_text(size = 32),  # 坐标轴标题字体大小
        axis.text = element_text(size = 32))   # 坐标轴刻度字体大小

# 保存所有图片
ggsave("Figure2B.png", plot = p1, width = 8, height = 8, dpi = 300)
