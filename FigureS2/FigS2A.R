# 清理工作空间
rm(list = ls())
gc()

# 安装并加载必要的包
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("grid")) install.packages("grid")

library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(grid)

# 1. 读取数据
clinical_data <- read.csv("../Figure2/Clinical.csv", stringsAsFactors = FALSE)
sample_data <- read.csv("../Figure2/Sample.csv", stringsAsFactors = FALSE)

# 2. 数据预处理
# 检查数据维度
cat("Clinical数据维度:", dim(clinical_data), "\n")
cat("Sample数据维度:", dim(sample_data), "\n")

# 去除Clinical数据中Molecular Classification为空的行
clinical_filtered <- clinical_data[!is.na(clinical_data$Molecular.Classification) & 
                                     clinical_data$Molecular.Classification != "", ]

cat("去除空值后的Clinical数据维度:", dim(clinical_filtered), "\n")
cat("Molecular Classification分布:\n")
table(clinical_filtered$Molecular.Classification)

# 3. 准备PCA数据
# 检查列名
cat("Sample数据的列名（前10个）:\n")
print(head(colnames(sample_data), 10))

# 检查第一列的名称（代谢物名称列）
first_col_name <- colnames(sample_data)[1]
cat("第一列名称:", first_col_name, "\n")

# 处理重复的代谢物名称 - 删除后出现的重复项
cat("原始代谢物数量:", nrow(sample_data), "\n")

# 使用第一列作为代谢物名称（不管列名是什么）
metabolite_names <- sample_data[, 1]
duplicate_metabolites <- duplicated(metabolite_names)
cat("发现重复代谢物数量:", sum(duplicate_metabolites), "\n")

if(sum(duplicate_metabolites) > 0) {
  duplicated_metabolite_names <- metabolite_names[duplicate_metabolites]
  cat("重复的代谢物名称:", paste(unique(duplicated_metabolite_names), collapse = ", "), "\n")
  
  # 删除重复的代谢物（保留第一个出现的）
  sample_data <- sample_data[!duplicate_metabolites, ]
  metabolite_names <- metabolite_names[!duplicate_metabolites]
  cat("删除重复后的代谢物数量:", nrow(sample_data), "\n")
}

# 设置代谢物名称为行名
rownames(sample_data) <- metabolite_names
sample_expression <- sample_data[, -1]  # 去除第一列代谢物名

# 筛选有效样本的表达数据
valid_samples <- clinical_filtered$Sample
expression_matrix <- sample_expression[, valid_samples]

# 检查是否所有样本都找到了
missing_samples <- setdiff(valid_samples, colnames(sample_expression))
if(length(missing_samples) > 0) {
  cat("警告：以下样本在表达数据中未找到:", missing_samples, "\n")
}

cat("用于PCA分析的数据维度:", dim(expression_matrix), "\n")

# 4. 数据标准化和转置（样本为行，代谢物为列）
# 转置矩阵使样本成为行，代谢物成为列
expression_t <- t(expression_matrix)

# 去除含有NA值的代谢物列
expression_clean <- expression_t[, colSums(is.na(expression_t)) == 0]
cat("去除NA后的数据维度:", dim(expression_clean), "\n")

# 5. 执行PCA分析
pca_result <- prcomp(expression_clean, center = TRUE, scale. = TRUE)

# 获取PCA结果摘要
cat("\nPCA分析摘要:\n")
summary(pca_result)

# 计算方差贡献率
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
cat("\n前10个主成分的方差贡献率:\n")
for(i in 1:min(10, length(variance_explained))) {
  cat(sprintf("PC%d: %.2f%%\n", i, variance_explained[i] * 100))
}

# 6. 准备绘图数据
pca_data <- data.frame(
  Sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  stringsAsFactors = FALSE
)

# 添加分类信息
pca_data <- merge(pca_data, clinical_filtered[, c("Sample", "Label", "Molecular.Classification")], 
                  by = "Sample", all.x = TRUE)

# 将Molecular Classification转换为因子，并设置生物学标签
pca_data$Molecular.Classification <- as.factor(pca_data$Molecular.Classification)

# 创建生物学标签映射
pca_data$Biological_Classification <- factor(pca_data$Molecular.Classification,
                                             levels = c("1", "2", "3"),
                                             labels = c("HR+HER-", "HER2+", "TNBC"))

# 7. 绘制PCA图 - 专注于Molecular Classification
# 设置颜色（更明显的对比色）
colors <- c("HR+HER-" = "#ff1e3f", "HER2+" = "#1eb0ff", "TNBC" = "#33A02C")

# PC1 vs PC2 散点图 - 带置信椭圆
p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Biological_Classification)) +
  # 添加95%置信椭圆
  stat_ellipse(aes(fill = Biological_Classification), 
               level = 0.95, alpha = 0.2, geom = "polygon") +
  # 添加样本点
  geom_point(size = 3, alpha = 0.8) +
  # 设置颜色和填充
  scale_color_manual(values = colors, name = "Molecular Classification") +
  scale_fill_manual(values = colors, name = "Molecular Classification") +
  labs(
    # title = "",
    # subtitle = "",
    x = paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2] * 100, 2), "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.position = "top"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)),
         fill = "none")

ggsave("FigureS2A.png", plot = p1, 
       width = 10, height = 8, dpi = 300)
