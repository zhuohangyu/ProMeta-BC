# 清理工作空间
rm(list = ls())
gc()


# 加载必要的库
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# 读取数据文件
clinical_data <- read.csv("Clinical.csv", stringsAsFactors = FALSE)
sample_data <- read.csv("Sample.csv", stringsAsFactors = FALSE, 
                        check.names = FALSE)
vip_data <- read.csv("./OPLS_DA_Results/VIP_Ranking.csv", stringsAsFactors = FALSE)

# 数据预处理
print("原始数据维度:")
print(paste("Clinical data:", nrow(clinical_data), "rows"))
print(paste("Sample data:", nrow(sample_data), "rows,", ncol(sample_data), "cols"))

# 检查列名
sample_cols <- colnames(sample_data)
print("Sample data前10个列名:")
print(sample_cols[1:10])

print("Clinical data前10个样本ID:")
print(clinical_data$Sample[1:10])

# 设置代谢物名称为行名
rownames(sample_data) <- sample_data$Metabolite
sample_data <- sample_data[, -1]  # 移除第一列

# 检查样本ID匹配情况
sample_cols_clean <- colnames(sample_data)
clinical_samples <- clinical_data$Sample

# 查找匹配的样本
common_samples <- intersect(clinical_samples, sample_cols_clean)
print(paste("匹配的样本数量:", length(common_samples)))

if (length(common_samples) == 0) {
  # 如果完全不匹配，可能是列名格式问题
  print("尝试修复列名格式...")
  # R可能将特殊字符转换为点
  sample_cols_original <- make.names(clinical_samples)
  common_samples <- intersect(sample_cols_original, sample_cols_clean)
  print(paste("修复后匹配的样本数量:", length(common_samples)))
}

# 只保留匹配的样本
if (length(common_samples) > 0) {
  # 过滤clinical数据，只保留在sample数据中存在的样本
  clinical_data_filtered <- clinical_data[clinical_data$Sample %in% common_samples, ]
  
  # 重新排列sample数据列顺序以匹配clinical数据
  sample_data <- sample_data[, common_samples]
  
  print(paste("过滤后的数据维度:"))
  print(paste("Clinical data:", nrow(clinical_data_filtered), "rows"))
  print(paste("Sample data:", nrow(sample_data), "rows,", ncol(sample_data), "cols"))
  
  # 更新clinical_data
  clinical_data <- clinical_data_filtered
} else {
  
  stop("无法找到匹配的样本ID，请检查数据格式")
}

# 创建分组信息
groups <- clinical_data$Label
benign_samples <- clinical_data$Sample[groups == "Benign"]
malignant_samples <- clinical_data$Sample[groups == "Malignant"]

# 确保分组样本在sample_data中存在
benign_samples <- intersect(benign_samples, colnames(sample_data))
malignant_samples <- intersect(malignant_samples, colnames(sample_data))

print(paste("良性样本数量:", length(benign_samples)))
print(paste("恶性样本数量:", length(malignant_samples)))

# 检查是否有足够的样本进行分析
if (length(benign_samples) < 3 || length(malignant_samples) < 3) {
  stop("每组样本数量不足，无法进行统计分析")
}

# 差异分析函数
perform_differential_analysis <- function(data, group1_samples, group2_samples) {
  results <- data.frame(
    Metabolite = rownames(data),
    Log2FC = numeric(nrow(data)),
    P_value = numeric(nrow(data)),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(data)) {
    # 获取两组数据
    group1_values <- as.numeric(data[i, group1_samples])
    group2_values <- as.numeric(data[i, group2_samples])
    
    # 移除NA值
    group1_values <- group1_values[!is.na(group1_values)]
    group2_values <- group2_values[!is.na(group2_values)]
    
    # 计算均值
    mean1 <- mean(group1_values)
    mean2 <- mean(group2_values)
    
    # 计算Log2 Fold Change (Malignant vs Benign)
    # 因为数据已经是log转换的，直接相减即可
    log2fc <- mean2 - mean1
    
    # 进行t检验
    if (length(group1_values) >= 3 && length(group2_values) >= 3) {
      t_test <- t.test(group2_values, group1_values)
      p_value <- t_test$p.value
    } else {
      p_value <- NA
    }
    
    results$Log2FC[i] <- log2fc
    results$P_value[i] <- p_value
  }
  
  # 计算-log10(p值)
  results$Neg_log10_P <- -log10(results$P_value)
  
  # 移除无效的p值
  results <- results[!is.na(results$P_value), ]
  
  return(results)
}

# 执行差异分析
diff_results <- perform_differential_analysis(sample_data, benign_samples, malignant_samples)

# 合并VIP分数
# 首先检查代谢物名称匹配情况
print("检查VIP数据匹配情况...")
print(paste("差异分析结果:", nrow(diff_results), "个代谢物"))
print(paste("VIP数据:", nrow(vip_data), "个代谢物"))

# 直接按代谢物名称匹配
volcano_data <- merge(diff_results, vip_data, 
                      by.x = "Metabolite", by.y = "Metabolite", all.x = TRUE)

# 检查匹配结果
matched_count <- sum(!is.na(volcano_data$VIP_Score))
print(paste("成功匹配VIP分数的代谢物数量:", matched_count))

# 如果匹配数量较少，尝试模糊匹配
if (matched_count < nrow(diff_results) * 0.5) {
  print("尝试模糊匹配...")
  
  # 创建清理版本的名称用于匹配
  volcano_data$Metabolite_clean <- gsub("[^A-Za-z0-9]", "", volcano_data$Metabolite)
  vip_data$Metabolite_clean <- gsub("[^A-Za-z0-9]", "", vip_data$Metabolite)
  
  # 重新合并
  volcano_data_temp <- merge(volcano_data[, c("Metabolite", "Log2FC", "P_value", "Neg_log10_P", "Metabolite_clean")], 
                             vip_data[, c("Metabolite", "VIP_Score", "Metabolite_clean")], 
                             by = "Metabolite_clean", all.x = TRUE)
  
  # 如果模糊匹配效果更好，使用模糊匹配结果
  matched_count_fuzzy <- sum(!is.na(volcano_data_temp$VIP_Score))
  if (matched_count_fuzzy > matched_count) {
    volcano_data <- volcano_data_temp
    volcano_data$Metabolite <- volcano_data$Metabolite.x
    print(paste("模糊匹配成功，匹配数量:", matched_count_fuzzy))
  }
}

# 处理缺失的VIP分数
missing_vip_count <- sum(is.na(volcano_data$VIP_Score))
if (missing_vip_count > 0) {
  print(paste("为", missing_vip_count, "个代谢物设置默认VIP分数"))
  volcano_data$VIP_Score[is.na(volcano_data$VIP_Score)] <- 0.1
}

# 设置显著性阈值
p_threshold <- 0.05
fc_threshold <- 0.4  # 降低Log2FC阈值到0.5

# 创建显著性分类
volcano_data$Significance <- "Not Significant"
volcano_data$Significance[volcano_data$P_value < p_threshold & volcano_data$Log2FC > fc_threshold] <- "Up-regulated"
volcano_data$Significance[volcano_data$P_value < p_threshold & volcano_data$Log2FC < -fc_threshold] <- "Down-regulated"

# 统计显著差异代谢物
sig_summary <- table(volcano_data$Significance)
print("Differential metabolites summary:")
print(sig_summary)

# 标识需要标记的显著代谢物（VIP分数高的前10个）
volcano_data <- volcano_data[order(volcano_data$VIP_Score, decreasing = TRUE), ]
top_metabolites <- head(volcano_data[volcano_data$Significance != "Not Significant", ], 10)

# 确保代谢物名称列存在
if ("Metabolite.x" %in% colnames(volcano_data)) {
  metabolite_col <- "Metabolite.x"
} else if ("Metabolite.y" %in% colnames(volcano_data)) {
  metabolite_col <- "Metabolite.y"  
} else {
  metabolite_col <- "Metabolite"
}

print(paste("使用代谢物名称列:", metabolite_col))

# 计算上下调代谢物数量
up_count   <- sum(volcano_data$Significance == "Up-regulated")
down_count <- sum(volcano_data$Significance == "Down-regulated")

label_text <- paste0("Up: ", up_count, " Metabolites  |  Down: ", down_count, " Metabolites")


volcano_plot <- ggplot(volcano_data, aes(x = Log2FC, y = Neg_log10_P, 
                                         size = VIP_Score, color = Significance)) +
  geom_point(alpha = 0.56) +
  scale_size_continuous(range = c(1, 8), name = "VIP Score") +
  scale_color_manual(values = c("Up-regulated" = "#FF6B6B", 
                                "Down-regulated" = "#4ECDC4", 
                                "Not Significant" = "#95A5A6"),
                     labels = c("Up-regulated" = "Up",
                                "Down-regulated" = "Down", 
                                "Not Significant" = "Not Sig"),
                     name = "Regulation") +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
             linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_hline(yintercept = -log10(p_threshold), 
             linetype = "dashed", color = "gray50", alpha = 0.7) +
  labs(
    # title = "Volcano Plot of Differential Metabolites",
       # subtitle = "Point size represents VIP score",
       x = "Log2 Fold Change (Malignant vs Benign)",
       y = "-Log10(P-value)") +
  theme_minimal() +
  
  # ✅ 添加顶部统计信息
  annotate(
    "text",
    x = -0.8,
    y = max(volcano_data$Neg_log10_P, na.rm = TRUE) + 2,
    label = label_text,
    size = 9,
    fontface = "bold"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 30),        
        axis.title.y = element_text(size = 30),        
        axis.text.x = element_text(size = 28),         
        axis.text.y = element_text(size = 28),         
        legend.title = element_text(size = 28),        
        legend.text = element_text(size = 26),         
        legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5)),    # 调整颜色图例中点的大小
         size = guide_legend(override.aes = list(color = "black")))  # 调整大小图例的显示

# 添加显著代谢物标签（如果有的话）
if (nrow(top_metabolites) > 0) {
  volcano_plot <- volcano_plot +
    geom_text_repel(data = top_metabolites,
                    aes_string(label = metabolite_col),
                    size = 8.5,
                    box.padding = 0.5,
                    point.padding = 0.5,
                    force = 2,
                    max.overlaps = 10,
                    show.legend = FALSE)
}



# 显示图形
print(volcano_plot)

# 保存图形
ggsave("Figure5B.png", volcano_plot, width = 10, height = 8, dpi = 300)

# 输出显著差异代谢物列表
significant_metabolites <- volcano_data[volcano_data$Significance != "Not Significant", 
                                        c(metabolite_col, "Log2FC", "P_value", "VIP_Score", "Significance")]
significant_metabolites <- significant_metabolites[order(significant_metabolites$P_value), ]

# 重命名列以保持一致性
colnames(significant_metabolites)[1] <- "Metabolite"

print("Significant differential metabolites list (sorted by p-value):")
if (nrow(significant_metabolites) > 0) {
  print(head(significant_metabolites, 20))
} else {
  print("No significant differential metabolites found")
}

# 保存结果
write.csv(volcano_data, "./DEGs/volcano_plot_data.csv", row.names = FALSE)
write.csv(significant_metabolites, "./DEGs/significant_metabolites.csv", row.names = FALSE)

# 创建VIP分数分布图
vip_dist_plot <- ggplot(volcano_data, aes(x = VIP_Score, fill = Significance)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  scale_fill_manual(values = c("Up-regulated" = "#FF6B6B", 
                               "Down-regulated" = "#4ECDC4", 
                               "Not Significant" = "#95A5A6"),
                    name = "Regulation") +
  labs(title = "Distribution of VIP Scores by Regulation Status",
       x = "VIP Score",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(vip_dist_plot)

# 保存VIP分布图
ggsave("./DEGs/vip_distribution.png", vip_dist_plot, width = 10, height = 6, dpi = 300)

# 创建Log2FC分布图
fc_dist_plot <- ggplot(volcano_data, aes(x = Log2FC, fill = Significance)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  scale_fill_manual(values = c("Up-regulated" = "#FF6B6B", 
                               "Down-regulated" = "#4ECDC4", 
                               "Not Significant" = "#95A5A6"),
                    name = "Regulation") +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
             linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Log2 Fold Changes",
       x = "Log2 Fold Change (Malignant vs Benign)",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(fc_dist_plot)

# 保存Log2FC分布图
ggsave("./DEGs/log2fc_distribution.png", fc_dist_plot, width = 10, height = 6, dpi = 300)

# 输出分析总结
cat("\n=== Analysis Summary ===\n")
cat("Total metabolites:", nrow(volcano_data), "\n")
cat("Up-regulated metabolites:", sum(volcano_data$Significance == "Up-regulated"), "\n")
cat("Down-regulated metabolites:", sum(volcano_data$Significance == "Down-regulated"), "\n")
cat("Non-significant metabolites:", sum(volcano_data$Significance == "Not Significant"), "\n")
cat("Mean VIP score:", round(mean(volcano_data$VIP_Score, na.rm = TRUE), 3), "\n")
cat("VIP score range:", round(min(volcano_data$VIP_Score, na.rm = TRUE), 3), "-", 
    round(max(volcano_data$VIP_Score, na.rm = TRUE), 3), "\n")
cat("Log2FC threshold:", fc_threshold, "\n")
cat("P-value threshold:", p_threshold, "\n")

# 详细的上下调分析
cat("\n=== Regulation Analysis Details ===\n")
up_regulated <- volcano_data[volcano_data$Significance == "Up-regulated", ]
down_regulated <- volcano_data[volcano_data$Significance == "Down-regulated", ]

if (nrow(up_regulated) > 0) {
  cat("Up-regulated metabolites (Top 5 by VIP score):\n")
  up_top5 <- head(up_regulated[order(up_regulated$VIP_Score, decreasing = TRUE), ], 5)
  for (i in 1:nrow(up_top5)) {
    cat(sprintf("  %s: Log2FC=%.3f, P-value=%.2e, VIP=%.3f\n", 
                up_top5[i, metabolite_col], up_top5$Log2FC[i], 
                up_top5$P_value[i], up_top5$VIP_Score[i]))
  }
  cat(sprintf("  Log2FC range: %.3f to %.3f\n", 
              min(up_regulated$Log2FC), max(up_regulated$Log2FC)))
  cat(sprintf("  Mean VIP score: %.3f\n", mean(up_regulated$VIP_Score, na.rm = TRUE)))
}

if (nrow(down_regulated) > 0) {
  cat("\nDown-regulated metabolites (Top 5 by VIP score):\n")
  down_top5 <- head(down_regulated[order(down_regulated$VIP_Score, decreasing = TRUE), ], 5)
  for (i in 1:nrow(down_top5)) {
    cat(sprintf("  %s: Log2FC=%.3f, P-value=%.2e, VIP=%.3f\n", 
                down_top5[i, metabolite_col], down_top5$Log2FC[i], 
                down_top5$P_value[i], down_top5$VIP_Score[i]))
  }
  cat(sprintf("  Log2FC range: %.3f to %.3f\n", 
              min(down_regulated$Log2FC), max(down_regulated$Log2FC)))
  cat(sprintf("  Mean VIP score: %.3f\n", mean(down_regulated$VIP_Score, na.rm = TRUE)))
}