###############################################
# Figure 2A: Benign vs Stage 0–I vs Malignant
# 目标：
# 1) 图中三组同时展示：
#    - Benign
#    - Stage 0–I（恶性早期，0/0A/0B/I/IA/IB）
#    - Malignant（所有恶性，包含早期+非早期）
# 2) 仅做两次统计比较（不做 ANOVA）：
#    - Stage 0–I vs Benign
#    - Malignant vs Benign
###############################################

# 1. 清理工作空间 --------------------------------
rm(list = ls())
gc()

# 2. 加载必要的包 --------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# 3. 读取数据 ------------------------------------
clinical_data <- read.csv("Clinical.csv", stringsAsFactors = FALSE)
sample_data   <- read.csv("Sample.csv",   stringsAsFactors = FALSE)

# 统一样本名称格式：将 "-" 替换为 "_"
clinical_data$Sample  <- gsub("-", "_", clinical_data$Sample)
colnames(sample_data) <- gsub("-", "_", colnames(sample_data))

cat("Clinical 数据维度:", dim(clinical_data), "\n")
cat("Sample   数据维度:", dim(sample_data),   "\n")

cat("标签分布:\n")
print(table(clinical_data$Label))

# 检查是否有 Stage 列
if (!"Stage" %in% colnames(clinical_data)) {
  stop("Clinical.csv 中未找到列 'Stage'，请检查临床信息文件。")
}

# 4. 构建样本映射 --------------------------------
sample_label_map <- setNames(clinical_data$Label, clinical_data$Sample)
sample_stage_map <- setNames(clinical_data$Stage, clinical_data$Sample)

# Sample.csv 第一列默认是蛋白/基因名，其余列是样本
sample_columns <- colnames(sample_data)[-1]

# 只保留 clinical 中有标签的样本
valid_samples <- sample_columns[sample_columns %in% names(sample_label_map)]
cat("有效匹配样本数:", length(valid_samples), "\n")

# 5. 计算每个样本检出的蛋白数 --------------------
calculate_protein_counts <- function(data, samples) {
  protein_counts <- sapply(samples, function(sample_name) {
    sample_values <- data[[sample_name]]
    # 非 NA 且非 0 的数量
    sum(!is.na(sample_values) & sample_values != 0, na.rm = TRUE)
  })
  return(protein_counts)
}

protein_counts <- calculate_protein_counts(sample_data, valid_samples)

# 6. 构建基础结果数据框（每个样本只出现一次） -----
result_df <- data.frame(
  Sample        = names(protein_counts),
  Protein_Count = as.numeric(protein_counts),
  Label         = sample_label_map[names(protein_counts)],
  Stage_raw     = sample_stage_map[names(protein_counts)],
  stringsAsFactors = FALSE
)

# 标准化 Stage 表达形式：去掉 "Stage " 前缀，统一大写，去空格
result_df$Stage_clean <- toupper(trimws(gsub("^Stage\\s*", "", result_df$Stage_raw)))

# 定义“更早期”分期：0 / 0A / 0B / I / IA / IB
early_stages <- c("0", "0A", "0B", "I", "IA", "IB")

# 标记是否恶性/是否早期恶性
result_df <- result_df %>%
  mutate(
    is_malignant = (Label == "Malignant"),
    is_early     = (Label == "Malignant" & Stage_clean %in% early_stages)
  )

# ---- 关键：构建绘图数据 plot_df（允许重复行）----
# 目的是让 Stage 0–I 的样本同时出现在：
#   - “Stage 0–I”
#   - “Malignant”
# 注意：我们不做 ANOVA，只与 Benign 比较，因此这种展示方式可接受
benign_df <- result_df %>%
  filter(Label == "Benign") %>%
  mutate(Group = "Benign")

early_df <- result_df %>%
  filter(is_early) %>%
  mutate(Group = "Stage 0–I ")

allmal_df <- result_df %>%
  filter(is_malignant) %>%
  mutate(Group = "Malignant ")

plot_df <- bind_rows(benign_df, early_df, allmal_df)

plot_df$Group <- factor(
  plot_df$Group,
  levels = c("Benign", "Stage 0–I ", "Malignant ")
)

cat("三组（绘图数据）样本数量（注意早期恶性会重复出现）:\n")
print(table(plot_df$Group))

# 7. 统计摘要（用于显示每组 n 等；n 用去重 Sample） ----------------
summary_stats <- plot_df %>%
  group_by(Group) %>%
  summarise(
    n_unique        = n_distinct(Sample),
    mean_proteins   = mean(Protein_Count),
    median_proteins = median(Protein_Count),
    sd_proteins     = sd(Protein_Count),
    min_proteins    = min(Protein_Count),
    max_proteins    = max(Protein_Count),
    .groups         = "drop"
  )

cat("\n=== 统计摘要（按 Group；n 为 unique sample）===\n")
print(summary_stats)

# 8. 仅做两次比较（不做 ANOVA） -------------------
# Stage 0–I vs Benign
# Malignant vs Benign
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))

p1 <- wilcox.test(early_df$Protein_Count, benign_df$Protein_Count)$p.value
p2 <- wilcox.test(allmal_df$Protein_Count, benign_df$Protein_Count)$p.value

p1_label <- fmt_p(p1)
p2_label <- fmt_p(p2)

cat("\nWilcoxon p (Stage 0–I vs Benign):", p1, " -> ", p1_label, "\n")
cat("Wilcoxon p (Malignant vs Benign):", p2, " -> ", p2_label, "\n")

# 9. 绘图：小提琴 + 箱线 + 蜂群点 -----------------
y_max <- max(plot_df$Protein_Count, na.rm = TRUE)
group_levels <- levels(plot_df$Group)

# 顶部 n（去重 sample）
n_vec <- summary_stats$n_unique[match(group_levels, summary_stats$Group)]

# 只画两条比较括号
y_base <- y_max * 1.28
y_step <- y_max * 0.10

pairwise_plot_df <- data.frame(
  Group1  = c("Stage 0–I ", "Malignant "),
  Group2  = c("Benign",     "Benign"),
  y       = c(y_base,       y_base + y_step),
  p_label = c(p1_label,     p2_label),
  stringsAsFactors = FALSE
)

p <- ggplot(plot_df, aes(x = Group, y = Protein_Count)) +
  # 小提琴图：分布
  geom_violin(
    aes(fill = Group),
    alpha = 0.6,
    trim  = FALSE,
    width = 0.7,
    color = "black",
    size  = 0.3
  ) +
  # 箱线图：中位数 & 四分位
  geom_boxplot(
    aes(fill = Group),
    alpha         = 0.8,
    width         = 0.25,
    outlier.size  = 1.2,
    outlier.alpha = 0.7,
    color         = "black",
    size          = 0.3
  ) +
  # 散点图（优先用蜂群点；未安装 ggbeeswarm 则退回 jitter）
  {
    if (require(ggbeeswarm, quietly = TRUE)) {
      ggbeeswarm::geom_quasirandom(
        aes(color = Group),
        size  = 0.8,
        alpha = 0.6,
        width = 0.18
      )
    } else {
      geom_jitter(
        aes(color = Group),
        width  = 0.08,
        height = 0,
        size   = 0.8,
        alpha  = 0.6
      )
    }
  } +
  # 颜色设置
  scale_fill_manual(
    values = c(
      "Benign"         = "#87CEEB",
      "Stage 0–I "     = "#FFB6C1",
      "Malignant " = "#CD5C5C"
    ),
    name = ""
  ) +
  scale_color_manual(
    values = c(
      "Benign"         = "#4682B4",
      "Stage 0–I "     = "#E9967A",
      "Malignant " = "#8B0000"
    ),
    name = ""
  ) +
  labs(
    x = "",
    y = "The number of proteins",
  ) +
  coord_cartesian(ylim = c(0, y_max * 1.5)) +
  theme_minimal() +
  theme(
    axis.text.x     = element_text(size = 22, color = "black"),
    axis.text.y     = element_text(size = 24, color = "black"),
    axis.title.y    = element_text(size = 28, face = "bold", color = "black"),
    plot.title      = element_text(size = 30, face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle   = element_text(size = 22, hjust = 0.5, color = "darkblue"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.3),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(20, 20, 20, 20),
    axis.line          = element_line(color = "black", size = 0.5)
  ) +
  # 顶部标注各组样本数 n（unique sample）
  annotate(
    "text",
    x     = seq_along(group_levels),
    y     = rep(y_max * 1.15, length(group_levels)),
    label = paste0("n = ", n_vec),
    size  = 6,
    fontface = "bold",
    color    = "darkblue"
  )

# 10. 添加两条比较括号与P值 -----------------------
# 计算 x 位置
pairwise_plot_df$x1 <- as.numeric(factor(pairwise_plot_df$Group1, levels = group_levels))
pairwise_plot_df$x2 <- as.numeric(factor(pairwise_plot_df$Group2, levels = group_levels))

p <- p +
  # 括号横线
  geom_segment(
    data = pairwise_plot_df,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    size = 0.6,
    color = "black"
  ) +
  # 括号竖线（左）
  geom_segment(
    data = pairwise_plot_df,
    aes(x = x1, xend = x1, y = y, yend = y - y_step*0.3),
    inherit.aes = FALSE,
    size = 0.6,
    color = "black"
  ) +
  # 括号竖线（右）
  geom_segment(
    data = pairwise_plot_df,
    aes(x = x2, xend = x2, y = y, yend = y - y_step*0.3),
    inherit.aes = FALSE,
    size = 0.6,
    color = "black"
  ) +
  # P 值文本
  geom_text(
    data = pairwise_plot_df,
    aes(x = (x1 + x2)/2, y = y + y_step*0.35, label = paste0("P = ", p_label)),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  )

# 显示图形
print(p)

# 11. 保存图形 & 数据 -----------------------------
ggsave("Figure2A.png", plot = p, width = 6.5, height = 6, dpi = 300)

# 建议：导出“原始每个样本只出现一次”的数据（不重复）
write.csv(result_df, "Figure2A_samples.csv", row.names = FALSE)

cat("\n=== 分析完成（Benign / Stage 0–I / Malignant；仅与 Benign 比较；无 ANOVA）===\n")
cat("图形已保存为: Figure2A.png\n")
cat("基础结果（unique样本）已保存为: Figure2A_base_unique_samples.csv\n")
cat("绘图数据（含重复）已保存为: Figure2A_plotdata_with_duplicates.csv\n\n")
cat("图形说明：\n")
cat("1. 小提琴图：显示各组蛋白检出数的分布形态。\n")
cat("2. 箱线图：展示中位数及四分位范围。\n")
cat("3. 散点/蜂群点：展示每个样本的具体数值。\n")
cat("4. 顶部文本：显示各组样本量 n（按 Sample 去重）。\n")
cat("5. 括号与P值：仅比较 Stage 0–I vs Benign；Malignant vs Benign（Wilcoxon）。\n")
