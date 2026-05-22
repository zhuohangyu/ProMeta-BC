# 清理环境
rm(list = ls())
gc()

# 加载必要的包
suppressMessages({
  library(tidyverse)
  library(openxlsx)
  library(ggpubr)
  library(ggthemes)
  library(ggrepel)
  library(UpSetR)
  library(RColorBrewer)
  library(ggnewscale)
  library(viridis)
  library(scales)
  library(cowplot)
  library(gridExtra)
  library(ComplexHeatmap)
  library(circlize)  # 添加此包以支持colorRamp2函数
})


# 设置主题和颜色
custom_theme <- theme_minimal() +
  theme(
    text = element_text(size = 24),
    plot.title = element_text(size = 26, hjust = 0.5, face = "bold", margin = margin(b = 20)),
    plot.subtitle = element_text(size = 24, hjust = 0.5, color = "grey40"),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 22),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# 自定义颜色方案
colors_regulation <- c("up" = "#E74C3C", "down" = "#3498DB", "not sig" = "grey70")
colors_comparison <- c(
  "Malignant vs Benign" = "#8E44AD", 
  "HR+ vs Benign" = "#2980B9", 
  "HER2+ vs Benign" = "#27AE60",
  "TNBC vs Benign" = "#F39C12"
)

#----------------------------------------------------------------------------------
#  Step 1: 数据加载和预处理
#----------------------------------------------------------------------------------

cat("开始数据加载和预处理...\n")

# 读取上传的差异分析结果
data <- read.csv('Fold_Change.csv', sep = ",", stringsAsFactors = FALSE)

# 查看数据结构
cat("数据维度:", dim(data), "\n")
cat("列名:", colnames(data), "\n")

# 数据清理：移除空行和缺失值行
data <- data %>% 
  filter(!is.na(Gene.name) & Gene.name != "") %>%
  filter(!is.na(log2FC_ALL))

cat("清理后数据行数:", nrow(data), "\n")

#----------------------------------------------------------------------------------
#  Step 2: 数据转换
#----------------------------------------------------------------------------------

cat("进行数据转换...\n")

# 创建长格式数据
long_data <- data %>%
  select(Gene.name, Protein.accession, Protein.description,
         log2FC_ALL, P.Value_ALL, Regulation_ALL,
         log2FC_HR., P.Value_HR., Regulation_HR.,
         log2FC_HER., P.Value_HER., Regulation_HER.,
         log2FC_TNBC, P.Value_TNBC, Regulation_TNBC) %>%
  pivot_longer(
    cols = c(log2FC_ALL, log2FC_HR., log2FC_HER., log2FC_TNBC),
    names_to = "Comparison_raw",
    values_to = "logFC"
  ) %>%
  mutate(
    P.Value = case_when(
      Comparison_raw == "log2FC_ALL" ~ P.Value_ALL,
      Comparison_raw == "log2FC_HR." ~ P.Value_HR.,
      Comparison_raw == "log2FC_HER." ~ P.Value_HER.,
      Comparison_raw == "log2FC_TNBC" ~ P.Value_TNBC,
      TRUE ~ NA_real_
    ),
    Regulation = case_when(
      Comparison_raw == "log2FC_ALL" ~ Regulation_ALL,
      Comparison_raw == "log2FC_HR." ~ Regulation_HR.,
      Comparison_raw == "log2FC_HER." ~ Regulation_HER.,
      Comparison_raw == "log2FC_TNBC" ~ Regulation_TNBC,
      TRUE ~ NA_character_
    ),
    Comparison = case_when(
      Comparison_raw == "log2FC_ALL" ~ "Malignant vs Benign",
      Comparison_raw == "log2FC_HR." ~ "HR+ vs Benign", 
      Comparison_raw == "log2FC_HER." ~ "HER2+ vs Benign",
      Comparison_raw == "log2FC_TNBC" ~ "TNBC vs Benign",
      TRUE ~ Comparison_raw
    )
  ) %>%
  filter(!is.na(logFC) & !is.na(P.Value)) %>%
  rename(Gene = Gene.name) %>%
  mutate(
    logP = -log10(P.Value),
    # 添加adj.P.Value以便后续分析
    adj.P.Value = p.adjust(P.Value, method = "fdr"),
    # 使用更严格的筛选标准
    group = case_when(
      P.Value < 0.05 & logFC > 0.263 ~ "up",
      P.Value < 0.05 & logFC < -0.263 ~ "down", 
      TRUE ~ "not sig"
    ),
    # 创建因子以控制显示顺序
    Comparison = factor(Comparison, 
                        levels = c("Malignant vs Benign", "HR+ vs Benign", 
                                   "HER2+ vs Benign", "TNBC vs Benign"))
  )

# 数据摘要
cat("转换后的数据摘要:\n")
print(table(long_data$Comparison))
print(table(long_data$group))

#----------------------------------------------------------------------------------
#  Step 3: 美化版多比较点图 (Multi-comparison Dot Plot)
#----------------------------------------------------------------------------------

cat("创建美化版多比较点图...\n")

# 计算各个比较组中的上下调基因数量
stats_summary <- long_data %>%
  group_by(Comparison, group) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = group, values_from = count, values_fill = 0)

upregulated_count <- stats_summary$up
downregulated_count <- stats_summary$down
names(upregulated_count) <- stats_summary$Comparison
names(downregulated_count) <- stats_summary$Comparison

# 计算logFC范围
col <- long_data %>%  
  group_by(Comparison) %>%  
  summarise(
    min_logFC = min(logFC, na.rm = TRUE),            
    max_logFC = max(logFC, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Comparison)

# 选择要标注的基因：每个比较组的上调前5和下调前5
top_up_genes <- long_data %>%
  filter(group == "up") %>%
  group_by(Comparison) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 5) %>%
  ungroup()

top_down_genes <- long_data %>%
  filter(group == "down") %>%
  group_by(Comparison) %>%
  arrange(logFC) %>%  # 升序，获取最负的值
  slice_head(n = 5) %>%
  ungroup()

# 合并上调和下调的顶级基因
top_genes <- bind_rows(top_up_genes, top_down_genes)

# 创建标题文本
title_text <- data.frame(
  Comparison = names(upregulated_count),
  upregulated = upregulated_count,
  downregulated = downregulated_count,
  label = paste("Up:", upregulated_count, "\n Down:", downregulated_count)
)

# 让中间标签分三行：XXX \n vs \n Benign
title_text <- title_text %>%
  mutate(label3 = gsub(" vs ", "\nvs\n", Comparison))


# 美化版多比较点图
p1 <- ggplot() +  
  # 绘制背景柱
  geom_col(data = col,
           aes(x = Comparison, y = min_logFC),
           fill = "grey90", alpha = 0.5, width = 0.85) +
  
  # 绘制点图，使用更好的视觉效果
  geom_jitter(data = subset(long_data, abs(logFC) > 0.2),
              aes(x = Comparison, y = logFC, color = group),
              size = 1.2, alpha = 0.7, width = 0.35, stroke = 0) +
  
  # 基因标签 - 分颜色标注上调和下调基因
  geom_text_repel(data = filter(top_genes, group == "up"),
                  aes(x = Comparison, y = logFC, label = Gene),
                  size = 4.5, max.overlaps = 20,
                  color = colors_regulation["up"], 
                  fontface = "bold",
                  force = 0.2,
                  bg.color = "white", bg.r = 0.1) +
  geom_text_repel(data = filter(top_genes, group == "down"),
                  aes(x = Comparison, y = logFC, label = Gene),
                  size = 4.5, max.overlaps = 20,
                  color = colors_regulation["down"], 
                  fontface = "bold",
                  force = 0.2,
                  bg.color = "white", bg.r = 0.1) +
  
  # 比较组标签背景
  geom_tile(
    data = title_text,
    aes(x = Comparison, y = 0),
    fill = colors_comparison[as.character(title_text$Comparison)],
    height = 0.5,
    color = "white",
    alpha = 0.9,
    size = 1
  ) +
  
  # 分割线
  geom_hline(yintercept = c(-0.263, 0.263), linetype = "dashed", 
             color = "black", alpha = 0.5, size = 0.5) +
  
  # 颜色设置 - 保留上调下调图例
  scale_color_manual(values = colors_regulation, name = "") +
  # 隐藏填充颜色图例
  scale_fill_manual(values = colors_comparison, guide = "none")+
  
  # 设置图例中点的大小
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  
  # 新的颜色图层用于文本 - 隐藏图例
  ggnewscale::new_scale_color() +
  scale_color_manual(values = c("Malignant vs Benign" = "white", 
                                "HR+ vs Benign" = "white", 
                                "HER2+ vs Benign" = "white",
                                "TNBC vs Benign" = "black"),
                     guide = "none") +  # 添加这行隐藏文本颜色图例
  
  # 统计数量标签
  geom_text(data = title_text, 
            aes(x = Comparison, y = max(long_data$logFC, na.rm = TRUE) + 0.5, 
                label = label), 
            size = 5, fontface = "bold", hjust = 0.5) +
  
  # 比较组名称
  geom_text(
    data = title_text,
    aes(x = Comparison, y = 0, label = label3),
    size = 4.5,
    fontface = "bold",
    lineheight = 0.9,
    color = "white"
  ) +
  
  # 坐标轴和标签
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.1, 0.15))
  ) +
  
  labs(
    y = "log2 Fold Change",
  ) +
  
  custom_theme +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "top",
    legend.box = "horizontal"
  )

# 保存多比较点图
ggsave("Figure2C.png", p1, 
       width = 6, height = 8, dpi = 300, bg = "white")
