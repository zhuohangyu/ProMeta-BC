# 清理环境
rm(list = ls())
gc()

# Load required packages
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(scales)
library(ggrepel)
library(cowplot)
library(stringr)
library(grid)
library(gridExtra)

# Create comprehensive dataset focusing on plasma proteomics and early breast cancer detection
enrichment_data <- data.frame(
  Term = c(
    "Focal adhesion", "Bacterial invasion of epithelial cells",
    "integrin-mediated signaling pathway", "cell adhesion",
    
    # 免疫与炎症相关通路
    "Leukocyte transendothelial migration", "Regulation of actin cytoskeleton",
    
    # 蛋白结合功能 (血浆蛋白相互作用)
    "protein binding", "enzyme binding", "unfolded protein binding",
    "protein-containing complex binding", "integrin binding",
    
    # 细胞骨架相关 (与外泌体释放相关)
    "cytoskeletal protein binding", "actin binding", "actin filament binding",
    
    # GTPase相关 (信号转导)
    "GTP binding", "GTPase activity", "G protein activity",
    
    # 生物过程调控
    "homotypic cell-cell adhesion", "positive regulation of cellular protein localization",
    "regulation of protein localization", "regulation of cellular localization",
    
    # 额外的重要通路
    "regulation of substrate adhesion-dependent cell spreading",
    "positive regulation of lamellipodium organization", # 细胞迁移相关
    
    # 代谢相关通路
    "Glycolysis / Gluconeogenesis", "L-lactate dehydrogenase activity",
    "purine nucleotide metabolic process", "monosaccharide metabolic process",
    "hexose catabolic process", "phosphatidylserine 1-acylhydrolase activity"
  ),
  
  Category = c(
    # 细胞组分
    "KEGG", "KEGG", "Biological process", "Biological process",
    
    # KEGG通路
    "KEGG", "KEGG",
    
    # 分子功能
    "Molecular function", "Molecular function", "Molecular function",
    "Molecular function", "Molecular function",
    
    # 细胞骨架
    "Molecular function", "Molecular function", "Molecular function",
    
    # GTPase
    "Molecular function", "Molecular function", "Molecular function",
    
    # 生物过程
    "Biological process", "Biological process", "Biological process", 
    "Biological process",
    
    "Biological process", "Biological process",
    
    # 代谢相关通路
    "KEGG", "Molecular function",
    "Biological process", "Biological process",
    "Biological process", "Molecular function"
  ),
  
  Regulation = c(
    "up-regulated", "up-regulated", "up-regulated", "up-regulated",
    "up-regulated", "up-regulated",
    "up-regulated", "up-regulated", "up-regulated", "up-regulated", "up-regulated",
    "up-regulated", "up-regulated", "up-regulated",
    "up-regulated", "up-regulated", "up-regulated",
    "up-regulated", "up-regulated", "up-regulated", "up-regulated",
    "up-regulated", "up-regulated",
    "up-regulated", "up-regulated",
    "up-regulated", "up-regulated",
    "up-regulated", "down-regulated"
    
  ),
  
  PValue = c(
    1.35e-15, 1.33e-11, 4.09e-12, 1.12e-9,
    
    # 免疫通路
    2.61e-14, 6.4e-15,
    
    # 蛋白结合
    1.41e-11, 1.04e-8, 1.48e-9, 3.87e-10, 2.15e-10,
    
    # 细胞骨架
    2.88e-11, 1.2e-10, 1.03e-10,
    
    # GTPase
    1.49e-13, 4.06e-11, 4.06e-11,
    
    # 生物过程
    1.01e-15, 8.25e-9, 3.13e-11, 1.68e-10,
    
    1.36e-8, 2e-8,
    
    #代谢通路
    1.49E-07, 0.000105,
    1.05E-05, 8.92E-06,
    8.79E-06, 0.000193
    
  ),
  
  FoldEnrich = c(
    # 与原数据保持一致
    5.72, 8.18, 9.03, 3.59,
    7.41, 5.22,
    1.54, 2.24, 6.27, 2.8, 3.21,
    3.32, 4.89, 5.15,
    4.24, 11.7, 11.7,
    14.39, 4.31, 4.16, 3.67,
    9.64, 12.84,
    6.74, 27,
    3.27, 5.74,
    9.2, 25.13
    
  ),
  
  GeneCount = c(
    # 基因数量
    33, 16, 10, 58,
    21, 34,
    267, 74, 14, 42, 25,
    35, 19, 18,
    55, 9, 9,
    8, 20, 26, 33,
    7, 5,
    12, 3,
    18, 10,
    7, 3
  ),
  
  stringsAsFactors = FALSE
)


# Calculate -log10(p-value)
enrichment_data$NegLogP <- -log10(enrichment_data$PValue)

# Create significance categories for better visualization
enrichment_data$Significance <- case_when(
  enrichment_data$PValue < 1e-20 ~ "Extremely Significant (p < 1e-20)",
  enrichment_data$PValue < 1e-10 ~ "Highly Significant (p < 1e-10)", 
  enrichment_data$PValue < 1e-5 ~ "Significant (p < 1e-5)",
  TRUE ~ "Less Significant"
)

# Create enrichment categories
enrichment_data$EnrichmentLevel <- case_when(
  enrichment_data$FoldEnrich >= 10 ~ "Very High (≥10x)",
  enrichment_data$FoldEnrich >= 5 ~ "High (5-10x)",
  enrichment_data$FoldEnrich >= 2 ~ "Moderate (2-5x)",
  TRUE ~ "Low (<2x)"
)

# Enhanced color palette with gradients
category_colors_enhanced <- c(
  "KEGG" = "#E31A1C", 
  "Biological process" = "#1F78B4", 
  "Cellular component" = "#33A02C", 
  "Molecular function" = "#FF7F00"
)

# ====================================================================
# Enhanced Volcano-style Scatter Plot with Multiple Improvements
# ====================================================================

# 确保已经加载了必要的库，例如 tidyverse (包含 dplyr)
# library(tidyverse) 
# library(dplyr)

# 1. 定义需要强制纳入的代谢相关通路列表 (Term 名字)
metabolic_terms_to_include <- c(
  "Glycolysis / Gluconeogenesis",
  "L-lactate dehydrogenase activity",
  "purine nucleotide metabolic process",
  "monosaccharide metabolic process",
  "hexose catabolic process",
  "phosphatidylserine 1-acylhydrolase activity"
)

# --- 2. 筛选具有平衡代表性的通路（Category Top 3） ---
# 注意：你的原代码注释写的是 Top 6，但代码是 slice_head(n = 3)。
# 我将使用 n = 3，如你的代码所示。
top_pathways_by_category <- enrichment_data %>%
  group_by(Category) %>%
  arrange(PValue) %>%
  slice_head(n = 3) %>%
  ungroup()

# --- 3. 筛选整体显著性最高的通路（Overall Top 5） ---
top_overall <- enrichment_data %>%
  arrange(PValue) %>%
  slice_head(n = 5)

# --- 4. 筛选需要强制纳入的代谢通路 ---
# 从原始数据中根据 Term 名称过滤出这些行
metabolic_pathways <- enrichment_data %>%
  filter(Term %in% metabolic_terms_to_include)

# --- 5. 合并所有结果并去重 ---
top_pathways <- bind_rows(
  top_pathways_by_category,
  top_overall,
  metabolic_pathways # 纳入强制代谢通路
) %>%
  # 根据 Term 列去重，保留第一次出现的完整行数据
  distinct(Term, .keep_all = TRUE)

# 查看最终筛选出的用于标记的通路
# print(top_pathways_for_labeling)

# Create the main scatter plot
enhanced_scatter <- ggplot(enrichment_data, aes(x = FoldEnrich, y = NegLogP)) +
  
  # Add gradient background regions
  annotate("rect", xmin = 0, xmax = 2, ymin = 0, ymax = max(enrichment_data$NegLogP), 
           fill = "lightblue", alpha = 0.05) +
  annotate("rect", xmin = 2, xmax = 5, ymin = 0, ymax = max(enrichment_data$NegLogP), 
           fill = "lightyellow", alpha = 0.05) +
  annotate("rect", xmin = 5, xmax = max(enrichment_data$FoldEnrich), ymin = 0, ymax = max(enrichment_data$NegLogP), 
           fill = "lightcoral", alpha = 0.05) +
  
  # Add significance threshold lines with enhanced styling
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", 
             color = "#D73027", 
             size = 1.2,
             alpha = 0.8) +
  geom_hline(yintercept = -log10(0.01), 
             linetype = "dotdash", 
             color = "#D73027", 
             size = 1,
             alpha = 0.6) +
  geom_vline(xintercept = 2, 
             linetype = "dashed", 
             color = "#4575B4", 
             size = 1.2,
             alpha = 0.8) +
  geom_vline(xintercept = 5, 
             linetype = "dotted", 
             color = "#4575B4", 
             size = 1,
             alpha = 0.6) +
  
  # Add main data points with enhanced aesthetics
  geom_point(aes(size = GeneCount, 
                 fill = Category,
                 color = Category,
                 alpha = NegLogP), 
             shape = 21,
             stroke = 1.5) +
  
  # Highlight extremely significant pathways with glow effect
  geom_point(data = filter(enrichment_data, PValue < 1e-20),
             aes(size = GeneCount + 20, 
                 fill = Category),
             shape = 21,
             color = "white",
             stroke = 2,
             alpha = 0.3) +
  
  # Add enhanced labels for top pathways from each category
  geom_text_repel(data = top_pathways,
                  aes(label = str_wrap(Term, 35),
                      color = Category),
                  size = 7,
                  fontface = "bold",
                  box.padding = 0.6,
                  point.padding = 0.4,
                  segment.color = "gray30",
                  segment.size = 0.6,
                  segment.linetype = "solid",
                  segment.alpha = 0.7,
                  max.overlaps = 40,
                  min.segment.length = 0.1,
                  bg.color = "white",
                  bg.r = 0.12,
                  force = 2,
                  force_pull = 0.1) +
  
  # Enhanced scales
  scale_size_continuous(name = "Gene Count", 
                        range = c(3, 16),
                        breaks = c(10, 50, 100, 200, 300),
                        guide = guide_legend(override.aes = list(color = "black", 
                                                                 stroke = 1,
                                                                 alpha = 0.8))) +
  
  scale_fill_manual(name = "Functional Category", 
                    values = category_colors_enhanced,
                    breaks = names(category_colors_enhanced),
                    labels = c(
                      "KEGG" = "KEGG",
                      "Biological process" = "Biological\nprocess",
                      "Cellular component" = "Cellular\ncomponent",
                      "Molecular function" = "Molecular\nfunction"
                    ),
                    guide = guide_legend(override.aes = list(size = 8, 
                                                             stroke = 1.5,
                                                             alpha = 0.8))) +
  
  scale_color_manual(name = "Functional Category", 
                     values = category_colors_enhanced,
                     breaks = names(category_colors_enhanced),
                     labels = c(
                       "KEGG" = "KEGG",
                       "Biological process" = "Biological\nprocess",
                       "Cellular component" = "Cellular\ncomponent",
                       "Molecular function" = "Molecular\nfunction"
                     )) +
  
  scale_alpha_continuous(name = "Statistical\nSignificance", 
                         range = c(0.4, 1),
                         guide = guide_legend(override.aes = list(size = 6,
                                                                  color = "black"))) +
  
  scale_x_continuous(trans = "log2", 
                     breaks = c(1, 2, 4, 8, 16, 32),
                     labels = c("1", "2", "4", "8", "16", "32"),
                     expand = expansion(mult = c(0.05, 0.15))) +
  
  scale_y_continuous(limits = c(0, max(enrichment_data$NegLogP) * 1.15),
                     expand = expansion(mult = c(0.02, 0.05))) +
  
  # Enhanced theme and labels
  labs(
       # title = "Significant Pathways (Malignant vs Benign)",
       x = "Fold Enrichment (log2scale)",
       y = "-log10 p-value") +
  
  # Professional theme
  theme_minimal(base_size = 22, base_family = "Arial") +
  theme(
    # Plot titles
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, 
                              margin = margin(b = 10), color = "gray20"),
    plot.subtitle = element_text(size = 24, hjust = 0.5, 
                                 margin = margin(b = 15), color = "gray40"),
    plot.caption = element_text(size = 24, color = "gray50", 
                                margin = margin(t = 15), hjust = 0),
    
    # Axes
    axis.text = element_text(size = 28, color = "gray30"),
    axis.title = element_text(size = 32, face = "bold", color = "gray20"),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    
    # Legend
    legend.title = element_text(size = 24, face = "bold", color = "gray20"),
    legend.text = element_text(size = 24, color = "gray30"),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.spacing = unit(0.5, "cm"),
    legend.margin = margin(l = 10),
    
    # Grid and background
    panel.grid.minor = element_line(color = "gray95", size = 0.3),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    
    # Border
    panel.border = element_rect(color = "gray80", fill = NA, size = 0.8)
  ) +
  
  # Custom legend arrangement
  guides(
    fill = guide_legend(title = "Functional Category", 
                        order = 1,
                        override.aes = list(size = 24, stroke = 1.5, alpha = 0.8)),
    size = guide_legend(title = "Gene Count", 
                        order = 2,
                        override.aes = list(fill = "gray60", color = "black", 
                                            stroke = 1, alpha = 0.8)),
    alpha = guide_legend(title = "-log10 p-value", 
                         order = 3,
                         override.aes = list(size = 24, fill = "gray60", color = "black")),
    color = "none"
  )


# Add annotation
enhanced_scatter_with_stats <- enhanced_scatter
  # annotation_custom(
  #   grob = textGrob(stats_text, 
  #                   x = 0.02, y = 0.98, 
  #                   just = c("left", "top"),
  #                   gp = gpar(fontsize = 13, 
  #                             col = "gray40",
  #                             fontfamily = "Arial")),
  #   xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  # )

# Save main enhanced scatter plot
ggsave("Figure2D.png", 
       plot = enhanced_scatter_with_stats, 
       width = 15, height = 12, dpi = 300, bg = "white")
