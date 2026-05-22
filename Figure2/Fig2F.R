# =============================================================================
# 分子分型特异性通路可视化分析 - 高级美化版本
# =============================================================================
# 清理环境
rm(list = ls())
gc()

library(tidyverse)
library(ggplot2)
library(scales)
library(ggtext)
library(patchwork)
library(extrafont)

# 2. 定义文件路径
file_info <- list(
  "HER2+" = "pathway_HER2+.csv",
  "HR+"   = "pathway_HR+.csv",
  "TNBC"  = "pathway_TNBC.csv",
  "ALL"   = "pathway_all.csv"
)

# 3. 定义列索引和新列名
COL_INDICES <- c(1, 3, 5, 6, 7)
NEW_NAMES <- c("Category", "Term", "FoldEnrich", "PValue", "FDR")

# 4. 读取所有数据
all_enrichment_df <- tibble()

for (group_name in names(file_info)) {
  file_path <- file_info[[group_name]]
  if(file.exists(file_path)) {
    data <- read_csv(file_path, show_col_types = FALSE) %>%
      select(all_of(COL_INDICES)) %>%
      setNames(NEW_NAMES) %>%
      mutate(Group = group_name)
    
    data <- data %>%
      mutate(
        Significant_P = PValue < 0.05,
        Significant_FDR = FDR < 0.05
      )
    
    all_enrichment_df <- bind_rows(all_enrichment_df, data)
  }
}

all_enrichment_df_unique <- all_enrichment_df %>%
  distinct(Term, Group, .keep_all = TRUE)

# 5. 定义预选的特异性通路列表
her2_terms <- c(
  "JAK-STAT signaling pathway",
  "Signaling pathways regulating pluripotency of stem cells",
  "phosphatidylinositol biosynthetic process",
  "phosphatidylinositol metabolic process",
  "Oxidative phosphorylation"
)

hr_terms <- c(
  "response to estradiol",
  "regulation of intracellular estrogen receptor signaling pathway",
  "chromatin organization",
  "mTOR signaling pathway",
  "amino acid transmembrane transport"
)


tnbc_terms <- c(
  "PD-L1 expression and PD-1 checkpoint pathway in cancer",
  "Cytokine-cytokine receptor interaction",
  "Th1 and Th2 cell differentiation",
  "Glutathione metabolism",
  "Choline metabolism in cancer"
)

target_map <- bind_rows(
  tibble(Term = her2_terms, Specific_Subtype = "HER2+"),
  tibble(Term = hr_terms,   Specific_Subtype = "HR+"),
  tibble(Term = tnbc_terms, Specific_Subtype = "TNBC")
)

# 6. 筛选和格式化数据
plot_data <- all_enrichment_df_unique %>%
  filter(Term %in% target_map$Term) %>%
  left_join(target_map, by = "Term")

# 7. 美化通路名称
format_pathway_name <- function(x) {
  # 1) 先做常规美化
  y <- str_to_title(x)

  # 2) 你的其他规则（保留）
  y <- str_replace_all(y, "\\bErad\\b", "ERAD")
  y <- str_replace_all(y, "\\bGolgi\\b", "Golgi")
  y <- str_replace_all(y, "\\bLdl\\b", "LDL")

  # 3) 关键：强制把这些名字“改回原样”
  # ——注意：PD-L1 / PD-1 往往是在长句子里出现，所以要替换词片段，而不是整句
  y <- str_replace_all(y, fixed("Jak-Stat"), "JAK-STAT")
  y <- str_replace_all(y, fixed("Mtor"), "mTOR")
  y <- str_replace_all(y, fixed("Pd-L1"), "PD-L1")
  y <- str_replace_all(y, fixed("Pd-1"),  "PD-1")

  return(y)
}

plot_data <- plot_data %>%
  filter(Group != "ALL")
plot_data$Term
plot_data <- plot_data %>%
  mutate(
    Term_formatted = format_pathway_name(Term),
    Term_formatted = case_when(
      Term == "PD-L1 expression and PD-1 checkpoint pathway in cancer" ~
        "PD-L1 expression and PD-1\ncheckpoint pathway in cancer",
      
      Term == "regulation of intracellular estrogen receptor signaling pathway" ~
        "Regulation of intracellular estrogen \nreceptor signaling pathway",
      
      Term == "Signaling pathways regulating pluripotency of stem cells" ~
        "Signaling pathways regulating \npluripotency of stem cells",
      
      Term == "phosphatidylinositol biosynthetic process" ~
        "Phosphatidylinositol\nBiosynthetic Process",
      TRUE ~ Term_formatted
    )
  )




# 8. 设置因子排序
term_levels <- c(rev(tnbc_terms), rev(hr_terms), rev(her2_terms))
term_levels <- term_levels[term_levels %in% plot_data$Term]
term_levels_formatted <- format_pathway_name(term_levels)

term_levels_formatted <- case_when(
  term_levels == "PD-L1 expression and PD-1 checkpoint pathway in cancer" ~
    "PD-L1 expression and PD-1\ncheckpoint pathway in cancer",
  term_levels == "regulation of intracellular estrogen receptor signaling pathway" ~
    "Regulation of intracellular estrogen \nreceptor signaling pathway",
  term_levels == "Signaling pathways regulating pluripotency of stem cells" ~
    "Signaling pathways regulating \npluripotency of stem cells",
  term_levels == "phosphatidylinositol biosynthetic process" ~
    "Phosphatidylinositol\nBiosynthetic Process",
  TRUE ~ term_levels_formatted
)






plot_data <- plot_data %>%
  mutate(
    NegLog10P = -log10(PValue + 1e-10),
    NegLog10FDR = -log10(FDR + 1e-10),
    Group = factor(Group, levels = c("HR+", "HER2+", "TNBC")),
    Term_formatted = factor(Term_formatted, levels = term_levels_formatted),
    Specific_Subtype = factor(Specific_Subtype, levels = c("HER2+", "HR+", "TNBC")),
    
    Sig_label = case_when(
      PValue < 0.001 ~ "***",
      PValue < 0.01 ~ "**",
      PValue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    
    Significant = ifelse(PValue < 0.05, "Yes", "No")
  )

# 9. 高级配色方案
# 使用更精致的渐变色
viridis_modified <- c(
  "#440154", "#46337E", "#365C8D", "#277F8E", "#1FA187",
  "#4AC16D", "#9FDA3A", "#FDE725"
)

# 蓝-白-红高对比度渐变
blue_white_red_hd <- c(
  "#08519c", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef",
  "#ffffff",
  "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d"
)

# 专业科学期刊配色
scientific_gradient <- c(
  "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
)

subtype_colors <- c(
  "HER2+" = "#E64B35",  
  "HR+" = "#4DBBD5",    
  "TNBC" = "#00A087"    
)

# 10. 增强版主题
theme_enhanced <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "Arial") +
    theme(
      # 标题
      plot.title = element_text(
        size = 18, 
        face = "bold", 
        hjust = 0,
        margin = margin(b = 8),
        color = "gray20"
      ),
      plot.subtitle = element_text(
        size = 14, 
        hjust = 0, 
        color = "gray40",
        margin = margin(b = 15),
        lineheight = 1.2
      ),
      
      # 坐标轴
      axis.text.x = element_text(
        size = 18,
        color = "black",
        face = "bold",
        margin = margin(t = 5)
      ),
      axis.text.y = element_text(
        size = 18,
        color = "black",
        lineheight = 1.3,
        margin = margin(r = 5)
      ),
      axis.title = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.5, color = "black"),
      axis.ticks.y = element_blank(),
      
      # 网格线
      panel.grid.major.y = element_line(
        linewidth = 0.3, 
        color = "gray85",
        linetype = "solid"
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      
      # 背景
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # 图例
      legend.position = "right",
      legend.title = element_text(
        size = 18, 
        face = "bold",
        color = "gray20"
      ),
      legend.text = element_text(
        size = 17,
        color = "gray30"
      ),
      legend.key.size = unit(0.5, "cm"),
      legend.spacing.y = unit(0.4, "cm"),
      legend.box.spacing = unit(0.5, "cm"),
      legend.background = element_rect(
        fill = alpha("white", 0.8),
        color = "gray80",
        linewidth = 0.5
      ),
      legend.margin = margin(8, 8, 8, 8),
      
      # 边距
      plot.margin = margin(15, 15, 15, 40)
    )
}

# =============================================================================
# 11. 方案1: 高级美化版 - 带背景色块分隔
# =============================================================================

# 创建背景矩形数据
background_data <- plot_data %>%
  distinct(Term_formatted, Specific_Subtype) %>%
  group_by(Specific_Subtype) %>%
  summarize(
    ymin = min(as.numeric(Term_formatted)) - 0.5,
    ymax = max(as.numeric(Term_formatted)) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(
    fill_color = case_when(
      Specific_Subtype == "HER2+" ~ alpha("#E64B35", 0.05),
      Specific_Subtype == "HR+" ~ alpha("#4DBBD5", 0.05),
      Specific_Subtype == "TNBC" ~ alpha("#00A087", 0.05)
    )
  )

p_enhanced_v1 <- ggplot(plot_data, aes(x = Group, y = Term_formatted)) +
  # 添加背景色块
  geom_rect(
    data = background_data,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = Specific_Subtype),
    inherit.aes = FALSE,
    alpha = 0.08,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = subtype_colors,
    guide = "none"
  ) +
  
  # 添加主气泡
  ggnewscale::new_scale_fill() +
  
  geom_point(
    aes(size = FoldEnrich, fill = NegLog10P, alpha = Significant), 
    shape = 21, 
    color = "black", 
    stroke = 0.5
  ) +
  
  # # 添加亚型标签（带圆角背景）
  # geom_label(
  #   data = plot_data %>% 
  #     filter(Group == "ALL") %>%
  #     distinct(Term_formatted, Specific_Subtype),
  #   aes(x = 0.35, y = Term_formatted, 
  #       label = Specific_Subtype, 
  #       color = Specific_Subtype),
  #   hjust = 0.5,
  #   size = 3.2,
  #   fontface = "bold",
  #   label.padding = unit(0.25, "lines"),
  #   label.size = 0.5,
  #   fill = "white",
  #   show.legend = FALSE
  # ) +
  
  scale_color_manual(values = subtype_colors) +
  
  scale_alpha_manual(
    values = c("Yes" = 1.0, "No" = 0.35),
    name = "Significance",
    labels = c("Yes" = "P < 0.05", "No" = "P ≥ 0.05"),
    guide = guide_legend(
      order = 3, 
      override.aes = list(size = 6, fill = "gray50")
    )
  ) +
  
  scale_size_continuous(
    range = c(3, 13),
    name = "Fold enrichment",
    breaks = c(2, 4, 6, 8, 10),
    guide = guide_legend(
      order = 1, 
      override.aes = list(fill = "gray50", alpha = 1, color = "black"),
      label.position = "right"
    )
  ) +
  
  scale_fill_gradientn(
    colors = blue_white_red_hd,
    name = expression(bold(-log[10]*"("*italic(P)*"-value)")),
    breaks = pretty_breaks(n = 5),
    guide = guide_colorbar(
      order = 2,
      barwidth = 1.2,
      barheight = 10,
      frame.colour = "black",
      frame.linewidth = 0.5,
      ticks.colour = "black",
      ticks.linewidth = 0.5
    )
  ) +
  
  scale_x_discrete(expand = expansion(add = c(0.8, 0.3))) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  
  # labs(
  #   title = "Molecular Subtype-Specific Pathway Enrichment in Breast Cancer",
  #   subtitle = "Functionally distinct pathways characterizing HER2+, HR+, and TNBC subtypes"
  # ) +
  
  theme_enhanced() +
  theme(
  axis.text.x = element_text(
    color = c(
      subtype_colors["HR+"],
      subtype_colors["HER2+"],
      subtype_colors["TNBC"]
    )
  )
)


ggsave(
  "Figure2F.png",
  plot = p_enhanced_v1,
  width = 10,
  height = 12,
  dpi = 600,
  bg = "white"
)