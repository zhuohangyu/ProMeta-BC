# =============================================================================
# Pathway Enrichment Visualization - 模仿参考图片配色风格
# 蓝-白-橙-红渐变 + 灰色气泡大小图例
# =============================================================================

# 1. 加载必要的R包
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(scales)

# 2. 定义文件路径 (请根据实际路径修改)
file_info <- list(
  "HER2+" = "pathway_HER2+.csv",
  "HR+"   = "pathway_HR+.csv",
  "TNBC"  = "pathway_TNBC.csv",
  "ALL"   = "pathway_all.csv"
)

# 3. 定义列索引和新列名：把第 2 列 Regulation 加进来
COL_INDICES <- c(1, 2, 3, 5, 6, 7)
NEW_NAMES   <- c("Category", "Regulation", "Term", "FoldEnrich", "PValue", "FDR")

# 4. 读取和处理数据
all_enrichment_list <- list()
all_enrichment_df <- tibble()

message("--- 数据处理与筛选 ---")
for (group_name in names(file_info)) {
  file_path <- file_info[[group_name]]
  
  data <- read_csv(file_path, show_col_types = FALSE) %>%
    select(all_of(COL_INDICES)) %>%
    setNames(NEW_NAMES) %>%
    mutate(
      Group = group_name,
      Regulation = case_when(
        Regulation %in% c("up-regulated", "up", "+")     ~ "Up",
        Regulation %in% c("down-regulated", "down", "-") ~ "Down",
        TRUE ~ Regulation
      ),
      Term = case_when(
        Regulation == "Up"   ~ paste0(Term, "\n (up-regulated)"),
        Regulation == "Down" ~ paste0(Term, "\n (down-regulated)"),
        TRUE ~ Term
      )
    ) %>%
    # 🚫 去除所有 COVID 相关通路
    filter(!str_detect(Term, regex("COVID", ignore_case = TRUE)))
  
  
  # 筛选策略不变
  if (group_name == "ALL") {
    data_filtered <- data %>%
      filter(FDR < 0.05, Category != "Pfam")
  } else {
    data_filtered <- data %>%
      filter(PValue < 0.05, Category != "Pfam")
  }
  
  # 注意：此时 Term 已经区分 ↑ / ↓，交集就会按“同一通路 + 同一调控方向”来算
  all_enrichment_list[[group_name]] <- data_filtered$Term
  all_enrichment_df <- bind_rows(all_enrichment_df, data_filtered)
}


# 去重
all_enrichment_df_unique <- all_enrichment_df %>%
  distinct(Term, Group, .keep_all = TRUE)

# 检查筛选结果
message("筛选后的通路数量:")
print(sapply(all_enrichment_list, length))

# 5. 找出共性通路
common_terms <- Reduce(intersect, all_enrichment_list)
message(paste("共性通路数量:", length(common_terms)))

# =============================================================================
# 定义参考图片风格的主题和配色
# =============================================================================

# 自定义蓝-白-红配色 (模仿参考图片)
blue_white_red_colors <- c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", 
                           "#f7f7f7", 
                           "#fddbc7", "#f4a582", "#d6604d", "#b2182b")

# 定义参考图片风格主题
theme_reference <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      # 标题
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 14, color = "gray40", margin = margin(b = 10)),
      
      # 坐标轴
      axis.text.x = element_text(color = "black", size = 14, face = "bold", angle = 45, hjust = 1),
      axis.text.y = element_text(color = "black", size = 13),
      axis.title = element_text(face = "bold", size = 11),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.8),
      
      # 面板
      panel.grid.major.y = element_line(color = "#e0e0e0", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white"),
      
      # 图例
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 13),
      legend.key.size = unit(0.8, "lines"),
      legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
      legend.margin = margin(5, 5, 5, 5),
      
      # 边距
      plot.margin = margin(10, 10, 10, 10)
    )
}

# =============================================================================
# 准备共性通路数据
# =============================================================================

common_terms_sorted <- all_enrichment_df_unique %>%
  filter(Term %in% common_terms, Group == "ALL") %>%
  arrange(PValue) %>%
  pull(Term) %>%
  unique()


top_n_heatmap <- min(22, length(common_terms_sorted))
top_terms_heatmap <- common_terms_sorted[1:top_n_heatmap]

plot_data_B <- all_enrichment_df_unique %>%
  filter(Term %in% top_terms_heatmap) %>%
  mutate(
    NegLog10P = -log10(PValue + 1e-10),
    Term = factor(Term, levels = rev(top_terms_heatmap)),
    Group = factor(Group, levels = c("ALL", "HR+", "HER2+", "TNBC"))
  )

# 绘制方案B
plot_B <- ggplot(plot_data_B, aes(x = Group, y = Term)) +
  geom_point(aes(size = FoldEnrich, fill = NegLog10P), 
             shape = 21, color = "black", stroke = 0.6, alpha = 0.9) +
  scale_size_continuous(
    range = c(3, 14), 
    name = "Fold Enrichment",
    breaks = c(2, 4, 6, 8),
    guide = guide_legend(order = 1, override.aes = list(fill = "#999999"))
  ) +
  scale_fill_gradientn(
    colors = blue_white_red_colors,
    name = expression(-log[10](italic(p))),
    guide = guide_colorbar(order = 2, barwidth = 1.2, barheight = 10)
  ) +
  scale_x_discrete(expand = c(0.15, 0.15)) +
  
  # ✅关键：y轴标签去掉上下调后缀（仅影响显示，不影响数据）
  scale_y_discrete(labels = function(x) {
    x %>%
      str_replace_all("\n\\s*\\(up-regulated\\)", "") %>%
      str_replace_all("\n\\s*\\(down-regulated\\)", "")
  }) +
  
  labs(
    # title = "Top 20 Enriched Pathways Across Breast Cancer Subtypes",
    x = NULL,
    y = NULL
  ) +
  theme_reference() +
  theme(
    axis.text.y = element_text(size = 17),
    legend.spacing.y = unit(0.8, "cm")
  )


# 保存方案B
ggsave("Figure2E.png", plot = plot_B, width = 10, height = 7, dpi = 300, bg = "white")
write_csv(plot_data_B, "Figure2E.csv")
