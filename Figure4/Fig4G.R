library(tidyverse)
library(readr)
library(forcats)

# ===== 1. 读取 MSEA 数据 =====
df1 <- read_tsv("Table 4c.MSEA Enrichment of classification_1-vs-benign.xls")
df2 <- read_tsv("Table 4c.MSEA Enrichment of classification_2-vs-benign.xls")
df3 <- read_tsv("Table 4c.MSEA Enrichment of classification_3-vs-benign.xls")
df4 <- read_tsv("Table 4c.MSEA Enrichment of malignant-vs-benign.xls")

df1$Group <- "Class1 vs Benign"
df2$Group <- "Class2 vs Benign"
df3$Group <- "Class3 vs Benign"
df4$Group <- "Malignant vs Benign"

df_all <- bind_rows(df1, df2, df3, df4)

# ===== 2. 数据处理 =====
df_all <- df_all %>%
  mutate(
    p.value = ifelse(p.value > 1, 1, p.value),
    sig = ifelse(p.value < 0.05, "*", "")
  )

# 只保留 p 值 < 0.01 的通路
df_all <- df_all %>%
  filter(p.value < 0.1)

# 如果需要，仍然要求至少出现在 3 个 Group
term_count <- df_all %>%
  group_by(Terms) %>%
  summarise(n_group = n_distinct(Group), .groups = "drop")

keep_terms <- term_count %>%
  filter(n_group >= 3) %>%
  pull(Terms)

# 按 malignant-vs-benign 的 NES 排序
malignant_order <- df_all %>%
  filter(Group == "Malignant vs Benign", Terms %in% keep_terms) %>%
  arrange(NES) %>%
  pull(Terms)

df_plot <- df_all %>%
  filter(Terms %in% keep_terms) %>%
  mutate(
    Terms = factor(Terms, levels = malignant_order),
    Group = factor(Group, levels = c("Malignant vs Benign",
                                     "Class1 vs Benign",
                                     "Class2 vs Benign",
                                     "Class3 vs Benign"))
  ) %>%
  mutate(Group = fct_recode(Group,
                            "Malignant" = "Malignant vs Benign",
                            "HR+/HER-" = "Class1 vs Benign",
                            "HER+"     = "Class2 vs Benign",
                            "TNBC"     = "Class3 vs Benign"))

# ===== 3. 绘图 =====
p <- ggplot(df_plot, aes(x = NES, y = Terms)) +
  geom_point(aes(size = Count, fill = p.value),
             shape = 21, color = "black", stroke = 0.8, alpha = 0.9) +
  geom_text(aes(label = sig),
            nudge_x = -0.8, vjust = 0.8,
            size = 8, color = "red", fontface = "bold") +
  scale_fill_gradient(
    low = "red", high = "grey90",
    limits = c(0, 0.1),                     # 限制范围到 0–0.1
    breaks = c(0, 0.01, 0.05, 0.1),         # 设置刻度
    guide = guide_colorbar(
      barwidth = 12,
      barheight = 0.8,
      title.position = "left",
      title.vjust = 1.2
    )
  ) +
  scale_size(range = c(3, 6), breaks = c(2, 6, 10), name = "Mapping count") +
  facet_wrap(~Group, scales = "free_x", nrow = 1) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 14),
    axis.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 16),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.spacing.y = unit(0.1, "cm"),
    legend.margin = margin(b = -5),
    strip.text = element_text(face = "bold", size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(t = 0, r = 15, b = 10, l = 15)
  ) +
  labs(
    x = "NES (Normalized Enrichment Score)",
    y = "MSEA Pathway",
    fill = "p-value",
    title = "MSEA Pathway Enrichment"
  ) +
  coord_cartesian(clip = "off")

ggsave("Figure4G.png", p, width = 10, height = 8, dpi = 300)
