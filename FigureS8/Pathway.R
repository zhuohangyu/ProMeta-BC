library(dplyr)
library(ggplot2)

# 1) 读入通路映射表
all_path <- read.csv("All_pathway.csv",
                     stringsAsFactors = FALSE,
                     check.names = FALSE)

# 看一下结构（确认列名）
str(all_path)
head(all_path)

# 主要用到的列：
# - Metabolite name
# - Terms（通路名称）
# - Pathway ID（KEGG map 号，可选）
# 还有 Pathway level 1 / 2 可用来分层展示
# 2) 读入显著代谢物结果
meta_res <- read.csv("significant_metabolites.csv",
                     stringsAsFactors = FALSE,
                     check.names = FALSE)

# 统一列名，去掉多余列
colnames(meta_res)[1:5] <- c("Metabolite", "Log2FC", "P_value", "VIP_Score", "Significance")
meta_res <- meta_res[, 1:5]

head(meta_res)
## 背景：有通路注释的代谢物
bg_mets <- unique(all_path$`Metabolite name`)
N <- length(bg_mets)   # 背景大小

## 显著代谢物（可以根据需要调整阈值）
sig_mets <- meta_res %>%
  filter(P_value < 0.05,
         Significance %in% c("Up-regulated", "Down-regulated")) %>%
  pull(Metabolite) %>%
  unique()

## 显著代谢物中真正“有通路注释”的那部分
sig_in_bg <- intersect(sig_mets, bg_mets)
n <- length(sig_in_bg)  # 显著集合大小

cat("背景代谢物数 N =", N, "\n")
cat("显著代谢物数（有通路注释） n =", n, "\n")
# 为方便，先只保留我们需要的列
path_df <- all_path %>%
  select(`Metabolite name`, Terms, `Pathway ID`,
         `Pathway level 1`, `Pathway level 2`) %>%
  distinct()

# 按通路 Terms 分组，做超几何检验
enrich_res <- path_df %>%
  group_by(Terms, `Pathway ID`, `Pathway level 1`, `Pathway level 2`) %>%
  summarise(
    pathway_mets = list(unique(`Metabolite name`)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    # 该通路包含的代谢物数（在背景中）
    m_size = length(intersect(pathway_mets, bg_mets)),
    # 显著代谢物中命中该通路的数目
    k_hits = length(intersect(pathway_mets, sig_in_bg)),
    N_universe = N,
    n_sig      = n,
    # 超几何检验，右尾（至少 k_hits 个）
    p_value = ifelse(
      m_size > 0 & k_hits > 0,
      phyper(q = k_hits - 1, m = m_size, n = N_universe - m_size,
             k = n_sig, lower.tail = FALSE),
      1
    )
  ) %>%
  ungroup() %>%
  # 只保留至少命中 1 个显著代谢物的通路
  filter(k_hits > 0) %>%
  mutate(
    FDR = p.adjust(p_value, method = "fdr")
  ) %>%
  arrange(FDR)

enrich_res_out <- enrich_res %>%
  mutate(
    # 差异代谢物的命中名单：pathway_mets ∩ sig_in_bg
    hit_metabolites = sapply(pathway_mets, function(x) {
      paste(intersect(x, sig_in_bg), collapse = "; ")
    })
  ) %>%
  # 删除 list 列，保留字符版代谢物列表
  select(-pathway_mets)

# 再写出 CSV
write.csv(enrich_res_out,
          "Metabolite_Pathway_ORA_from_All_pathway.csv",
          row.names = FALSE)

top_n <- 20
plot_df <- enrich_res %>%
  head(top_n) %>%
  # 为了方便排序，把 Terms 变成 factor
  mutate(Terms = factor(Terms, levels = rev(unique(Terms))))

ggplot(plot_df,
       aes(x = Terms, y = -log10(p_value))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "",
    y = expression(-log[10](P[value])),
    title = "Pathway enrichment of significant metabolites"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

ggsave("Metabolite_Pathway_ORA_barplot_from_All_pathway.pdf",
       width = 7, height = 5)
