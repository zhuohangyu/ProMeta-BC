# 清理工作空间
rm(list = ls())
gc()
#----------------------------------------------------------------------------------
# Step 1: 读取表达矩阵和临床信息
#----------------------------------------------------------------------------------
library(tidyr)
library(tibble)
library(ggtext)

clinical_data <- read.csv("Clinical.csv", stringsAsFactors = FALSE)
protein_data <- read.csv("Combined_data_filtered.csv", stringsAsFactors = FALSE, check.names = FALSE)

# 提取表达矩阵 (列4开始是表达值)
protein_matrix <- protein_data[, 4:ncol(protein_data)]
rownames(protein_matrix) <- protein_data[,1]  # 第一列是Gene name
expr <- t(protein_matrix)                     # 样本 x 基因

# 样本顺序对齐
sample_names <- rownames(expr)
clinical_data_ordered <- clinical_data[match(sample_names, clinical_data$Sample), ]

#----------------------------------------------------------------------------------
# Step 2: 定义分组信息
#----------------------------------------------------------------------------------

# 分子分型 (1=HR+/HER-, 2=HER2+, 3=TNBC)
clinical_data_ordered$Subtype <- ifelse(clinical_data_ordered$Label == "Benign", 
                                        "Benign",
                                        ifelse(clinical_data_ordered$Molecular.Classification == 1, "HR+/HER-", 
                                               ifelse(clinical_data_ordered$Molecular.Classification == 2, "HER2+", 
                                                      ifelse(clinical_data_ordered$Molecular.Classification == 3, "TNBC", NA))))

table(clinical_data_ordered$Subtype)

# 定义 Malignant vs Benign
clinical_data_ordered$Malignant <- ifelse(clinical_data_ordered$Label == "Benign", "Benign", "Malignant")

#----------------------------------------------------------------------------------
# Step 3: 差异分析 (limma)
#----------------------------------------------------------------------------------

library(limma)

expr <- t(expr)   # 转成基因 x 样本

# (a) Malignant vs Benign
group1 <- factor(clinical_data_ordered$Malignant, levels=c("Benign","Malignant"))
design1 <- model.matrix(~0+group1)
colnames(design1) <- make.names(colnames(design1))

contr1 <- makeContrasts(Malignant_vs_Benign = group1Malignant - group1Benign,  levels=design1)

fit1 <- lmFit(expr, design1)
fit1 <- contrasts.fit(fit1, contr1)
fit1 <- eBayes(fit1)


library(dplyr)

# 1. 去掉 Subtype 为空或 NA 的样本
clinical_clean <- clinical_data_ordered %>%
  filter(!is.na(Subtype) & Subtype != "")
# 2. 找共同样本ID（假设临床里是 SampleID，expr 列名是样本）
common_ids <- intersect(colnames(expr), clinical_clean$Sample)

length(common_ids)      # 看看有多少个样本保留下来
# 3. 表达矩阵只保留共同样本列
expr2 <- expr[, common_ids]

# 4. 临床数据按 expr2 的列顺序重排
clinical_clean2 <- clinical_clean %>%
  filter(Sample %in% common_ids) %>%
  arrange(match(Sample, common_ids))

# 检查顺序是否一致（应该全是 TRUE）
all(colnames(expr2) == clinical_clean2$Sample)

# 5. 构建分组：Benign, HR+/HER-, HER2+, TNBC
group2 <- factor(
  clinical_clean2$Subtype,
  levels = c("Benign", "HR+/HER-", "HER2+", "TNBC")
)

design2 <- model.matrix(~0 + group2)
colnames(design2) <- make.names(colnames(design2))

# 看维度是否匹配
dim(expr2)     # 行：蛋白，列：样本
dim(design2)   # 行：样本，列：分组

# 6. 对比矩阵
contr2 <- makeContrasts(
  HR_vs_Benign   = group2HR..HER. - group2Benign,
  HER2_vs_Benign = group2HER2.    - group2Benign,
  TNBC_vs_Benign = group2TNBC     - group2Benign,
  levels = design2
)

# 7. 拟合模型（注意现在用 expr2）
fit2 <- lmFit(expr2, design2)
fit2 <- contrasts.fit(fit2, contr2)
fit2 <- eBayes(fit2)


# 提前保存 Gene name
gene_names <- protein_data[,2]
names(gene_names) <- rownames(protein_matrix)  # 确保对应关系

res_all <- topTable(fit1, coef="Malignant_vs_Benign", number=Inf, adjust="fdr")
res_all$Gene <- rownames(res_all)
res_all$Gene.name <- gene_names[res_all$Gene]   # 保留原始Gene nam
res_all$Comparison <- "Malignant vs Benign"

# (b) 三个亚型 vs Benign
res_list <- list(
  "HR+ vs Benign"   = topTable(fit2, coef="HR_vs_Benign", number=Inf, adjust="fdr"),
  "HER2+ vs Benign" = topTable(fit2, coef="HER2_vs_Benign", number=Inf, adjust="fdr"),
  "TNBC vs Benign"  = topTable(fit2, coef="TNBC_vs_Benign", number=Inf, adjust="fdr")
)

# 添加 Gene 和 Gene.name
for (nm in names(res_list)) {
  res_list[[nm]]$Gene <- rownames(res_list[[nm]])
  res_list[[nm]]$Gene.name <- gene_names[res_list[[nm]]$Gene]
  res_list[[nm]]$Comparison <- nm
}

#----------------------------------------------------------------------------------
# Step 4: 整理结果 & 转换为 long_data
#----------------------------------------------------------------------------------

library(dplyr)
library(tidyr)

# 把 Malignant vs Benign 加入 res_list
res_list$`Malignant vs Benign` <- res_all

# 添加基因名列
for (nm in names(res_list)) {
  res_list[[nm]]$Gene <- rownames(res_list[[nm]])
}
long_data <- bind_rows(
  lapply(names(res_list), function(comp){
    df <- res_list[[comp]]
    df %>%
      mutate(
        logP = -log10(adj.P.Val),
        group = case_when(
          adj.P.Val < 0.05 & logFC > 0.5 ~ "up",
          adj.P.Val < 0.05 & logFC < -0.5 ~ "down",
          TRUE ~ "not sig"
        )
      ) %>%
      dplyr::select(Gene, Gene.name, logFC, P.Value, adj.P.Val, logP, group, Comparison)
  })
)
# 保存
write.csv(long_data, "./DEGs/long_data_with_genename.csv", row.names=FALSE)

#----------------------------------------------------------------------------------
# Step 5: 设置Comparison顺序
#----------------------------------------------------------------------------------
comparison_levels <- c("Malignant vs Benign", "HR+ vs Benign", "HER2+ vs Benign", "TNBC vs Benign")

long_data$Comparison <- factor(long_data$Comparison, levels = comparison_levels)

# 如果有统计表，也要设置顺序
# （后面 stats_summary / title_text 都基于 long_data，所以顺序会保持一致）


#----------------------------------------------------------------------------------
# Step 6: 多组别可视化
#----------------------------------------------------------------------------------

library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(UpSetR)
library(circlize)
library(cowplot)

# 自定义颜色
colors_regulation <- c("up" = "#E74C3C", "down" = "#3498DB", "not sig" = "grey70")

#----------------------------------------------------------------------------------
# 美化版多比较点图（带中间方块标签）
#----------------------------------------------------------------------------------

library(ggrepel)
library(ggnewscale)
library(scales)

# 自定义颜色
colors_regulation <- c("up" = "#E74C3C", "down" = "#3498DB", "not sig" = "grey70")
colors_comparison <- c(
  "Malignant vs Benign" = "#C0392B", 
  "HR+ vs Benign"       = "#8E44AD", 
  "HER2+ vs Benign"     = "#27AE60",
  "TNBC vs Benign"      = "#F39C12"
)

# 各组上下调数量
stats_summary <- long_data %>%
  group_by(Comparison, group) %>%
  summarise(count=n(), .groups="drop") %>%
  tidyr::pivot_wider(names_from=group, values_from=count, values_fill=0)

# 选择每组前5上调和前5下调基因
top_up_genes <- long_data %>%
  filter(group=="up") %>%
  group_by(Comparison) %>%
  arrange(desc(logFC)) %>%
  slice_head(n=5) %>%
  ungroup()

# top_down_genes <- long_data %>%
#   filter(group=="down") %>%
#   group_by(Comparison) %>%
#   arrange(logFC) %>%
#   slice_head(n=2) %>%
#   ungroup()

# top_genes <- bind_rows(top_up_genes, top_down_genes)
top_genes <- top_up_genes

# 统计标签
# title_text <- stats_summary %>%
#   mutate(
#     label = paste0(
#       "<span style='color:#E74C3C;'>Up: ", up, "</span><br>",
#       "<span style='color:#3498DB;'>Down: ", down, "</span>"
#     )
#   )
title_text1 <- stats_summary %>%
  mutate(label = paste0("Up: ", up))
title_text2 <- stats_summary %>%
  mutate(label = paste0("Down: ", down))

# 中间方块高度（放在 y=0 附近）
tile_y <- 0  
ymin <- min(long_data$logFC, na.rm=TRUE) - 1
ymax <- max(long_data$logFC, na.rm=TRUE) + 1.8
# 美化 dotplot
p1 <- ggplot() +
  
  geom_jitter(
    data = long_data %>% filter(abs(logFC) > 0.2),
    aes(x = Comparison, y = logFC, color = group),
    size = 1.6, alpha = 0.75, width = 0.35
  ) +
  
  geom_tile(
    data = title_text1,
    aes(x = Comparison, y = tile_y, fill = Comparison),
    height = 1, width = 1, alpha = 0.75, color = "white"
  ) +
  
  # geom_text(
  #   data = title_text,
  #   aes(x = Comparison, y = tile_y, label = Comparison),
  #   size = 6,                      # ← 比原来更大
  #   fontface = "bold",
  #   color = "white"
  # ) +
  geom_text(
    data = title_text1,
    aes(
      x = Comparison,
      y = tile_y,
      label = gsub(" vs Benign", "", Comparison)   # ✅删除 vs Benign
    ),
    size = 8,
    fontface = "bold",
    color = "white"
  ) +

  
  geom_text_repel(
    data = filter(top_genes, group == "up"),
    aes(x = Comparison, y = logFC, label = Gene.name),
    size = 7,                      # ← 放大
    color = colors_regulation["up"],
    fontface = "bold",
    max.overlaps = 20
  ) +
  
  geom_text_repel(
    data = filter(top_genes, group == "down"),
    aes(x = Comparison, 
        y = logFC, 
        label = Gene.name),
    size = 7,                      # ← 放大
    color = colors_regulation["down"],
    fontface = "bold",
    max.overlaps = 20
  ) +
  
  geom_hline(
    yintercept = c(-0.5, 0.5),
    linetype = "dashed",
    color = "black",
    alpha = 0.6
  ) +
  
  scale_color_manual(values = colors_regulation, name = "") +
  scale_fill_manual(values = colors_comparison, guide = "none") +
  
  # geom_richtext(
  #   data = title_text,
  #   aes(
  #     x = Comparison,
  #     y = max(long_data$logFC, na.rm = TRUE) + 1.5,
  #     label = label
  #   ),
  #   size = 7,
  #   fontface = "bold",
  #   hjust = 0.5,
  #   fill = NA,          # ✅去掉背景框
  #   label.color = NA    # ✅去掉边框
  # )+

  geom_text(
    data = title_text1,
    aes(x = Comparison,
        y = max(long_data$logFC, na.rm = TRUE) + 1.2,
        label = label),
    size = 8,
    fontface = "bold",
    hjust = 0.5
  ) +
  geom_text(
    data = title_text2,
    aes(x = Comparison,
        y = min(long_data$logFC, na.rm = TRUE) - 0.2,
        label = label),
    size = 8,
    fontface = "bold",
    hjust = 0.5
  ) +


  scale_y_continuous(
    limits = c(ymin, ymax),
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(y = "log2 Fold Change") +
  
  theme_minimal(base_size = 28) +   # ← 全局文字整体放大
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    # legend.position = "top",
    # legend.text = element_text(size = 24)   # ← 图例字体加大
  )

# p1 <- p1 +
#   theme(
#     legend.box.margin = margin(t = 0, r = 0, b = -5, l = 0),   # 压缩外边距
#     legend.spacing.y  = unit(0.01, "cm"),                       # 行间距更紧
#     legend.key.height = unit(0.2, "cm")                       # 小圆点背景高度
#   )
# legend.key.height = unit(0.18, "cm")


ggsave("Figure3B.png", p1, width = 10, height = 8, dpi = 300, bg = "white")




# 2) 火山图 (每组单独) ------------------------------------------------------
create_volcano_plot <- function(data, comparison_name) {
  plot_data <- data %>% filter(Comparison == comparison_name)
  top_up <- plot_data %>% filter(group=="up") %>% arrange(desc(logFC)) %>% head(5)
  top_down <- plot_data %>% filter(group=="down") %>% arrange(logFC) %>% head(5)
  top_genes <- bind_rows(top_up, top_down)
  
  ggplot(plot_data, aes(x=logFC, y=logP)) +
    geom_point(aes(color=group), alpha=0.7, size=1.5) +
    scale_color_manual(values=colors_regulation) +
    geom_vline(xintercept=c(-0.263,0.263), linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    geom_text_repel(data=top_genes, aes(label=Gene), size=3, max.overlaps=15) +
    labs(title=paste("Volcano Plot:", comparison_name), x="log2 Fold Change", y="-log10(P-value)") +
    theme_minimal(base_size=14)
}

comparisons <- unique(long_data$Comparison)
volcano_plots <- lapply(comparisons, function(comp) create_volcano_plot(long_data, comp))

# 保存单图和组合图
for (i in seq_along(comparisons)) {
  ggsave(paste0("./DEGs/volcano_", gsub(" ", "_", comparisons[i]), ".png"),
         volcano_plots[[i]], width=8, height=6, dpi=300)
}
combined_volcano <- plot_grid(plotlist=volcano_plots, ncol=2)
ggsave("./DEGs/combined_volcano_plots.png", combined_volcano, width=14, height=10, dpi=300)

# 3) UpSet 图 ---------------------------------------------------------------
filtered_res <- long_data %>% filter(group != "not sig")
gene_sets <- filtered_res %>% group_by(Comparison) %>%
  summarise(genes = list(unique(Gene)), .groups="drop") %>% deframe()
all_genes <- unique(unlist(gene_sets))
upset_data <- data.frame(Gene=all_genes)
for (comp in names(gene_sets)) {
  upset_data[[comp]] <- ifelse(upset_data$Gene %in% gene_sets[[comp]], 1, 0)
}
png("./DEGs/upset_plot.png", width=10, height=6, units="in", res=300)
upset(upset_data[,-1], sets=names(gene_sets), order.by="freq")
dev.off()

#----------------------------------------------------------------------------------
# Step 7: 生成 Regulation 宽表（类似 Fold_Change.csv）
#         列：Gene / Gene.name / Regulation_ALL / Regulation_HR+ / Regulation_HER+ / Regulation_TNBC
#----------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(VennDiagram)
library(grid)

# 先把 long_data 的 group 转成 up / down / unchanged
reg_long <- long_data %>%
  mutate(Reg = case_when(
    group == "up"   ~ "up",
    group == "down" ~ "down",
    TRUE            ~ "unchanged"
  )) %>%
  select(Gene, Gene.name, Comparison, Reg) %>%
  distinct()

# 把 Comparison 映射为 Regulation 列名
# Malignant vs Benign -> Regulation_ALL
# HR+ vs Benign       -> Regulation_HR+
# HER2+ vs Benign     -> Regulation_HER+
# TNBC vs Benign      -> Regulation_TNBC
reg_long <- reg_long %>%
  mutate(RegCol = case_when(
    Comparison == "Malignant vs Benign" ~ "Regulation_ALL",
    Comparison == "HR+ vs Benign"       ~ "Regulation_HR+",
    Comparison == "HER2+ vs Benign"     ~ "Regulation_HER+",
    Comparison == "TNBC vs Benign"      ~ "Regulation_TNBC",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(RegCol)) %>%
  select(Gene, Gene.name, RegCol, Reg)

# 宽表：每个 Gene 一行，4 个 Regulation 列
reg_wide <- reg_long %>%
  pivot_wider(
    names_from  = RegCol,
    values_from = Reg,
    values_fill = "unchanged"
  )

# 保存成类似 Fold_Change 的文件
if (!dir.exists("./DEGs")) dir.create("./DEGs")
write.csv(reg_wide, "./DEGs/Regulation_for_venn.csv", row.names = FALSE)

#----------------------------------------------------------------------------------
# Step 8: 基于 Regulation 宽表筛选 fc_filt，并画 HR+ / HER+ / TNBC 的 Venn 图
#----------------------------------------------------------------------------------

# fc_filt：在 ALL / HR+ / HER+ / TNBC 中至少有一个不是 unchanged
fc_filt <- reg_wide
  # filter(
  #   Regulation_ALL  != "unchanged" 
  # )

# 三个亚型的“有变化”基因集合（严格来自 fc_filt）
genes_HR <- fc_filt %>%
  filter(`Regulation_HR+` != "unchanged") %>%
  pull(Gene.name) %>%
  unique()

genes_HER <- fc_filt %>%
  filter(`Regulation_HER+` != "unchanged") %>%
  pull(Gene.name) %>%
  unique()

genes_TNBC <- fc_filt %>%
  filter(`Regulation_TNBC` != "unchanged") %>%
  pull(Gene.name) %>%
  unique()

# 也可以用 Gene 而不是 Gene.name，看你后续 Venn 图里想显示什么
# 把名字换成 Gene 即可：
# pull(Gene)

# 检查数量
length(genes_HR)
length(genes_HER)
length(genes_TNBC)

# 画 Venn 图
venn_list <- list(
  "HR+"   = genes_HR,
  "HER2+" = genes_HER,
  "TNBC"  = genes_TNBC
)

grid.newpage()
venn_plot <- venn.diagram(
  x        = venn_list,
  filename = NULL,  # 返回 grob，不直接写文件
  fill     = c("#6BAED6", "#74C476", "#FD8D3C"),
  alpha    = 0.6,
  cex      = 1.4,
  cat.cex  = 1.4,
  cat.dist = 0.06,
  lwd      = 1.2
)
grid.draw(venn_plot)

# 若需保存为 PNG：
png("FigureS3C.png", width = 2000, height = 2000, res = 300)
grid.newpage()
grid.draw(venn_plot)
dev.off()


