# ============================================================================
# 通路活性相关性分析 - 精简版
# ME01/ME09/ME10 + Cluster1/Cluster9
# ============================================================================

library(GSVA)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(dplyr)

# ---- 1. 读取数据 ----
protein_data <- read.csv("Combined_data_filtered.csv", check.names = FALSE)
rownames(protein_data) <- make.unique(protein_data$`Gene name`)
protein_expr <- as.matrix(protein_data[, 4:ncol(protein_data)])
storage.mode(protein_expr) <- "numeric"

metab_expr <- as.matrix(read.csv("Sample.csv", row.names = 1, check.names = FALSE))

# ---- 2. 读取KEGG富集结果 ----
me01 <- read.csv("ME01_KEGG.csv")[, 1:4]
me09 <- read.csv("ME09_KEGG.csv")[, 1:4]
me10 <- read.csv("ME10_KEGG.csv")[, 1:4]
colnames(me01) <- colnames(me09) <- colnames(me10) <- c("Module", "Gene", "Terms", "Pathway_ID")
protein_kegg <- rbind(me01, me09, me10)

cluster1 <- read.csv("Cluster1.csv")
cluster9 <- read.csv("Cluster9.csv")
colnames(cluster1) <- colnames(cluster9) <- c("Cluster", "Metabolite", "Terms", "Pathway_ID")
metab_kegg <- rbind(cluster1, cluster9)

# ---- 3. 构建通路基因集合 ----
# 蛋白通路
prot_sets <- protein_kegg %>%
  filter(!is.na(Terms) & Terms != "NA") %>%
  group_by(Terms) %>%
  summarise(genes = list(unique(Gene)), .groups = "drop")
prot_gene_sets <- setNames(lapply(prot_sets$genes, function(x) 
  intersect(as.character(x), rownames(protein_expr))), prot_sets$Terms)
prot_gene_sets <- prot_gene_sets[sapply(prot_gene_sets, length) >= 1]

# 代谢物通路
metab_sets <- metab_kegg %>%
  filter(!is.na(Terms) & Terms != "NA") %>%
  group_by(Terms) %>%
  summarise(metabolites = list(unique(Metabolite)), .groups = "drop")
metab_compound_sets <- setNames(lapply(metab_sets$metabolites, function(x) 
  intersect(as.character(x), rownames(metab_expr))), metab_sets$Terms)
metab_compound_sets <- metab_compound_sets[sapply(metab_compound_sets, length) >= 1]

cat("蛋白通路:", length(prot_gene_sets), "\n")
cat("代谢物通路:", length(metab_compound_sets), "\n")

# ---- 4. 共同样本 & 数据预处理 ----
common_samples <- intersect(colnames(protein_expr), colnames(metab_expr))
cat("共同样本:", length(common_samples), "\n")

prot_mat <- protein_expr[, common_samples]
prot_mat[prot_mat == 0] <- NA
# 填充NA
for(i in 1:nrow(prot_mat)) {
  na_idx <- is.na(prot_mat[i, ])
  if(any(na_idx) && !all(na_idx)) prot_mat[i, na_idx] <- mean(prot_mat[i, !na_idx])
}
prot_mat <- prot_mat[rowSums(!is.na(prot_mat)) > ncol(prot_mat) * 0.3, ]

metab_mat <- metab_expr[, common_samples]

# ---- 5. GSVA计算 ----
cat("计算GSVA...\n")
gsva_prot <- gsva(prot_mat, prot_gene_sets, method = "gsva", kcdf = "Gaussian", 
                  min.sz = 1, verbose = FALSE)
gsva_metab <- gsva(metab_mat, metab_compound_sets, method = "gsva", kcdf = "Gaussian", 
                   min.sz = 1, verbose = FALSE)

# ---- 6. Spearman相关性 ----
cor_mat <- cor(t(gsva_prot), t(gsva_metab), method = "spearman", use = "pairwise.complete.obs")

# P值计算
pval_mat <- matrix(NA, nrow(cor_mat), ncol(cor_mat), dimnames = dimnames(cor_mat))
for(i in 1:nrow(gsva_prot)) {
  for(j in 1:nrow(gsva_metab)) {
    pval_mat[i,j] <- cor.test(gsva_prot[i,], gsva_metab[j,], method = "spearman", exact = FALSE)$p.value
  }
}

cat("\n相关性矩阵:\n")
print(round(cor_mat, 3))

# ---- 6.1 多重检验校正（BH-FDR）----
# 将矩阵打平成向量做统一校正（推荐：全局控制FDR）
pval_vec <- as.vector(pval_mat)
padj_vec <- p.adjust(pval_vec, method = "BH")

padj_mat <- matrix(padj_vec,
                   nrow = nrow(pval_mat),
                   ncol = ncol(pval_mat),
                   dimnames = dimnames(pval_mat))

cat("\nFDR(BH) 校正完成。\n")


# ---- 7. 保存结果 ----
write.csv(cor_mat,  "pathway_correlation_matrix.csv")
write.csv(pval_mat, "pathway_correlation_pvalue_raw.csv")
write.csv(padj_mat, "pathway_correlation_pvalue_FDR_BH.csv")
write.csv(gsva_prot, "gsva_protein_scores.csv")
write.csv(gsva_metab, "gsva_metabolite_scores.csv")

# 详细结果表
# 详细结果表
results <- expand.grid(
  Protein = rownames(cor_mat),
  Metabolite = colnames(cor_mat),
  stringsAsFactors = FALSE
)

results$Correlation <- mapply(function(p, m) cor_mat[p, m], results$Protein, results$Metabolite)
results$P_value_raw <- mapply(function(p, m) pval_mat[p, m], results$Protein, results$Metabolite)
results$FDR_BH      <- mapply(function(p, m) padj_mat[p, m], results$Protein, results$Metabolite)

# 用 FDR 打星号（你也可以改阈值）
results$Sig <- ifelse(results$FDR_BH < 0.001, "***",
                      ifelse(results$FDR_BH < 0.01, "**",
                             ifelse(results$FDR_BH < 0.05, "*", "")))

results <- results[order(-abs(results$Correlation)), ]
write.csv(results, "pathway_correlation_results_FDR_BH.csv", row.names = FALSE)


# ---- 8. 热图可视化 ----
my_col <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

# 显著性标记
sig_mat <- matrix(paste0(round(cor_mat, 2), 
                         ifelse(pval_mat < 0.001, "***", 
                                ifelse(pval_mat < 0.01, "**", 
                                       ifelse(pval_mat < 0.05, "*", "")))),
                  nrow = nrow(cor_mat), dimnames = dimnames(cor_mat))

pdf("pathway_correlation_heatmap.pdf", width = 12, height = 9)
pheatmap(cor_mat, color = my_col, breaks = seq(-1, 1, length.out = 101),
         display_numbers = sig_mat, number_color = "black", fontsize_number = 9,
         main = "Protein-Metabolite Pathway Correlation (GSVA + Spearman)",
         angle_col = 45, border_color = "grey80", cellwidth = 55, cellheight = 28)
dev.off()

# ggplot版本
plot_df <- melt(cor_mat, varnames = c("Protein", "Metabolite"), value.name = "Cor")
plot_df$Pval <- melt(pval_mat)$value
plot_df$Label <- paste0(round(plot_df$Cor, 2), 
                        ifelse(plot_df$Pval < 0.001, "***", 
                               ifelse(plot_df$Pval < 0.01, "**", 
                                      ifelse(plot_df$Pval < 0.05, "*", ""))))

p <- ggplot(plot_df, aes(Metabolite, Protein, fill = Cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Label), size = 3.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       limits = c(-1, 1), name = "Spearman r") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Pathway Activity Correlation\n* p<0.05, ** p<0.01, *** p<0.001",
       x = "Metabolite Pathways", y = "Protein Pathways")

ggsave("pathway_correlation_ggplot.pdf", p, width = 13, height = 9)
ggsave("pathway_correlation_ggplot.png", p, width = 13, height = 9, dpi = 300)

cat("\n分析完成！\n")
cat("显著相关 (p<0.05):\n")
print(results[results$P_value < 0.05, ], row.names = FALSE)