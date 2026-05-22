# 清理工作空间
rm(list = ls())
gc()

# 安装和加载必需的包
required_packages <- c("WGCNA", "flashClust", "ggplot2", "dplyr", "ComplexHeatmap", 
                       "circlize", "ggpubr", "tibble", "corrplot")

for (pkg in required_packages) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 设置WGCNA参数
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# =============================================================================
# 1. 数据读取和预处理
# =============================================================================

# 读取蛋白表达数据
protein_data <- read.csv("Combined_data_filtered.csv", header = TRUE, row.names = 1, check.names = FALSE)
print(paste("原始蛋白表达数据维度:", dim(protein_data)))

# 读取临床数据
clinical_data <- read.csv("Clinical.csv", header = TRUE)
print(paste("临床数据维度:", dim(clinical_data)))

# 查看Label分布
print("Label分布:")
print(table(clinical_data$Label))

# 筛选Malignant样本用于网络构建
malignant_samples <- clinical_data[clinical_data$Label == "Malignant", "Sample"]
print(paste("恶性样本数量:", length(malignant_samples)))

# 提取蛋白信息列
protein_info <- protein_data[, 1:2]  # Gene name 和 Protein description
colnames(protein_info) <- c("Gene_name", "Protein_description")

# 提取表达数据（跳过前两列蛋白信息）
expression_data <- protein_data[, -(1:2)]

# 筛选恶性样本的表达数据用于网络构建
sample_columns_malignant <- intersect(colnames(expression_data), malignant_samples)
print(paste("匹配的恶性样本数量:", length(sample_columns_malignant)))

malignant_expression <- expression_data[, sample_columns_malignant]
print(paste("筛选后的恶性样本数据维度:", dim(malignant_expression)))

# 转置数据（WGCNA需要样本为行，基因为列）
datExpr <- as.data.frame(t(malignant_expression))
print(paste("转置后数据维度:", dim(datExpr)))

# 为列添加蛋白名称
colnames(datExpr) <- rownames(protein_info)

# =============================================================================
# 2. 数据质量控制
# =============================================================================

# 检查缺失值
print("缺失值统计:")
print(sum(is.na(datExpr)))

# 如果有缺失值，可以选择删除或填充
if (sum(is.na(datExpr)) > 0) {
  # 删除含有过多缺失值的基因
  goodGenes <- goodGenes(datExpr, verbose = 3)
  datExpr <- datExpr[, goodGenes]
  protein_info <- protein_info[goodGenes, ]
}

# 检查离群样本
gsg <- goodSamplesGenes(datExpr, verbose = 3)
print(paste("好的样本:", gsg$allOK))

if (!gsg$allOK) {
  # 移除有问题的基因和样本
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(paste("移除基因:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
  }
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(paste("移除样本:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
  }
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  protein_info <- protein_info[gsg$goodGenes, ]
}

print(paste("质控后数据维度:", dim(datExpr)))

# 样本聚类检测离群样本
sampleTree <- hclust(dist(datExpr), method = "average")
pdf("./WGCNA/01_Sample_Clustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# =============================================================================
# 3. 软阈值选择
# =============================================================================

# 选择软阈值范围
powers <- c(c(1:10), seq(from = 12, to = 24, by = 2))

# 计算网络拓扑指标
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 绘制软阈值选择图
png("./WGCNA/02_Soft_Threshold_Selection.png", 
    width = 3000, height = 1800, res = 300)  # 高清输出
par(mfrow = c(1, 2))
cex1 <- 0.85

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     labels = powers, cex = cex1, col = "red")

dev.off()


# 选择最佳软阈值
if (is.na(sft$powerEstimate)) {
  # 手动选择第一个R^2大于0.85的power值
  fit_indices <- sft$fitIndices
  powers_values <- fit_indices[, 1]
  r_squared <- -sign(fit_indices[, 3]) * fit_indices[, 2]
  softPower <- powers_values[which(r_squared > 0.85)[1]]
  if (is.na(softPower)) {
    softPower <- 6  # 如果仍然没有找到合适的值，使用默认值6
  }
  cat("Warning: 使用第一个R^2>0.85的软阈值:", softPower, "\n")
} else {
  softPower <- sft$powerEstimate
}
cat("选择的软阈值:", softPower, "\n")
cor <- WGCNA::cor
# =============================================================================
# 4. 一步法网络构建和模块识别
# =============================================================================

# 使用blockwiseModules进行网络构建
net <- blockwiseModules(datExpr, power = softPower,
                        TOMType = "unsigned", 
                        corType = "pearson",    # 强制使用Pearson相关性
                        minModuleSize = 5,        # 最小模块大小
                        reassignThreshold = 0, 
                        mergeCutHeight = 0.001,    # 很低的合并阈值
                        numericLabels = TRUE, 
                        pamRespectsDendro = FALSE,
                        deepSplit = 4,            # 最深分割
                        saveTOMs = FALSE,
                        maxBlockSize = 5000,
                        verbose = 1)

cor <- stats::cor
# 查看模块数量
print(paste("识别到的模块数量:", length(unique(net$colors))))
print("模块大小分布:")
print(table(net$colors))

# =============================================================================
# 5. 模块可视化
# =============================================================================

# 转换数字标签为颜色标签
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# 创建模块标签映射
unique_labels <- unique(net$colors)
color_mapping <- labels2colors(unique_labels)
label_color_df <- data.frame(ModuleLabel = unique_labels, ModuleColor = color_mapping)

# 创建蛋白-模块映射表
protein_color_df <- data.frame(Protein = colnames(datExpr), 
                               ModuleLabel = paste("ME", sprintf("%02d", net$colors), sep=""), 
                               ModuleColor = paste("ME", moduleColors, sep=""))
write.csv(protein_color_df, "./WGCNA/protein_module_colors_with_numbers.csv", row.names = FALSE)

# 绘制基因聚类树和模块颜色
# 绘制基因聚类树和模块颜色（高清PNG版本）
png("./WGCNA/03_Gene_Dendrogram_and_Modules.png",
    width = 3000, height = 2250, res = 300)  # 推荐分辨率与比例

plotDendroAndColors(
  geneTree,
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

dev.off()


# =============================================================================
# 6. 【新增】模块间相似性分析
# =============================================================================

# 过滤掉grey模块
filteredGenes <- moduleColors != "grey"
datExpr_filtered <- datExpr[, filteredGenes]
moduleColors_filtered <- moduleColors[filteredGenes]

# 计算模块特征基因
MEs0 <- moduleEigengenes(datExpr_filtered, moduleColors_filtered)$eigengenes
MEs <- orderMEs(MEs0)

# 计算模块间相关性矩阵
moduleTraitCor_inter <- cor(MEs, use = "pairwise.complete.obs")

# 更新行名和列名为标准化的模块标签
ME_names <- names(MEs)
module_color_to_name <- data.frame(
  ModuleColor = sub("^ME", "", ME_names),
  ModuleName = ME_names
)

new_label_color_df <- merge(module_color_to_name, label_color_df, by = "ModuleColor", all.x = TRUE)
original_order <- match(module_color_to_name$ModuleColor, new_label_color_df$ModuleColor)
new_label_color_df1 <- new_label_color_df[original_order, ]
new_label_color_df1$ModuleLabel <- paste("ME", sprintf("%02d", new_label_color_df1$ModuleLabel), sep = "")

rownames(moduleTraitCor_inter) <- new_label_color_df1$ModuleLabel[match(rownames(moduleTraitCor_inter), new_label_color_df1$ModuleName)]
colnames(moduleTraitCor_inter) <- new_label_color_df1$ModuleLabel[match(colnames(moduleTraitCor_inter), new_label_color_df1$ModuleName)]

# 按标签排序
moduleTraitCor_inter <- moduleTraitCor_inter[order(rownames(moduleTraitCor_inter)), order(colnames(moduleTraitCor_inter))]

# 使用ComplexHeatmap绘制模块间相似性热图
moduleTree <- hclust(as.dist(1 - moduleTraitCor_inter), method = "average")

# 绘制模块相似性热图（高清PNG）
png("./WGCNA/05_Module_Similarity_Heatmap.png",
    width = 3000, height = 2400, res = 300)  # 高清分辨率适用于投稿

Heatmap(
  moduleTraitCor_inter, 
  name = "Pearson \nCorrelation \nCoefficient", 
  col = colorRamp2(c(-1, -0.5, 0, 0.5, 1),
                   c("#1976D2", "#76C7FF", "#FFFFFF", "#FF987A", "#E53935")),
  cluster_rows = moduleTree, 
  cluster_columns = moduleTree, 
  row_names_side = "left", 
  column_names_side = "top", 
  show_row_dend = TRUE, 
  show_column_dend = TRUE,  # 显示行列聚类树
  heatmap_legend_param = list(
    title = "Pearson \nCorrelation \nCoefficient",
    at = c(-1, 0, 1)
  ),
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (abs(moduleTraitCor_inter[i, j]) > 0.75) grid.text("***", x, y, gp = gpar(fontsize = 10))
    else if (abs(moduleTraitCor_inter[i, j]) > 0.5) grid.text("**", x, y, gp = gpar(fontsize = 10))
    else if (abs(moduleTraitCor_inter[i, j]) > 0.3) grid.text("*", x, y, gp = gpar(fontsize = 10))
  }
)

dev.off()

# 保存模块间相关性矩阵
write.csv(moduleTraitCor_inter, "./WGCNA/module_similarity_matrix.csv", row.names = TRUE)

# =============================================================================
# 7. 【新增】模块与良性/恶性肿瘤类型的相关性分析
# =============================================================================

# 准备包含所有样本（良性+恶性）的数据进行表型相关性分析
all_samples <- clinical_data$Sample
sample_columns_all <- intersect(colnames(expression_data), all_samples)
print(paste("所有可用样本数量:", length(sample_columns_all)))

# 提取所有样本的表达数据
all_expression <- expression_data[, sample_columns_all]
datExpr_all <- as.data.frame(t(all_expression))
colnames(datExpr_all) <- rownames(protein_info)

# 质量控制（针对所有样本）
gsg_all <- goodSamplesGenes(datExpr_all, verbose = 3)
if (!gsg_all$allOK) {
  datExpr_all <- datExpr_all[gsg_all$goodSamples, gsg_all$goodGenes]
}

# 筛选与网络构建时相同的基因
common_genes <- intersect(colnames(datExpr_filtered), colnames(datExpr_all))
datExpr_all_filtered <- datExpr_all[, common_genes]
datExpr_filtered_common <- datExpr_filtered[, common_genes]
moduleColors_common <- moduleColors_filtered[match(common_genes, colnames(datExpr_filtered))]

# 使用恶性样本构建的网络，计算所有样本的模块特征基因
MEs_all <- moduleEigengenes(datExpr_all_filtered, moduleColors_common)$eigengenes
MEs_all <- orderMEs(MEs_all)

# 准备临床表型数据
datTraits_all <- clinical_data[clinical_data$Sample %in% rownames(datExpr_all_filtered), ]
datTraits_all <- datTraits_all[match(rownames(datExpr_all_filtered), datTraits_all$Sample), ]

# 创建设计矩阵（良性 vs 恶性）
datTraits_all$Label <- factor(datTraits_all$Label, levels = c("Benign", "Malignant"))
design <- model.matrix(~0 + datTraits_all$Label)
colnames(design) <- c("Benign", "Malignant")

# 计算模块与表型的相关性
nGenes_all <- ncol(datExpr_all_filtered)
nSamples_all <- nrow(datExpr_all_filtered)

moduleTraitCor <- cor(MEs_all, design, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples_all)

# 更新行名为标准化的模块标签
ME_names_all <- names(MEs_all)
module_color_to_name_all <- data.frame(
  ModuleColor = sub("^ME", "", ME_names_all),
  ModuleName = ME_names_all
)

new_label_color_df_all <- merge(module_color_to_name_all, label_color_df, by = "ModuleColor", all.x = TRUE)
original_order_all <- match(module_color_to_name_all$ModuleColor, new_label_color_df_all$ModuleColor)
new_label_color_df_all1 <- new_label_color_df_all[original_order_all, ]
new_label_color_df_all1$ModuleLabel <- paste("ME", sprintf("%02d", new_label_color_df_all1$ModuleLabel), sep = "")

rownames(moduleTraitCor) <- new_label_color_df_all1$ModuleLabel[match(rownames(moduleTraitCor), new_label_color_df_all1$ModuleName)]
rownames(moduleTraitPvalue) <- new_label_color_df_all1$ModuleLabel[match(rownames(moduleTraitPvalue), new_label_color_df_all1$ModuleName)]

# 按行名排序
moduleTraitCor <- moduleTraitCor[order(rownames(moduleTraitCor)), , drop = FALSE]
moduleTraitPvalue <- moduleTraitPvalue[order(rownames(moduleTraitPvalue)), , drop = FALSE]

# 转置矩阵（性状在行，模块在列）
moduleTraitCor <- t(moduleTraitCor)
moduleTraitPvalue <- t(moduleTraitPvalue)

# 绘制模块-表型相关性热图
# 模块-性状相关性图（ellipse + 手动画色条）导出为PNG
png("./WGCNA/06_Module_Trait_Correlation_Plot.png",
    width = 3000, height = 800, res = 300)  # 对应原来 15×3 的比例

layout(matrix(c(1, 2), nrow = 2), heights = c(1.5, 1))
par(mar = c(2, 2, 2, 10))

# 绘制相关性热图
corrplot::corrplot(
  moduleTraitCor, 
  method = "ellipse",
  col = colorRampPalette(c("#1976D2", "#76C7FF", "#FFFFFF", "#FF987A", "#E53935"))(200),
  type = "full",
  addCoef.col = NULL,
  tl.col = "black",
  tl.srt = 45,
  cl.pos = "n",                 # 不用默认legend
  p.mat = moduleTraitPvalue,
  sig.level = -1,
  insig = "n"
)

# 添加显著性标记
nRows <- nrow(moduleTraitCor)
nCols <- ncol(moduleTraitCor)
for (i in 1:nRows) {
  for (j in 1:nCols) {
    y_pos <- nRows - i + 1  # corrplot 的 y 轴是反的
    if (moduleTraitPvalue[i, j] < 0.001) {
      text(j, y_pos, "***", cex = 1.5, col = "black")
    } else if (moduleTraitPvalue[i, j] < 0.01) {
      text(j, y_pos, "**", cex = 1.5, col = "black")
    } else if (moduleTraitPvalue[i, j] < 0.05) {
      text(j, y_pos, "*", cex = 1.5, col = "black")
    }
  }
}

# 手动添加颜色图例
par(mar = c(2, 2, 2, 2))
color.bar <- function(lut, min, max, nticks = 5,
                      title = "Pearson Correlation Coefficient") {
  scale = (length(lut) - 1) / (max - min)
  plot(c(0, 1), c(0, 1), type = 'n', bty = 'n',
       xaxt = 'n', yaxt = 'n',
       xlab = '', ylab = '', main = title)
  
  x_start <- 0.35
  x_end   <- 0.65
  
  axis(1,
       at = seq(x_start, x_end, length.out = nticks),
       labels = round(seq(min, max, length.out = nticks), 2),
       las = 1, line = -0.5, cex.axis = 0.8)
  
  for (i in 1:(length(lut) - 1)) {
    x = x_start + (i - 1) * (x_end - x_start) / (length(lut) - 1)
    rect(x, 0.3,
         x + (x_end - x_start) / length(lut),
         0.7,
         col = lut[i], border = NA)
  }
}

color.bar(
  colorRampPalette(c("#1976D2", "#76C7FF", "#FFFFFF", "#FF987A", "#E53935"))(200),
  min = -1, max = 1
)

dev.off()


# =============================================================================
# 8. 【新增】模块得分分布可视化和统计检验
# =============================================================================

# 设置颜色
color.bin <- c("Malignant" = "#00599F", "Benign" = "#d80700")

# 准备绘图数据
plot.data <- MEs_all
module_colors <- sub("ME", "", colnames(plot.data))
new_module_labels <- vector("integer", length(module_colors))

for (i in seq_along(module_colors)) {
  color <- module_colors[i]
  new_module_labels[i] <- label_color_df$ModuleLabel[match(color, label_color_df$ModuleColor)]
}
new_module_labels <- paste("ME", sprintf("%02d", new_module_labels), sep = "")
colnames(plot.data) <- new_module_labels
sorted_colnames <- sort(colnames(plot.data))
plot.data <- plot.data[, sorted_colnames]

# 保存模块特征基因矩阵
write.csv(plot.data, file = "./WGCNA/ModuleEigengenes_All_Samples.csv", row.names = TRUE)

# 进行统计检验
plot.info <- NULL
module.name <- colnames(plot.data)
plot.stat <- data.frame(Module = module.name,
                        Wilcox.P = NA,
                        Wilcox.Padj = NA)

# Wilcoxon检验
for (i in 1:ncol(plot.data)) {
  sub <- data.frame(Sample = rownames(plot.data),
                    Score = as.numeric(plot.data[, i]),
                    ScoreScale = scale(as.numeric(plot.data[, i])),
                    Module = colnames(plot.data)[i],
                    Type = datTraits_all$Label)
  plot.stat$Wilcox.P[i] <- wilcox.test(sub$ScoreScale ~ sub$Type)$p.value
  plot.info <- rbind(plot.info, sub)
}

# 调整p值
plot.stat <- plot.stat %>% mutate(Wilcox.Padj = p.adjust(Wilcox.P, method = "fdr"))

# 确保Type为因子类型
plot.info <- plot.info %>% mutate(Type = factor(as.character(Type), levels = c("Malignant", "Benign")))

# 计算统计检验结果用于绘图
stat_test <- compare_means(ScoreScale ~ Type, data = plot.info, method = "wilcox.test", group.by = "Module")

# 绘制模块得分分布图
p <- ggbarplot(plot.info, x = "Module", y = "ScoreScale",
               color = "Type", fill = "Type",
               palette = color.bin, width = 0.5, size = 0,
               add = c("mean_se"), add.params = list(width = 0.5),
               order = module.name,
               position = position_dodge(0.6),
               xlab = "", ylab = "Module Score of Protein")

p <- p + theme_minimal()
p <- p + theme(
  axis.text.x = element_text(size = 20, color = "black", angle = 90, hjust = 1, vjust = 0.5),
  axis.text.y = element_text(size = 20, color = "black"),
  axis.title.y = element_text(size = 22),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 22)
)
p <- p + geom_hline(yintercept = 0, color = "black")
p <- p + geom_text(data = stat_test, aes(x = Module, y = 0.8, label = p.format), angle = 90, size = 7)

ggsave("./WGCNA/07_Module_Tumor_Type_Comparison.pdf", p, width = 30, height = 10)

# =============================================================================
# 9. 【新增】合并相关性和p值结果
# =============================================================================

# 将模块-表型相关性数据转换为数据框并合并
moduleTraitCor_df <- as.data.frame(moduleTraitCor)
moduleTraitPvalue_df <- as.data.frame(moduleTraitPvalue)

colnames(moduleTraitCor_df) <- paste0(colnames(moduleTraitCor_df), "_cor")
colnames(moduleTraitPvalue_df) <- paste0(colnames(moduleTraitPvalue_df), "_pvalue")

combined_df <- moduleTraitCor_df %>%
  rownames_to_column(var = "Trait") %>%
  left_join(moduleTraitPvalue_df %>% rownames_to_column(var = "Trait"), by = "Trait")

write.csv(combined_df, "./WGCNA/combined_module_trait_correlation.csv", row.names = FALSE)

# =============================================================================
# 10. 【新增】输出各模块的KME值
# =============================================================================

# 计算每个基因与其所在模块的KME值
file.remove('./WGCNA/All_Gene_KME.txt')
for(module in substring(colnames(MEs), 3)) {
  if(module == "grey") next
  ME <- as.data.frame(MEs[, paste("ME", module, sep="")])
  colnames(ME) <- module
  datModExpr <- datExpr_filtered[, moduleColors_filtered == module]
  if(ncol(datModExpr) > 0) {  # 确保模块不为空
    datKME <- signedKME(datModExpr, ME)
    datKME <- cbind(datKME, rep(module, nrow(datKME)))
    write.table(datKME, quote = FALSE, row.names = TRUE, append = TRUE, 
                file = "./WGCNA/All_Gene_KME.txt", col.names = FALSE)
  }
}

# =============================================================================
# 11. 导出网络到Cytoscape
# =============================================================================

# 计算拓扑重叠矩阵
TOM <- TOMsimilarityFromExpr(datExpr_filtered, power = softPower)
probes <- colnames(datExpr_filtered)
dimnames(TOM) <- list(probes, probes)

# 导出网络文件
cyt <- exportNetworkToCytoscape(TOM,
                                edgeFile = "./WGCNA/CytoscapeInput-edges-malignant.txt",
                                nodeFile = "./WGCNA/CytoscapeInput-nodes-malignant.txt",
                                weighted = TRUE,
                                threshold = 0.005,  # 可以根据需要调整阈值
                                nodeNames = probes,
                                nodeAttr = moduleColors_filtered)

# =============================================================================
# 12. 导出最终结果
# =============================================================================

# 导出模块信息
module_info <- data.frame(
  Protein_accession = rownames(protein_info)[match(colnames(datExpr_filtered), rownames(protein_info))],
  Gene_name = protein_info$Gene_name[match(colnames(datExpr_filtered), rownames(protein_info))],
  Protein_description = protein_info$Protein_description[match(colnames(datExpr_filtered), rownames(protein_info))],
  Module_Color = moduleColors_filtered,
  Module_Label = paste("ME", sprintf("%02d", net$colors[filteredGenes]), sep=""),
  stringsAsFactors = FALSE
)

write.csv(module_info, "./WGCNA/Module_Assignment_Enhanced.csv", row.names = FALSE)

# 导出统计检验结果
write.csv(plot.stat, "./WGCNA/Module_Tumor_Type_Statistics.csv", row.names = FALSE)

# =============================================================================
# 14. 【新增】额外分析：模块富集性状分析
# =============================================================================

cat("\n=== 模块与表型相关性详细分析 ===\n")

# 识别与恶性肿瘤显著相关的模块
malignant_cor <- moduleTraitCor["Malignant", ]
malignant_pval <- moduleTraitPvalue["Malignant", ]

significant_malignant <- names(malignant_cor)[malignant_pval < 0.05]
if(length(significant_malignant) > 0) {
  cat("与恶性肿瘤显著相关的模块:\n")
  for(mod in significant_malignant) {
    cat(sprintf("%s: r=%.3f, p=%.3f\n", mod, malignant_cor[mod], malignant_pval[mod]))
  }
} else {
  cat("没有发现与恶性肿瘤显著相关的模块\n")
}

# 识别与良性肿瘤显著相关的模块
benign_cor <- moduleTraitCor["Benign", ]
benign_pval <- moduleTraitPvalue["Benign", ]

significant_benign <- names(benign_cor)[benign_pval < 0.05]
if(length(significant_benign) > 0) {
  cat("\n与良性肿瘤显著相关的模块:\n")
  for(mod in significant_benign) {
    cat(sprintf("%s: r=%.3f, p=%.3f\n", mod, benign_cor[mod], benign_pval[mod]))
  }
} else {
  cat("没有发现与良性肿瘤显著相关的模块\n")
}

# 保存详细的相关性分析结果
detailed_results <- data.frame(
  Module = colnames(moduleTraitCor),
  Malignant_Correlation = malignant_cor,
  Malignant_Pvalue = malignant_pval,
  Benign_Correlation = benign_cor,
  Benign_Pvalue = benign_pval,
  Malignant_Significant = malignant_pval < 0.05,
  Benign_Significant = benign_pval < 0.05
)

write.csv(detailed_results, "./WGCNA/Detailed_Module_Trait_Analysis.csv", row.names = FALSE)

# =============================================================================
# 15. 【新增】为整合分析准备额外的输出变量
# =============================================================================

cat("\n=== 准备WGCNA整合分析所需的变量 ===\n")

# 1. 创建 wgcna_modules - WGCNA模块结果
# ================================================================
wgcna_modules <- data.frame(
  Protein = colnames(datExpr_filtered),
  Module_Color = moduleColors_filtered,
  Module_Label = paste("ME", sprintf("%02d", net$colors[filteredGenes]), sep=""),
  stringsAsFactors = FALSE
)

# 添加蛋白质注释信息
protein_annotations <- protein_info[match(wgcna_modules$Protein, rownames(protein_info)), ]
wgcna_modules$Gene_name <- protein_annotations$Gene_name
wgcna_modules$Protein_description <- protein_annotations$Protein_description

# 添加模块大小信息
module_sizes <- table(wgcna_modules$Module_Color)
wgcna_modules$Module_Size <- module_sizes[wgcna_modules$Module_Color]

kme_raw <- read.table(
  "./WGCNA/All_Gene_KME.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

# 第 10 步里写入顺序是：行名(Protein) + kME值 + 模块标签(module)
# 所以现在三列分别是 Protein, kME, Module_Color
colnames(kme_raw) <- c("Protein", "kME", "Module_Color")

# 合并 kME 到 wgcna_modules
wgcna_modules <- wgcna_modules %>%
  left_join(kme_raw[, c("Protein", "kME")], by = "Protein")

# 初始化 hub 标记
wgcna_modules$Is_Hub <- FALSE


# 为每个模块标记top 10%的hub基因（基于 kME）
if ("kME" %in% colnames(wgcna_modules)) {
  wgcna_modules$Is_Hub <- FALSE
  
  for (module in unique(wgcna_modules$Module_Color)) {
    if (module != "grey") {  # 排除grey模块
      module_proteins <- wgcna_modules[wgcna_modules$Module_Color == module, ]
      
      if (nrow(module_proteins) > 0 && !all(is.na(module_proteins$kME))) {
        # 取 top 10% 或至少 1 个
        n_hubs <- max(1, ceiling(nrow(module_proteins) * 0.1))
        
        top_idx <- order(module_proteins$kME,
                         decreasing = TRUE,
                         na.last = TRUE)[1:n_hubs]
        
        # 找回全表中的行号
        global_idx <- which(wgcna_modules$Module_Color == module)[top_idx]
        wgcna_modules$Is_Hub[global_idx] <- TRUE
      }
    }
  }
}

module_proteins <- subset(wgcna_modules, Module_Color == "turquoise")  # 举例
head(module_proteins$kME)
table(module_proteins$Is_Hub)
# =============================================================================
# 18. 关键模块的 hub 蛋白识别与导出
# =============================================================================

cat("\n=== 关键模块的 hub 蛋白识别 ===\n")

# 1）定义关键模块（与 Malignant 或 Benign 显著相关）
key_modules <- unique(c(significant_malignant, significant_benign))
key_modules <- key_modules[!is.na(key_modules)]

if (length(key_modules) == 0) {
  cat("当前没有与 Malignant 或 Benign 显著相关的模块，无法定义关键模块。\n")
} else {
  cat("关键模块（与 Malignant 或 Benign 显著相关）:\n")
  print(key_modules)
  
  # 这些是 "ME01" 这样的标签，你的 Module_Label 也是这个格式
  # 先做一个 Module_Label -> Module_Color 的映射
  label_to_color <- wgcna_modules %>%
    dplyr::distinct(Module_Label, Module_Color)
  
  key_module_colors <- label_to_color$Module_Color[
    match(key_modules, label_to_color$Module_Label)
  ]
  key_module_colors <- key_module_colors[!is.na(key_module_colors)]
  
  # 2）从 wgcna_modules 中筛选关键模块里的 hub 蛋白
  key_module_hubs <- wgcna_modules %>%
    dplyr::filter(Module_Color %in% key_module_colors,
                  Is_Hub == TRUE) %>%
    dplyr::arrange(Module_Label, dplyr::desc(kME))
  
  # 3）导出详细列表：每一行是一个 hub 蛋白
  write.csv(key_module_hubs,
            "./WGCNA/Key_Module_Hub_Proteins.csv",
            row.names = FALSE)
  
  cat("✓ 已将关键模块 hub 蛋白列表导出到 ./WGCNA/Key_Module_Hub_Proteins.csv\n")
}

# # 计算每个蛋白质与其模块的连通性 (kME)
# if (exists("MEs") && ncol(MEs) > 0) {
#   # 确保MEs和datExpr_filtered样本顺序一致
#   common_samples <- intersect(rownames(MEs), rownames(datExpr_filtered))
#   MEs_aligned <- MEs[common_samples, ]
#   datExpr_aligned <- datExpr_filtered[common_samples, ]
#   
#   # 计算kME值
#   kME_values <- signedKME(datExpr_aligned, MEs_aligned)
#   
#   # 为每个蛋白质添加其所在模块的kME值
#   wgcna_modules$kME <- NA
#   for (i in 1:nrow(wgcna_modules)) {
#     protein <- wgcna_modules$Protein[i]
#     module_ME <- paste0("ME", wgcna_modules$Module_Color[i])
#     
#     if (protein %in% rownames(kME_values) && module_ME %in% colnames(kME_values)) {
#       wgcna_modules$kME[i] <- kME_values[protein, module_ME]
#     }
#   }
# }
# 
# # 添加hub基因标记（基于kME值）
# if ("kME" %in% colnames(wgcna_modules)) {
#   # 为每个模块标记top 10%的hub基因
#   wgcna_modules$Is_Hub <- FALSE
#   for (module in unique(wgcna_modules$Module_Color)) {
#     if (module != "grey") {  # 排除grey模块
#       module_proteins <- wgcna_modules[wgcna_modules$Module_Color == module, ]
#       if (nrow(module_proteins) > 0 && !all(is.na(module_proteins$kME))) {
#         # 取top 10%或至少1个作为hub
#         n_hubs <- max(1, ceiling(nrow(module_proteins) * 0.1))
#         top_proteins <- order(module_proteins$kME, decreasing = TRUE, na.last = TRUE)[1:n_hubs]
#         wgcna_modules$Is_Hub[wgcna_modules$Module_Color == module][top_proteins] <- TRUE
#       }
#     }
#   }
# }

cat("✓ wgcna_modules 创建完成，包含", nrow(wgcna_modules), "个蛋白质，", length(unique(wgcna_modules$Module_Color)), "个模块\n")

# 2. 更新 significant_results - 基于WGCNA的差异分析
# ================================================================
# 创建基于模块-表型关联的显著性结果
wgcna_significant_results <- data.frame(
  Feature = colnames(moduleTraitCor),
  Feature_Type = "WGCNA_Module",
  stringsAsFactors = FALSE
)

# 添加与恶性肿瘤的关联
if ("Malignant" %in% rownames(moduleTraitCor)) {
  wgcna_significant_results$Malignant_Correlation <- moduleTraitCor["Malignant", ]
  wgcna_significant_results$Malignant_Pvalue <- moduleTraitPvalue["Malignant", ]
  wgcna_significant_results$Malignant_Significant <- moduleTraitPvalue["Malignant", ] < 0.05
}

# 添加与良性肿瘤的关联  
if ("Benign" %in% rownames(moduleTraitCor)) {
  wgcna_significant_results$Benign_Correlation <- moduleTraitCor["Benign", ]
  wgcna_significant_results$Benign_Pvalue <- moduleTraitPvalue["Benign", ]
  wgcna_significant_results$Benign_Significant <- moduleTraitPvalue["Benign", ] < 0.05
}

# 创建综合显著性标记
sig_cols <- grep("_Significant$", colnames(wgcna_significant_results), value = TRUE)
if (length(sig_cols) > 0) {
  wgcna_significant_results$Any_Significant <- rowSums(wgcna_significant_results[, sig_cols], na.rm = TRUE) > 0
}

# 添加模块大小和hub基因数量
for (i in 1:nrow(wgcna_significant_results)) {
  module_name <- wgcna_significant_results$Feature[i]
  # 从模块标签获取颜色
  module_color <- sub("^ME", "", module_name)
  if (module_color %in% wgcna_modules$Module_Color) {
    module_data <- wgcna_modules[wgcna_modules$Module_Color == module_color, ]
    wgcna_significant_results$Module_Size[i] <- nrow(module_data)
    wgcna_significant_results$Hub_Count[i] <- sum(module_data$Is_Hub, na.rm = TRUE)
  }
}

cat("✓ wgcna_significant_results 创建完成，包含", nrow(wgcna_significant_results), "个模块\n")

# 3. 更新 correlation_matrix - 模块间相关性
# ================================================================
# 使用之前计算的模块间相关性矩阵
wgcna_correlation_matrix <- moduleTraitCor_inter

cat("✓ wgcna_correlation_matrix 更新完成，", nrow(wgcna_correlation_matrix), "x", ncol(wgcna_correlation_matrix), "矩阵\n")

# 4. 更新 p_adjusted - 基于WGCNA的调整p值
# ================================================================
# 创建调整后的p值矩阵
wgcna_p_adjusted <- moduleTraitPvalue

# 应用FDR校正
wgcna_p_adjusted_flat <- as.vector(wgcna_p_adjusted)
wgcna_p_adjusted_flat <- p.adjust(wgcna_p_adjusted_flat, method = "fdr")
wgcna_p_adjusted <- matrix(wgcna_p_adjusted_flat, 
                           nrow = nrow(moduleTraitPvalue), 
                           ncol = ncol(moduleTraitPvalue))
rownames(wgcna_p_adjusted) <- rownames(moduleTraitPvalue)
colnames(wgcna_p_adjusted) <- colnames(moduleTraitPvalue)

cat("✓ wgcna_p_adjusted 创建完成，", nrow(wgcna_p_adjusted), "x", ncol(wgcna_p_adjusted), "矩阵\n")

# 5. 更新 target_pair - WGCNA目标配对
# ================================================================
wgcna_target_pair <- data.frame(
  Source_Type = character(0),
  Source_ID = character(0),
  Target_Type = character(0),
  Target_ID = character(0),
  Relationship = character(0),
  Strength = numeric(0),
  P_value = numeric(0),
  stringsAsFactors = FALSE
)

# 添加模块间的显著关联
for (i in 1:(nrow(wgcna_correlation_matrix)-1)) {
  for (j in (i+1):ncol(wgcna_correlation_matrix)) {
    cor_val <- wgcna_correlation_matrix[i, j]
    if (abs(cor_val) > 0.3) {  # 模块间相关性阈值
      wgcna_target_pair <- rbind(wgcna_target_pair, data.frame(
        Source_Type = "WGCNA_Module",
        Source_ID = rownames(wgcna_correlation_matrix)[i],
        Target_Type = "WGCNA_Module",
        Target_ID = colnames(wgcna_correlation_matrix)[j],
        Relationship = ifelse(cor_val > 0, "Positive_Correlation", "Negative_Correlation"),
        Strength = abs(cor_val),
        P_value = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 添加模块与表型的显著关联
for (i in 1:nrow(wgcna_p_adjusted)) {
  for (j in 1:ncol(wgcna_p_adjusted)) {
    if (!is.na(wgcna_p_adjusted[i, j]) && wgcna_p_adjusted[i, j] < 0.05) {
      cor_val <- moduleTraitCor[i, j]
      wgcna_target_pair <- rbind(wgcna_target_pair, data.frame(
        Source_Type = "WGCNA_Module",
        Source_ID = rownames(wgcna_p_adjusted)[i],
        Target_Type = "Tumor_Type",
        Target_ID = colnames(wgcna_p_adjusted)[j],
        Relationship = ifelse(cor_val > 0, "Positive_Association", "Negative_Association"),
        Strength = abs(cor_val),
        P_value = wgcna_p_adjusted[i, j],
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 添加hub基因与模块的关联
hub_genes <- wgcna_modules[wgcna_modules$Is_Hub == TRUE, ]
if (nrow(hub_genes) > 0) {
  for (i in 1:nrow(hub_genes)) {
    wgcna_target_pair <- rbind(wgcna_target_pair, data.frame(
      Source_Type = "Hub_Protein",
      Source_ID = hub_genes$Protein[i],
      Target_Type = "WGCNA_Module",
      Target_ID = paste0("ME", hub_genes$Module_Color[i]),
      Relationship = "Hub_Membership",
      Strength = abs(hub_genes$kME[i]),
      P_value = NA,
      stringsAsFactors = FALSE
    ))
  }
}

cat("✓ wgcna_target_pair 创建完成，包含", nrow(wgcna_target_pair), "个配对关系\n")



# =============================================================================
# 16. 保存所有WGCNA变量（包括新增的整合分析变量）
# =============================================================================

# 保存完整的WGCNA分析结果
save(net, MEs, moduleColors, geneTree, datExpr_filtered, moduleTraitCor, moduleTraitPvalue,
     wgcna_modules, wgcna_significant_results, wgcna_correlation_matrix, wgcna_p_adjusted, wgcna_target_pair,
     clinical_data, protein_info, softPower,
     file = "./WGCNA/wgcna_analysis_complete_with_integration_vars.RData")

cat("✓ 已保存完整的WGCNA分析结果（包含整合分析变量）到 wgcna_analysis_complete_with_integration_vars.RData\n")

# 单独保存整合分析所需的核心变量
save(wgcna_modules, wgcna_significant_results, wgcna_correlation_matrix, wgcna_p_adjusted, wgcna_target_pair,
     file = "./WGCNA/wgcna_integration_variables.RData")

cat("✓ 已保存WGCNA整合分析变量到 wgcna_integration_variables.RData\n")

# =============================================================================
# 17. 创建统一命名的整合变量（用于跨方法整合）
# =============================================================================

# 重命名变量以便与Mfuzz结果整合
significant_results <- wgcna_significant_results  # 统一命名
correlation_matrix <- wgcna_correlation_matrix    # 统一命名  
p_adjusted <- wgcna_p_adjusted                   # 统一命名
target_pair <- wgcna_target_pair                 # 统一命名

# 保存统一命名的变量
save(wgcna_modules, significant_results, correlation_matrix, p_adjusted, target_pair,
     file = "./WGCNA/wgcna_unified_integration_variables.RData")