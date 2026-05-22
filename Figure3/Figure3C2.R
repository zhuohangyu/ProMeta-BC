# 清理工作空间
rm(list = ls())
gc()

# 简化版蛋白质模块富集分析
# =====================================

# 安装和加载必要的包
if (!require("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!require("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!require("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# 读取数据
data <- read.csv("./WGCNA/protein_module_colors_with_numbers.csv")
print(paste("总蛋白质数：", nrow(data)))
print("模块分布：")
print(table(data$ModuleLabel))

# 按模块分组并进行富集分析
modules <- unique(data$ModuleLabel)
modules <- modules[modules != "ME00"]  # 排除灰色模块
modules <- c("ME10")
for (module in modules) {
  cat("\n=== 分析模块", module, "===\n")
  
  # 获取该模块的蛋白质
  module_proteins <- data$Protein[data$ModuleLabel == module]
  cat("蛋白质数量：", length(module_proteins), "\n")
  
  # ID转换
  tryCatch({
    gene_ids <- bitr(module_proteins, 
                     fromType = "UNIPROT", 
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db,
                     drop = TRUE)
    
    if (nrow(gene_ids) > 3) {  # 至少需要3个基因
      
      # GO生物过程分析
      go_bp <- enrichGO(gene = gene_ids$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        pvalueCutoff = 0.2,
                        qvalueCutoff = 0.5,
                        readable = TRUE)
      
      # KEGG通路分析
      kegg <- enrichKEGG(gene = gene_ids$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 0.2,
                         qvalueCutoff = 0.5)
      
      # 打印结果摘要
      cat("GO-BP 显著条目数：", nrow(go_bp@result), "\n")
      cat("KEGG 显著通路数：", nrow(kegg@result), "\n")
      
      # 绘制图表
      if (nrow(go_bp@result) > 0) {
        p1 <- dotplot(go_bp, showCategory = 10, 
                      title = paste(module, "GO Biological Process"))
        print(p1)
        
        # 保存图片
        ggsave(paste0(module, "_GO_BP.png"), p1, width = 10, height = 8)
      }
      
      if (nrow(kegg@result) > 0) {
        p2 <- dotplot(kegg, showCategory = 10, 
                      title = paste(module, "KEGG Pathway"))
        print(p2)
        
        # 保存图片
        ggsave(paste0(module, "_KEGG.png"), p2, width = 10, height = 8)
      }
      
      # 保存结果表格
      if (nrow(go_bp@result) > 0) {
        write.csv(go_bp@result, paste0(module, "_GO_BP_results.csv"), row.names = FALSE)
      }
      if (nrow(kegg@result) > 0) {
        write.csv(kegg@result, paste0(module, "_KEGG_results.csv"), row.names = FALSE)
      }
      
    } else {
      cat("基因数量太少，跳过富集分析\n")
    }
    
  }, error = function(e) {
    cat("ID转换失败：", e$message, "\n")
  })
}

cat("\n分析完成！\n")