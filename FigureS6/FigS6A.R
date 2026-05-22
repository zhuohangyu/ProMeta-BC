## =========================================================
## VennDiagram 高级美化版
## =========================================================
library(tidyverse)
library(VennDiagram)
library(grid)
library(RColorBrewer)

# 读取数据并筛选（保持你原来的代码）
fc <- read.csv("DEGs.csv", check.names = FALSE)
gene_col <- colnames(fc)[1]

fc_filt <- fc %>%
  filter(`regulation_malignant-vs-benign` != "unchanged")

genes_HR   <- fc_filt %>%
  filter(!is.na(`regulation_HR+`), `regulation_HR+` != "unchanged") %>%
  pull(!!sym(gene_col)) %>%
  unique()

genes_HER  <- fc_filt %>%
  filter(!is.na(`regulation_HER+`), `regulation_HER+` != "unchanged") %>%
  pull(!!sym(gene_col)) %>%
  unique()

genes_TNBC <- fc_filt %>%
  filter(!is.na(`regulation_TNBC`), `regulation_TNBC` != "unchanged") %>%
  pull(!!sym(gene_col)) %>%
  unique()

venn_list <- list(
  "HR+"   = genes_HR,
  "HER2+" = genes_HER,
  "TNBC"  = genes_TNBC
)

## 美化版 Venn 图
grid.newpage()
venn_plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  
  # 填充颜色 - 使用柔和的学术配色
  fill = c("#E57373", "#81C784", "#64B5F6"),  # 红、绿、蓝
  
  # 透明度
  alpha = 0.5,
  
  # 圆圈边框
  lwd = 2.5,                    # 边框粗细
  col = c("#C62828", "#388E3C", "#1976D2"),  # 边框颜色（深色版）
  lty = "solid",                # 实线
  
  # 区域内数字样式
  cex = 2,                      # 数字大小
  fontface = "bold",            # 粗体
  fontfamily = "sans",          # 字体
  
  # 类别标签样式
  cat.cex = 2.2,                # 标签字体大小
  cat.fontface = "bold",        # 标签粗体
  cat.fontfamily = "sans",
  cat.col = c("#C62828", "#388E3C", "#1976D2"),  # 标签颜色
  
  # 标签位置微调
  cat.pos = c(-25, 25, 180),    # 角度位置
  cat.dist = c(0.055, 0.055, 0.055),  # 距离
  
  # 旋转
  rotation.degree = 0,
  
  # 其他美化
  margin = 0.1                  # 边距
)

grid.draw(venn_plot)

# 保存高分辨率图片
png("FigureS6A.png", 
    width = 3000, height = 3000, res = 400)
grid.newpage()
grid.draw(venn_plot)
dev.off()