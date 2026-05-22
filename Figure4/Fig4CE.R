# 清理工作空间
rm(list = ls())
gc()


## ===== 1) 加载包 =====
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(tibble)
library(grid)

## ===== 2) 读入数据 =====
deg <- read_csv("DEGs.csv")
clin <- read_csv("Clinical.csv")

## ===== 3) 筛选代谢物（四个比较都不是 unchanged）=====
deg_filt <- deg %>%
  filter(`regulation_malignant-vs-benign` != "unchanged",
         `regulation_HR+`  != "unchanged",
         `regulation_HER+` != "unchanged",
         `regulation_TNBC` != "unchanged")

## ===== 4) 构建表达矩阵（行=代谢物，列=样本）=====
## 找到 Category1 的列号
cat_col <- which(colnames(deg_filt) == "Category1")




## ===== 6) 行信息（log2FC, FDR, VIP）=====
## ===== 代谢物排序 =====
deg_filt <- deg_filt %>%
  arrange(Category1, `log2fc_malignant-vs-benign`)

## 表达矩阵：从 Category1 后一列开始都是样本
expr_mat <- deg_filt %>%
  select(`Metabolite name`, (cat_col+1):ncol(.)) %>%
  column_to_rownames("Metabolite name") %>%
  as.matrix()
## ===== Z-score 标准化 (行方向) =====
expr_mat <- t(scale(t(expr_mat)))

## ===== 5) 列分组（根据 Clinical.csv）=====
group <- clin$`Molecular Classification`
names(group) <- clin$Sample
group <- group[colnames(expr_mat)]

## ⭐ 数值 → 文字 + 处理 NA
group <- dplyr::recode(
  group,
  `1` = "HR+",
  `2` = "HER2+",
  `3` = "TNBC",
  .default = "Benign",        # 包括 NA 自动变 Benign
  .missing = "Benign"
)

## ⭐ 明确顺序
group <- factor(group, levels = c("Benign","HR+","HER2+","TNBC"))

length(group)
ncol(expr_mat)

## 行信息（和之前一样，但要用排好序的 deg_filt）
row_df <- deg_filt %>%
  transmute(
    feature = `Metabolite name`,
    Class   = Category1,
    log2FC_0 = `log2fc_malignant-vs-benign`, 
    FDR_0    = `fdr_malignant-vs-benign`, 
    VIP_0    = `VIP_malignant-vs-benign`,
    log2FC_1 = `log2fc_HR+`,   FDR_1 = `fdr_HR+`,   VIP_1 = `VIP_HR+`,
    log2FC_2 = `log2fc_HER+`,  FDR_2 = `fdr_HER+`,  VIP_2 = `VIP_HER+`,
    log2FC_3 = `log2fc_TNBC`,  FDR_3 = `fdr_TNBC`,  VIP_3 = `VIP_TNBC`
  ) %>% column_to_rownames("feature")



## ===== 7) 颜色映射 =====
# 表达矩阵动态阈值
col_heat <- colorRamp2(c(-2,  0, 2), c("#47bdeb", "white", "#f9632c"))

# 类别颜色
class_cols <- c(
  "Amino Acid"="#4db6ac",
  "Carbohydrate"="#fff59d",
  "Lipid"="#b39ddb",
  "Other"="#ef9a9a",
  "Xenobiotics"="#90caf9"
)

# FDR 颜色
neglog10 <- function(x) -log10(x + .Machine$double.eps)
fdr_all <- range(neglog10(c(row_df$FDR_0,row_df$FDR_1,row_df$FDR_2,row_df$FDR_3)), na.rm=TRUE)
col_fdr <- colorRamp2(seq(fdr_all[1], fdr_all[2], length.out=4),
                      c("#f1b6da","#de77ae","#c51b7d","#7a0177"))

# VIP 缩放
scale_vip_fun <- function(vip_vec, min_pt=1.5, max_pt=2.8){
  vr <- range(vip_vec, na.rm=TRUE)
  if(diff(vr)==0) return(rep((min_pt+max_pt)/2, length(vip_vec)))
  (vip_vec - vr[1])/(vr[2]-vr[1]) * (max_pt-min_pt) + min_pt
}

## ===== 8) 左侧类别色块 =====
## 行分组因子
row_split_factor <- factor(row_df$Class,
                           levels=c("Amino Acid","Carbohydrate","Lipid",
                                    "Xenobiotics", "Other"))
row_lvls <- levels(row_split_factor)

left_block <- rowAnnotation(
  ClassBlock = anno_block(
    gp = gpar(fill = unname(class_cols[row_lvls]), col = NA),
    labels = row_lvls,
    labels_gp = gpar(col = "black", fontface = 2, cex = 0.9),
    labels_rot = 0   # 🔑 横向文字
  ),
  width = unit(25, "mm")   # 🔑 增加宽度
)

## ===== 9) 右侧四个 lollipop =====
pt0 <- scale_vip_fun(row_df$VIP_0); col0 <- col_fdr(neglog10(row_df$FDR_0))
pt1 <- scale_vip_fun(row_df$VIP_1); col1 <- col_fdr(neglog10(row_df$FDR_1))
pt2 <- scale_vip_fun(row_df$VIP_2); col2 <- col_fdr(neglog10(row_df$FDR_2))
pt3 <- scale_vip_fun(row_df$VIP_3); col3 <- col_fdr(neglog10(row_df$FDR_3))
## 1) 把 VIP 映射到点的直径（mm）
vip_size0 <- unit(scale_vip_fun(row_df$VIP_0, min_pt = 1.5, max_pt = 4.0), "mm")
vip_size1 <- unit(scale_vip_fun(row_df$VIP_1, min_pt = 1.5, max_pt = 4.0), "mm")
vip_size2 <- unit(scale_vip_fun(row_df$VIP_2, min_pt = 1.5, max_pt = 4.0), "mm")
vip_size3 <- unit(scale_vip_fun(row_df$VIP_3, min_pt = 1.5, max_pt = 4.0), "mm")

## 2) 右侧注释改写：用 size 控制点大小，gp 只控制颜色
anno_right <- rowAnnotation(
  "ALL" = anno_points(
    row_df$log2FC_0,
    pch = 16,
    size = vip_size0,
    gp   = gpar(col = col0),
    axis = TRUE,
    axis_param = list(
      at = c(-1, 0, 1),
      labels = c("-1","0","1"),
      gp = gpar(fontsize = 14)
    ),
    width = unit(20, "mm")
  ),
  
  "HR+" = anno_points(
    row_df$log2FC_1,
    pch = 16,
    size = vip_size1,
    gp   = gpar(col = col1),
    axis = TRUE,
    axis_param = list(
      at = c(-1, 0, 1),
      labels = c("-1","0","1"),
      gp = gpar(fontsize = 14)
    ),
    width = unit(20, "mm")
  ),
  
  "HER2+" = anno_points(
    row_df$log2FC_2,
    pch = 16,
    size = vip_size2,
    gp   = gpar(col = col2),
    axis = TRUE,
    axis_param = list(
      at = c(-1, 0, 1),
      labels = c("-1","0","1"),
      gp = gpar(fontsize = 14)
    ),
    width = unit(20, "mm")
  ),
  
  "TNBC" = anno_points(
    row_df$log2FC_3,
    pch = 16,
    size = vip_size3,
    gp   = gpar(col = col3),
    axis = TRUE,
    axis_param = list(
      at = c(-1, 0, 1),
      labels = c("-1","0","1"),
      gp = gpar(fontsize = 14)
    ),
    width = unit(20, "mm")
  ),
  
  annotation_name_side = "top",
  annotation_name_rot  = 0,   # ✅标题横着显示
  annotation_name_offset = unit(3, "mm"),
  annotation_name_gp = gpar(fontsize = 14),
  gap = unit(2, "mm")   # ✅气泡图之间的缝隙（推荐 3–6mm）
)


## ===== 10) 图例 =====
lgd_fdr <- Legend(
  title = expression(-log[10](FDR)),
  col_fun = col_fdr,
  at = floor(seq(fdr_all[1], fdr_all[2], length.out=5)), 
  direction = "vertical",
  title_gp  = gpar(fontsize = 18, fontface = "bold"),
  labels_gp = gpar(fontsize = 16)
)
lgd_vip <- Legend(
  title = "VIP (size)",
  type = "points", pch=16,
  size = unit(c(1.5,2.2,2.8),"mm"),
  labels = c("1","2","3"), 
  direction = "vertical",
  legend_gp = gpar(col="grey20"),
  title_gp  = gpar(fontsize = 18, fontface = "bold"),
  labels_gp = gpar(fontsize = 16)
)

## ===== 11) 主热图 =====
ht <- Heatmap(
  expr_mat,
  name="Z-score",
  col=col_heat,
  cluster_rows=FALSE, cluster_columns=FALSE,
  row_split=row_split_factor,
  column_split=group,
  show_row_names=TRUE,
  row_names_side="left",
  row_names_gp = gpar(fontsize = 16),
  row_names_max_width = unit(12, "cm"),  # 🔑 预留更多空间给名字
  show_column_names=FALSE,
  left_annotation=left_block,
  right_annotation=anno_right,
  row_title = NULL,               # 🔑 删除 row_split 的文字标题
  heatmap_legend_param=list(
    at = seq(-2, 2, by = 1),          # 🔁 整数刻度
    labels = seq(-2, 2, by = 1),
    direction = "vertical",
    title_gp  = gpar(fontsize = 18, fontface = "bold"),
    labels_gp = gpar(fontsize = 16)
  )
)
length(group)
dim(expr_mat)

png("Figure4E.png", width = 5500, height = 3000, res = 300)  # 300DPI，适合论文

## ===== 12) 绘制热图 =====
draw(ht,
     heatmap_legend_list=list(lgd_fdr, lgd_vip),
     merge_legends=FALSE,
     annotation_legend_side="right", 
     heatmap_legend_side="right")

## ===== 13) 给四个 lollipop 每个 split 都加虚线 =====
anns <- c("ALL","HR+","HER2+","TNBC")
n_slices <- length(levels(row_split_factor))

for (ann in anns) {
  for (s in 1:n_slices) {
    decorate_annotation(ann, slice = s, {
      ## 竖直虚线
      grid.lines(x = unit(c(0.5,0.5),"npc"),
                 y = unit(c(0,1),"npc"),
                 gp = gpar(lty=3,col="grey40"))
      
      ## log2FC 值
      xvals <- switch(ann,
                      "ALL" = row_df$log2FC_0,
                      "HR+"   = row_df$log2FC_1,
                      "HER2+" = row_df$log2FC_2,
                      "TNBC"  = row_df$log2FC_3)
      
      xr <- range(xvals, na.rm=TRUE)
      xscaled <- (xvals - xr[1])/(xr[2]-xr[1])
      yscaled <- seq(0,1,length.out=length(xvals))
    })
  }
}

dev.off()   # 关闭设备，图片写入文件


## ===== 14) Malignant vs Benign 火山图 (美化版) =====
library(ggplot2)
library(ggrepel)  # 用于标注代谢物名称

# 构建火山图数据
vol_data <- deg %>%
  mutate(
    negLog10FDR = -log10(`fdr_malignant-vs-benign` + .Machine$double.eps),
    Regulation  = factor(`regulation_malignant-vs-benign`,
                         levels = c("down", "unchanged", "up"))
  )

# 定义阈值
fc_cutoff  <- 1.0   # log2FC阈值
fdr_cutoff <- 0.05  # FDR阈值

# 统计上调、下调数量
n_up   <- sum(vol_data$Regulation == "up", na.rm = TRUE)
n_down <- sum(vol_data$Regulation == "down", na.rm = TRUE)

# 配色方案（更专业的科学期刊配色）
col_reg <- c(
  "down"      = "#47bdeb",  # 深蓝色
  "unchanged" = "#CCCCCC",  # 灰色
  "up"        = "#f9632c"   # 深红色
)

# 选择 top 显著代谢物进行标注（可根据需要调整）
top_n_each <- 5

top_up <- vol_data %>%
  filter(Regulation == "up") %>%
  arrange(`fdr_malignant-vs-benign`) %>%     # FDR 越小越显著
  slice_head(n = top_n_each)

top_down <- vol_data %>%
  filter(Regulation == "down") %>%
  arrange(`fdr_malignant-vs-benign`) %>%
  slice_head(n = top_n_each)

vol_data_label <- bind_rows(top_up, top_down)

p_volcano <- ggplot(
  vol_data,
  aes(x = `log2fc_malignant-vs-benign`,
      y = negLog10FDR,
      color = Regulation,
      size = `VIP_malignant-vs-benign`)
) +
  geom_hline(yintercept = -log10(fdr_cutoff), 
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_point(alpha = 0.7, shape = 21, stroke = 0.3, 
             aes(fill = Regulation), color = "white") +
  geom_text_repel(
    data = vol_data_label,
    aes(label = `Metabolite name` ),
    size = 5,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.3,
    segment.color = "grey50",
    max.overlaps = 20,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  
  ## 配色
  scale_fill_manual(values = col_reg, name = "") +
  scale_color_manual(values = col_reg, guide = "none") +
  scale_size_continuous(
    name = "VIP",
    range = c(0.5, 4),
    breaks = c(1, 2, 3, 4)
  ) +
  
  ## ⭐ 这里控制两个图例上下各一行
  guides(
    fill = guide_legend(
      order = 1,           # Regulation 在上面
      # title = "Regulation",
      override.aes = list(
        shape = 21,
        size  = 4,
        color = "white"
      )
    ),
    size = guide_legend(
      order = 2,           # VIP 在下面
      title = "VIP",
      override.aes = list(
        shape = 21,
        fill  = "grey60",
        color = "white"
      )
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.1))
  ) +
  labs(
    x = "log2 Fold Change",
    y = "log10 Pvalue",
    subtitle = sprintf("Metabolite: Up %d | Down %d", n_up, n_down)
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    plot.subtitle = element_text(size = 22, hjust = 0.5, color = "grey30"),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(color = "black", size = 20),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    
    legend.position = "bottom",
    legend.box = "vertical",              # ⭐ 关键：让两个图例上下排列
    
    # legend.spacing.y = unit(0.15, "cm"),   # ⭐上下图例间距变紧
    # legend.key.height = unit(0.4, "cm"),  # ⭐每个条目更紧凑
    legend.title = element_text(face = "bold", size = 19),
    legend.text  = element_text(size = 19),
    legend.key.size = unit(0.8, "lines"),
    legend.background = element_rect(fill = "white", color = "grey80"),
    
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
    panel.grid.minor = element_blank(), 
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  filename = "Figure4C.png",
  plot     = p_volcano,
  width    = 5.2,
  height   = 6,
  dpi      = 300,
  bg       = "white"
)
