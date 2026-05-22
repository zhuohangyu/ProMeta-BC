# 清理工作空间
rm(list = ls())
gc()

# 加载必要的包
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(grid)

# 读取数据
clinical_data <- read.csv("Clinical.csv", stringsAsFactors = FALSE)

# 计算BMI
clinical_data$BMI <- clinical_data$BMI

# 处理Clinical stage，将IA、IB等转换为I，IIA、IIB等转换为II，以此类推
clinical_data$Clinical.stage <- gsub("[ABC]$", "", clinical_data$Clinical.stage)

# 处理TNM分期，去除前缀（如果有的话）
clinical_data$T <- gsub("^T", "", clinical_data$T)
clinical_data$N <- gsub("^N", "", clinical_data$N)
clinical_data$M <- gsub("^M", "", clinical_data$M)

# 按Cohort排序数据
cohort_order <- c("Cohort_1", "Cohort_2", "Cohort_3")
clinical_data$Cohort <- factor(clinical_data$Cohort, levels = cohort_order)

# 重要：确保数据按照Cohort的因子水平排序
clinical_data <- clinical_data %>% 
  arrange(match(Cohort, cohort_order))

# 确保所有分类变量都是字符型
clinical_data$Label <- as.character(clinical_data$Label)
clinical_data$Clinical.stage <- as.character(clinical_data$Clinical.stage)
clinical_data$Molecular.Classification <- as.character(clinical_data$Molecular.Classification)
clinical_data$T <- as.character(clinical_data$T)
clinical_data$N <- as.character(clinical_data$N)
clinical_data$M <- as.character(clinical_data$M)

# 验证排序
cat("\n========== 验证Cohort排序 ==========\n")
cat("前10个样本的Cohort:", as.character(clinical_data$Cohort[1:10]), "\n")
cat("每个Cohort的样本数:\n")
print(table(clinical_data$Cohort))

# ========== 改进的配色方案 ==========

# 获取Molecular Classification的唯一值
mol_class_unique <- sort(unique(as.character(clinical_data$Molecular.Classification)))

# Cohort颜色 - 使用更专业柔和的配色（深度递进）
cohort_colors <- structure(c("#2C3E50", "#34495E", "#5D6D7E"),
                           names = cohort_order)

# Label颜色 - 更直观的医学配色：红色表示恶性，蓝色表示良性
label_colors <- c("Malignant" = "#DC2626",    # 深红色 - 恶性
                  "Benign" = "#2563EB")       # 深蓝色 - 良性


# Age颜色渐变 - 使用优雅的紫色到橙色渐变
age_col_fun <- colorRamp2(c(min(clinical_data$Age, na.rm = TRUE), 
                            quantile(clinical_data$Age, 0.25, na.rm = TRUE),
                            quantile(clinical_data$Age, 0.5, na.rm = TRUE),
                            quantile(clinical_data$Age, 0.75, na.rm = TRUE),
                            max(clinical_data$Age, na.rm = TRUE)), 
                          c("#F8F9FA", "#AED6F1", "#5DADE2", "#2E86AB", "#1B4F72"))

# BMI颜色渐变 - 使用健康指标相关的绿色到红色渐变
bmi_col_fun <- colorRamp2(c(min(clinical_data$BMI, na.rm = TRUE), 18.5, 25, 30, 
                            max(clinical_data$BMI, na.rm = TRUE)), 
                          c("#85C1E9", "#27AE60", "#F1C40F", "#E67E22", "#C0392B"))

# Molecular Classification颜色 - 使用青绿色系渐变
mol_class_colors <- c("1" = "#D5F4E6",  # 浅绿色
                      "2" = "#58D68D",  # 中绿色  
                      "3" = "#239B56")  # 深绿色

# Clinical stage颜色 - 使用渐进的紫色系
stage_colors <- c("0" = "#FADBD8", 
                  "I" = "#F1948A", 
                  "II" = "#E74C3C", 
                  "III" = "#C0392B",
                  "IV" = "#922B21")

# TNM分期颜色优化
# T分期颜色 - 使用青色系
t_colors <- c("0" = "#EBF5FB",
              "1" = "#AED6F1", 
              "2" = "#5DADE2",
              "3" = "#2E86AB",
              "4" = "#1B4F72",
              "is" = "#D6EAF8",
              "x" = "#BDC3C7")

# N分期颜色 - 使用橙色系
n_colors <- c("0" = "#FEF9E7",
              "1" = "#FCF3CF", 
              "2" = "#F7DC6F",
              "3" = "#F4D03F",
              "x" = "#BDC3C7")

# M分期颜色 - 使用红色系
m_colors <- c("0" = "#FADBD8",
              "1" = "#E74C3C",
              "x" = "#BDC3C7")

# 验证颜色分配
cat("\n========== 改进后颜色分配 ==========\n")
for(cohort in cohort_order) {
  cat(cohort, ":", cohort_colors[cohort], "\n")
}

cat("\n========== Label颜色分配 ==========\n")
for(label in names(label_colors)) {
  cat(label, ":", label_colors[label], "\n")
}

# 创建热图注释（移除Cohort行，只保留列标题）
ha <- HeatmapAnnotation(
  # 移除Cohort行，只在column_split中使用
  Label = clinical_data$Label,
  
  # 连续变量 - 调整为标准高度，添加边框
  Age = anno_simple(clinical_data$Age, 
                    col = age_col_fun, 
                    height = unit(8, "mm"),  # 从12mm减少到8mm
                    border = TRUE),
  BMI = anno_simple(clinical_data$BMI, 
                    col = bmi_col_fun, 
                    height = unit(8, "mm"),  # 从12mm减少到8mm
                    border = TRUE),
  
  # 分期信息
  `Clinical Stage` = clinical_data$Clinical.stage,
  `Molecular Class` = as.character(clinical_data$Molecular.Classification),
  
  # TNM分期
  `T Stage` = clinical_data$T,
  `N Stage` = clinical_data$N,
  `M Stage` = clinical_data$M,
  
  # 设置颜色列表
  col = list(
    Label = label_colors,
    `Clinical Stage` = stage_colors,
    `Molecular Class` = mol_class_colors,
    `T Stage` = t_colors,
    `N Stage` = n_colors,
    `M Stage` = m_colors
  ),
  
  # 注释名称样式
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
  gap = unit(1, "mm"),  # 适中的间隙
  which = "column",
  # 为所有注释添加边框
  border = TRUE,
  # 统一所有注释的高度
  simple_anno_size = unit(8, "mm")
)

# 创建一个虚拟矩阵
dummy_mat <- matrix(runif(nrow(clinical_data), 0, 0.001), nrow = 1, ncol = nrow(clinical_data))
rownames(dummy_mat) <- ""
colnames(dummy_mat) <- clinical_data$Sample

# 验证列对应关系
cat("\n========== 验证列对应关系 ==========\n")
cat("矩阵列数:", ncol(dummy_mat), "\n")
cat("数据行数:", nrow(clinical_data), "\n")

# 创建主热图
ht <- Heatmap(dummy_mat,
              name = "dummy",
              top_annotation = ha,
              show_column_names = FALSE,
              show_row_names = FALSE,
              show_heatmap_legend = FALSE,
              # 使用factor确保顺序，同时设置列标题颜色
              column_split = factor(clinical_data$Cohort, levels = cohort_order),
              column_title = cohort_order,
              column_title_gp = gpar(fontsize = 14, fontface = "bold", 
                                     col = cohort_colors[cohort_order]),
              column_gap = unit(3, "mm"),
              height = unit(1, "mm"),
              col = colorRamp2(c(0, 0.001), c("white", "white")),
              cluster_columns = FALSE,
              column_order = 1:ncol(dummy_mat),
              # 只为主热图区域添加极细边框，不影响注释颜色
              border = TRUE,
              border_gp = gpar(col = "#F8F9FA", lwd = 0.1)  # 极淡的灰色，极细线条
)

# 创建单独的图例
# Age图例
age_legend <- Legend(col_fun = age_col_fun, 
                     title = "Age (years)",
                     title_gp = gpar(fontsize = 11, fontface = "bold"),
                     labels_gp = gpar(fontsize = 10),
                     legend_height = unit(3, "cm"),
                     direction = "vertical")

# BMI图例
bmi_legend <- Legend(col_fun = bmi_col_fun, 
                     title = "BMI (kg/m²)",
                     title_gp = gpar(fontsize = 11, fontface = "bold"),
                     labels_gp = gpar(fontsize = 10),
                     legend_height = unit(3, "cm"),
                     direction = "vertical")

# Cohort图例（用于说明列标题颜色）
cohort_legend <- Legend(at = cohort_order,
                        legend_gp = gpar(fill = cohort_colors),
                        title = "Cohort",
                        title_gp = gpar(fontsize = 11, fontface = "bold"),
                        labels_gp = gpar(fontsize = 10))

# 保存高分辨率PNG
png("FigureS1A.png", width = 9600, height = 3200, res = 600)
ht_draw <- draw(ht, 
                annotation_legend_side = "right",
                heatmap_legend_side = "right")

# 添加自定义图例
draw(age_legend, x = unit(0.95, "npc"), y = unit(0.8, "npc"), just = c("left", "top"))
draw(bmi_legend, x = unit(0.95, "npc"), y = unit(0.5, "npc"), just = c("left", "top"))
draw(cohort_legend, x = unit(0.95, "npc"), y = unit(0.2, "npc"), just = c("left", "top"))

dev.off()