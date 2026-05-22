rm(list = ls())

# 加载必要的库
library(dplyr)
library(ggplot2)
library(ggthemes)
library(readr)
library(openxlsx)
library(patchwork)
library(tidyr)
library(stringr)  # 用于字符串处理

# 读取数据文件
clinical_data <- read_csv("Clinical.csv")
module_data <- read_csv("./WGCNA/ModuleEigengenes_All_Samples.csv")

# 检查数据结构
print("临床数据结构:")
print(head(clinical_data))
print("模块数据结构:")
print(head(module_data))

# 处理模块数据 - 重命名第一列为SampleID并转置
colnames(module_data)[1] <- "SampleID"

# 保存样本ID
sample_ids <- module_data$SampleID

# 使用基础R方法设置行名并转置
rownames(module_data) <- module_data$SampleID
module_data_for_transpose <- module_data[, -1]  # 移除SampleID列
module_data_t <- t(module_data_for_transpose)
module_data_t <- as.data.frame(module_data_t)

# 确保转置后的列名是样本ID
colnames(module_data_t) <- sample_ids

print("转置后的模块数据维度:", dim(module_data_t))
# print("前几个样本ID:", head(colnames(module_data_t)))

# 获取模块名称
module_names <- rownames(module_data_t)
print(paste("模块数量:", length(module_names)))
print(paste("模块名称:", paste(module_names, collapse = ", ")))

# 重命名临床数据的Sample列为SampleID以便合并
clinical_data <- clinical_data %>%
  rename(SampleID = Sample)

# 检查原始数据的列名和BMI计算条件
print("原始临床数据列名:")
print(colnames(clinical_data))
print(paste("是否包含Height:", "Height" %in% colnames(clinical_data)))
print(paste("是否包含Weight:", "Weight" %in% colnames(clinical_data)))
print(paste("是否已包含BMI:", "BMI" %in% colnames(clinical_data)))

# ================= 全局字体放大 =================
base_big <- theme(
  axis.text.x   = element_text(size = 15, angle = 45, hjust = 1),
  axis.text.y   = element_text(size = 15),
  axis.title.x  = element_text(size = 16, face = "bold"),
  axis.title.y  = element_text(size = 16, face = "bold"),
  plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
  legend.title  = element_text(size = 14, face = "bold"),
  legend.text   = element_text(size = 13)
)

# 强制计算BMI
if ("Height" %in% colnames(clinical_data) && "Weight" %in% colnames(clinical_data)) {
  print("开始计算BMI...")
  
  # 计算前检查数据
  height_count <- sum(!is.na(clinical_data$Height))
  weight_count <- sum(!is.na(clinical_data$Weight))
  print(paste("Height非缺失值数量:", height_count))
  print(paste("Weight非缺失值数量:", weight_count))
  
  # 计算BMI
  clinical_data <- clinical_data %>%
    mutate(BMI = Weight / (Height/100)^2)
  
  # 检查BMI计算结果
  bmi_count <- sum(!is.na(clinical_data$BMI))
  print(paste("BMI计算完成，非缺失值数量:", bmi_count))
  print(paste("BMI前5个值:", paste(head(clinical_data$BMI, 5), collapse = ", ")))
  
  if ("BMI" %in% colnames(clinical_data)) {
    print("✓ BMI列已成功添加到clinical_data")
  } else {
    print("✗ BMI列添加失败")
  }
} else {
  print("✗ 缺少Height或Weight列，无法计算BMI")
}

print("更新后的临床数据列名:")
print(colnames(clinical_data))

# =============================================================================
# 数据预处理：标准化临床分期和TNM分期
# =============================================================================

print("\n=== 开始数据预处理 ===")

# 1. 标准化临床分期 (Clinical stage) 为 0-I-II-III-IV 格式
if ("Clinical stage" %in% colnames(clinical_data)) {
  print("处理Clinical stage...")
  
  # 查看原始数据分布
  print("原始Clinical stage分布:")
  print(table(clinical_data$`Clinical stage`, useNA = "ifany"))
  
  clinical_data <- clinical_data %>%
    mutate(`Clinical stage` = case_when(
      # 标准化各种0期表示
      str_detect(toupper(`Clinical stage`), "^0$|^STAGE\\s*0$|^TIS") ~ "0",
      # 标准化各种I期表示  
      str_detect(toupper(`Clinical stage`), "^I$|^STAGE\\s*I$|^1$|^IA$|^IB$|^IC$") ~ "I",
      # 标准化各种II期表示
      str_detect(toupper(`Clinical stage`), "^II$|^STAGE\\s*II$|^2$|^IIA$|^IIB$|^IIC$") ~ "II",
      # 标准化各种III期表示
      str_detect(toupper(`Clinical stage`), "^III$|^STAGE\\s*III$|^3$|^IIIA$|^IIIB$|^IIIC$") ~ "III",
      # 标准化各种IV期表示
      str_detect(toupper(`Clinical stage`), "^IV$|^STAGE\\s*IV$|^4$|^IVA$|^IVB$|^IVC$") ~ "IV",
      # 保持NA值
      is.na(`Clinical stage`) ~ NA_character_,
      # 其他值保持原样（可能需要手动检查）
      TRUE ~ as.character(`Clinical stage`)
    ))
  
  print("标准化后Clinical stage分布:")
  print(table(clinical_data$`Clinical stage`, useNA = "ifany"))
}

# 2. 处理TNM分期：排除包含X的值
tnm_vars <- c("T", "N", "M")

for (var in tnm_vars) {
  if (var %in% colnames(clinical_data)) {
    print(paste("处理", var, "分期..."))
    
    # 查看原始数据分布
    print(paste("原始", var, "分期分布:"))
    print(table(clinical_data[[var]], useNA = "ifany"))
    
    # 将包含X的值设为NA（排除）
    clinical_data <- clinical_data %>%
      mutate(!!sym(var) := case_when(
        # 排除包含X的值
        str_detect(toupper(!!sym(var)), "X") ~ NA_character_,
        # 保持其他值
        TRUE ~ as.character(!!sym(var))
      ))
    
    print(paste("处理后", var, "分期分布:"))
    print(table(clinical_data[[var]], useNA = "ifany"))
  }
}

# 3. 检查处理后的数据质量
print("\n=== 数据预处理完成 ===")
print("各变量的有效样本数:")
for (var in c("Clinical stage", "T", "N", "M")) {
  if (var %in% colnames(clinical_data)) {
    valid_count <- sum(!is.na(clinical_data[[var]]) & clinical_data[[var]] != "")
    print(paste(var, ":", valid_count, "个有效样本"))
  }
}

# 检查两个数据集的样本重叠情况
common_samples <- intersect(clinical_data$SampleID, colnames(module_data_t))
print(paste("共同样本数量:", length(common_samples)))

# 筛选有效数据
clinical_filtered <- clinical_data %>%
  filter(SampleID %in% common_samples)

module_filtered <- module_data_t[, colnames(module_data_t) %in% common_samples]

print(paste("筛选后临床数据样本数:", nrow(clinical_filtered)))
print(paste("筛选后模块数据样本数:", ncol(module_filtered)))

# =============================================================================
# 1. 连续变量关联分析 (年龄、BMI)
# =============================================================================

continuous_vars <- c("Age", "BMI")
print(paste("分析的连续变量:", paste(continuous_vars, collapse = ", ")))

# 初始化一个空的数据框，预定义列名
plot_info_continuous <- data.frame(
  Module = character(0),
  Factor = character(0),
  Type = character(0),
  Correlation = numeric(0),
  P_value = numeric(0),
  Test = character(0),
  stringsAsFactors = FALSE
)

for (var in continuous_vars) {
  print(paste("正在处理变量:", var))
  
  if (var %in% colnames(clinical_filtered)) {
    print(paste("  变量", var, "存在于数据中"))
    
    # 筛选非缺失值的样本
    valid_samples <- clinical_filtered %>%
      filter(!is.na(!!sym(var))) %>%
      pull(SampleID)
    
    print(paste("  有效样本数:", length(valid_samples)))
    
    if (length(valid_samples) > 10) {  # 确保有足够的样本
      print(paste("  样本数足够，开始分析"))
      
      plot_data <- clinical_filtered %>%
        filter(SampleID %in% valid_samples)
      
      module_subset <- module_filtered[, colnames(module_filtered) %in% valid_samples]
      print(paste("  模块数据维度:", nrow(module_subset), "x", ncol(module_subset)))
      
      for (i in 1:length(module_names)) {
        module_scores <- as.numeric(module_subset[module_names[i], plot_data$SampleID])
        clinical_values <- plot_data[[var]]
        
        # 计算Spearman相关性
        if (length(module_scores) == length(clinical_values) && 
            length(module_scores) > 2) {
          cor_result <- cor.test(module_scores, clinical_values, method = "spearman")
          
          # 创建新行数据框
          new_row <- data.frame(
            Module = module_names[i],
            Factor = var,
            Type = "Continuous",
            Correlation = as.numeric(cor_result$estimate),
            P_value = as.numeric(cor_result$p.value),
            Test = "Spearman correlation",
            stringsAsFactors = FALSE
          )
          plot_info_continuous <- rbind(plot_info_continuous, new_row)
          
          if (i == 1) print(paste("  成功添加", var, "的分析结果"))
        }
      }
    } else {
      print(paste("  样本数不足 (", length(valid_samples), "≤ 10), 跳过", var))
    }
  } else {
    print(paste("  警告: 变量", var, "不存在于数据中"))
    print(paste("  可用变量:", paste(colnames(clinical_filtered), collapse = ", ")))
  }
}

# 处理连续变量数据类型
if (nrow(plot_info_continuous) > 0) {
  plot_info_continuous <- plot_info_continuous %>%
    mutate(
      LogP = -log10(P_value),
      Padj = p.adjust(P_value, method = "fdr")
    ) %>%
    mutate(
      Tag = case_when(
        Correlation < 0 & P_value < 0.05 ~ "Negative",
        Correlation > 0 & P_value < 0.05 ~ "Positive",
        .default = "notSig"
      )
    )
}

# =============================================================================
# 2. 分类变量关联分析
# =============================================================================

categorical_vars <- c("Label", "Clinical stage", "Molecular Classification", "T", "N", "M")

# 初始化一个空的数据框，预定义列名
plot_info_categorical <- data.frame(
  Module = character(0),
  Factor = character(0),
  Type = character(0),
  Correlation = numeric(0),
  P_value = numeric(0),
  Test = character(0),
  stringsAsFactors = FALSE
)

for (var in categorical_vars) {
  print(paste("正在处理分类变量:", var))
  
  if (var %in% colnames(clinical_filtered)) {
    print(paste("  变量", var, "存在于数据中"))
    
    # 筛选非缺失值的样本（排除NA、空字符串、"X"相关值）
    valid_samples <- clinical_filtered %>%
      filter(!is.na(!!sym(var)) & !!sym(var) != "" & !str_detect(toupper(!!sym(var)), "X")) %>%
      pull(SampleID)
    
    print(paste("  有效样本数:", length(valid_samples)))
    
    if (length(valid_samples) > 10) {
      print(paste("  样本数足够，开始分析"))
      
      plot_data <- clinical_filtered %>%
        filter(SampleID %in% valid_samples)
      
      module_subset <- module_filtered[, colnames(module_filtered) %in% valid_samples]
      print(paste("  模块数据维度:", nrow(module_subset), "x", ncol(module_subset)))
      
      # 检查分类数量
      unique_categories <- unique(plot_data[[var]])
      n_categories <- length(unique_categories)
      print(paste("  分类数量:", n_categories, "类别:", paste(unique_categories, collapse = ", ")))
      
      if (n_categories > 1) {  # 至少要有2个分类才能分析
        for (i in 1:length(module_names)) {
          module_scores <- as.numeric(module_subset[module_names[i], plot_data$SampleID])
          categories <- plot_data[[var]]
          
          if (length(module_scores) == length(categories) && 
              length(module_scores) > 2) {
            
            # 根据分类数量选择统计检验
            if (n_categories == 2) {
              # 两组比较：Wilcoxon检验
              test_result <- wilcox.test(module_scores ~ factor(categories))
              test_name <- "Wilcoxon rank-sum test"
              test_stat <- NA  # Wilcoxon检验没有相关系数
            } else {
              # 多组比较：Kruskal-Wallis检验
              test_result <- kruskal.test(module_scores ~ factor(categories))
              test_name <- "Kruskal-Wallis test"
              test_stat <- NA  # K-W检验没有相关系数
            }
            
            # 创建新行数据框
            new_row <- data.frame(
              Module = module_names[i],
              Factor = var,
              Type = "Categorical",
              Correlation = test_stat,
              P_value = as.numeric(test_result$p.value),
              Test = test_name,
              stringsAsFactors = FALSE
            )
            plot_info_categorical <- rbind(plot_info_categorical, new_row)
            
            if (i == 1) print(paste("  成功添加", var, "的分析结果"))
          }
        }
      } else {
        print(paste("  分类数不足 (", n_categories, "≤ 1), 跳过", var))
      }
    } else {
      print(paste("  样本数不足 (", length(valid_samples), "≤ 10), 跳过", var))
    }
  } else {
    print(paste("  警告: 变量", var, "不存在于数据中"))
    print(paste("  可用变量:", paste(colnames(clinical_filtered), collapse = ", ")))
  }
}

# 处理分类变量数据类型
if (nrow(plot_info_categorical) > 0) {
  plot_info_categorical <- plot_info_categorical %>%
    mutate(
      LogP = -log10(P_value),
      Padj = p.adjust(P_value, method = "fdr")
    ) %>%
    mutate(
      Tag = ifelse(P_value < 0.05, "Significant", "notSig")
    )
}
# ===== 给连续 / 分类变量加一个截断后的 LogP_cc（0–3 用于共同图例） =====
if (nrow(plot_info_continuous) > 0) {
  plot_info_continuous <- plot_info_continuous %>%
    mutate(LogP_cc = pmin(LogP, 3))   # 0–3 之间
}

if (nrow(plot_info_categorical) > 0) {
  plot_info_categorical <- plot_info_categorical %>%
    mutate(LogP_cc = pmin(LogP, 3))   # 0–3 之间
}

# ========= 新增：只保留 Label，并改名为 Malignant =========
# 单独拎出 Label，并改名为 Malignant
plot_label <- plot_info_categorical %>%
  filter(Factor == "Label") %>%
  mutate(Factor = "Malignant")

# 去除 Label 之后剩余分类变量
# 去除 Label 之后剩余分类变量，并指定显示顺序
desired_order <- c("Clinical stage", "Molecular Classification", "T", "N", "M")

plot_categorical_other <- plot_info_categorical %>%
  filter(Factor != "Label") %>%
  mutate(
    Factor = factor(Factor, levels = desired_order)
  )



# =============================================================================
# 3. 合并结果并可视化
# =============================================================================

# 合并所有结果
all_results <- list()
if (nrow(plot_info_continuous) > 0) {
  plot_info_continuous$Dataset <- "Continuous"
  all_results$continuous <- plot_info_continuous
}
if (nrow(plot_info_categorical) > 0) {
  plot_info_categorical$Dataset <- "Categorical"
  all_results$categorical <- plot_info_categorical
}

# 合并所有数据
# ==== 合并所有数据后 ====
if (length(all_results) > 0) {
  combined_data <- do.call(rbind, all_results)
  
  combined_data$Module <- factor(combined_data$Module,
                                 levels = rev(module_names))
  
  # ===== 两套 size 图例 =====
  # 1) Continuous + Categorical 共同的 0–3 图例
  size_cc <- scale_size_continuous(
    name   = "-log10(P)",
    range  = c(1.5, 6),
    limits = c(0, 3),
    breaks = c(1, 2, 3),
    labels = c(1, 2, 3),
    guide  = guide_legend(override.aes = list(shape = 16))
  )
  
  # 2) Malignant 单独的高值图例（自动按实际范围给 3 个刻度，一般就是 25–30）
  label_min <- floor(min(plot_label$LogP-1, na.rm = TRUE))
  label_max <- ceiling(max(plot_label$LogP+1, na.rm = TRUE))
  label_breaks <- c(28, 29, 30, 31)
  
  size_label <- scale_size_continuous(
    name   = "-log10(P)",
    range  = c(1.5, 6),
    limits = c(label_min, label_max),
    breaks = label_breaks,
    labels = label_breaks,
    guide  = guide_legend(override.aes = list(shape = 16))
  )
  # ====================== 画图 ==========================
  if (nrow(plot_info_continuous) > 0) {
    p_continuous <- ggplot(plot_info_continuous, aes(Factor, Module)) +
      geom_point(aes(size = LogP_cc, colour = Tag), shape = 16) +
      size_cc +                                # 用共同 0–3 图例
      xlab("") + ylab("") +
      theme_base() +
      labs(title = "Continuous Variables") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.background = element_blank(),
        panel.border = element_blank()
      ) +
      scale_colour_manual(
        name   = "Correlation",
        values = c(Negative = "#009632",
                   Positive = "#8f00b7",
                   notSig   = "#CCCCCC"),
        guide = guide_legend(
          override.aes = list(size = 5)  # 图例里的点变大
        )
      )
  }
  
  if (nrow(plot_categorical_other) > 0) {
    p_categorical_other <- ggplot(plot_categorical_other, aes(Factor, Module)) +
      geom_point(aes(size = LogP_cc, colour = Tag), shape = 16) +
      size_cc +  
      scale_x_discrete(
        labels = c(
          "Clinical stage" = "Clinical stage",
          "Molecular Classification" = "Molecular\nClassification",
          "T" = "T",
          "N" = "N",
          "M" = "M"
        )
      ) +
      xlab("") + ylab("") +
      theme_base() +
      labs(title = "Categorical Variables") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.background = element_blank(),
        panel.border = element_blank()
      ) +
      scale_colour_manual(
        name   = "Significance",
        values = c(Significant = "#ed7a00",
                   notSig     = "#CCCCCC"),
        guide = guide_legend(
          override.aes = list(size = 5)  # 图例里的点变大
        )
      )
  }
  
  # Malignant 面板：用原始 LogP + 单独高值图例
  p_label <- ggplot(plot_label, aes(Factor, Module)) +
    geom_point(aes(size = LogP, colour = Tag), shape = 16) +
    size_label +                              # 单独高值图例（25–30 左右）
    xlab("") + ylab("") +
    theme_base() +
    labs(title = "Malignant") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      plot.background = element_blank(),
      panel.border = element_blank()
    ) +
    scale_colour_manual(
      name   = "Malignant",
      values = c(Significant = "#E53935",
                 notSig     = "#CCCCCC"),
      guide = guide_legend(
        override.aes = list(size = 5)  # 图例里的点变大
      )
    )
  p_label <- p_label + base_big
  p_continuous <- p_continuous + base_big
  p_categorical_other <- p_categorical_other + base_big
  
  # ====================== 合并图表 ==========================
  if (exists("p_continuous") && exists("p_categorical_other") && exists("p_label")) {
    final_plot <- (p_continuous | p_categorical_other | p_label) +
      plot_layout(guides = "collect", widths = c(1.0, 1.2, 1.0)) &
      theme(
        legend.position = "right",
        legend.title    = element_text(size = 16, face = "bold"),
        legend.text     = element_text(size = 15),
        plot.title      = element_text(size = 18, face = "bold", hjust = 0.5)
      )
  } else if (exists("p_continuous") && exists("p_categorical_other")) {
    final_plot <- (p_continuous | p_categorical_other) +
      plot_layout(guides = "collect")
  } else if (exists("p_continuous")) {
    final_plot <- p_continuous
  }
  
  # 保存结果图
  if (exists("final_plot")) {
    ggsave("Figure3D.png",
           final_plot, width = 12, height = 8, dpi = 300)
    cat("图表已保存为: module_analysis_results.png\n")
  }
  
  # 保存数据结果
  write.xlsx(combined_data, "./WGCNA_Patient/module_analysis_results.xlsx")
  cat("结果数据已保存为: module_analysis_results.xlsx\n")
  
  # 打印显著结果
  cat("\n=== 分析结果总结 ===\n")
  if (nrow(plot_info_continuous) > 0) {
    cat("连续变量关联分析:\n")
    significant_continuous <- plot_info_continuous %>%
      filter(P_value < 0.05) %>%
      arrange(P_value)
    print(significant_continuous)
  }
  
  if (nrow(plot_info_categorical) > 0) {
    cat("\n分类变量关联分析:\n")
    significant_categorical <- plot_info_categorical %>%
      filter(P_value < 0.05) %>%
      arrange(P_value)
    print(significant_categorical)
  }
  
} else {
  cat("没有找到有效的分析结果\n")
}

# =============================================================================
# 4. 创建热图分析 (类似原始代码中的分子分类分析)
# =============================================================================

# Label的热图分析
if ("Label" %in% colnames(clinical_filtered)) {
  print("开始创建Label热图分析...")
  
  valid_samples <- clinical_filtered %>%
    filter(!is.na(Label) & Label != "" & !str_detect(toupper(Label), "X")) %>%
    pull(SampleID)
  
  print(paste("Label热图分析有效样本数:", length(valid_samples)))
  
  if (length(valid_samples) > 10) {
    plot_data <- clinical_filtered %>%
      filter(SampleID %in% valid_samples)
    
    module_subset <- module_filtered[, colnames(module_filtered) %in% valid_samples]
    
    print(paste("Label唯一值:", paste(unique(plot_data$Label), collapse = ", ")))
    
    # 转换为长格式数据用于分析
    long_data <- expand.grid(
      Module = module_names,
      Label = unique(plot_data$Label),
      stringsAsFactors = FALSE
    )
    
    # 计算每个模块与每个标签的相关性
    correlation_results <- NULL
    for (module in module_names) {
      for (label in unique(plot_data$Label)) {
        module_scores <- as.numeric(module_subset[module, plot_data$SampleID])
        label_binary <- ifelse(plot_data$Label == label, 1, 0)
        
        if (length(unique(label_binary)) > 1) {
          cor_result <- cor.test(module_scores, label_binary, method = "pearson")
          
          result <- data.frame(
            Module = module,
            Label = label,
            Correlation = cor_result$estimate,
            P_value = cor_result$p.value
          )
          correlation_results <- rbind(correlation_results, result)
        }
      }
    }
    
    if (!is.null(correlation_results) && nrow(correlation_results) > 0) {
      # 添加显著性标记
      correlation_results <- correlation_results %>%
        mutate(Significance = ifelse(P_value < 0.05, "*", ""))
      
      # 设置模块因子顺序
      correlation_results$Module <- factor(correlation_results$Module, 
                                           levels = rev(module_names))
      
      # 创建热图
      p_heatmap <- ggplot(correlation_results, aes(x = Label, y = Module, fill = Correlation)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(
          low = "blue", mid = "white", high = "red", midpoint = 0,
          limits = c(-1, 1), name = "Correlation"
        ) +
        geom_text(aes(label = Significance), color = "black", size = 5) +
        labs(
          title = "Correlation Heatmap: Module vs Label",
          x = "Label",
          y = "Module"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          panel.grid = element_blank()
        )
      
      # 保存热图
      ggsave("./WGCNA_Patient/correlation_heatmap_label.pdf", p_heatmap, width = 8, height = 6, dpi = 300)
      print("标签相关性热图已保存为: correlation_heatmap_label.pdf")
      
      # 保存热图数据
      write.xlsx(correlation_results, "./WGCNA_Patient/correlation_heatmap_data.xlsx")
      print("热图数据已保存为: correlation_heatmap_data.xlsx")
    }
  }
}

print("\n=== 分析完成 ===")
print("已生成的文件:")
print("1. module_analysis_results.pdf - 主要分析结果图")
print("2. module_analysis_results.xlsx - 分析结果数据")
print("3. correlation_heatmap_label.pdf - 标签相关性热图")
print("4. correlation_heatmap_data.xlsx - 热图数据")

