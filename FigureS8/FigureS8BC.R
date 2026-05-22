# 增强版聚类代谢物可视化验证
# 采用mFuzz风格的美观可视化
# ==========================================

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(scales)

# 1. 定义聚类代谢物
# ==========================================
# 读取 cluster_1 到 cluster_10 的文件，生成 list
cluster_metabolites <- lapply(1:10, function(i) {
  # 文件名
  file_name <- paste0("./mFuzz/cluster_", i, "_metabolites_for_pathway_benign.txt")
  
  # 按行读取代谢物名
  metab <- readLines(file_name, warn = FALSE)
  
  return(metab)
})

# 给 list 命名
names(cluster_metabolites) <- paste0("Cluster_", 1:10)

# 查看结果
str(cluster_metabolites)


# 2. 读取和处理数据
# ==========================================
load_and_process_data <- function() {
  
  # 读取临床数据
  clinical_data <- read.csv("Clinical.csv", stringsAsFactors = FALSE, na.strings = c("", "NA"))
  colnames(clinical_data)[3] <- "Stage"
  
  # 处理临床分期
  clinical_data$Stage_processed <- case_when(
    is.na(clinical_data$Stage) | clinical_data$Stage == "" ~ "Benign",
    clinical_data$Label == "Benign" ~ "Benign",
    
    clinical_data$Stage %in% c(0, "0", "Stage0", "0期") ~ "0",
    
    clinical_data$Stage %in% c("IA", "IB", "I", "StageI", "1期") ~ "I",
    
    clinical_data$Stage %in% c("II", "IIA", "IIB", "2期") ~ "II",
    
    clinical_data$Stage %in% c("III", "IIIA", "IIIB", "IIIC", "3期") ~ "III",
    
    clinical_data$Stage %in% c("IV", "StageIV", "4期") ~ "IV",
    
    TRUE ~ as.character(clinical_data$Stage)
  )
  
  
  # 读取代谢物数据
  metabolite_data <- read.csv("Combined_data_filtered.csv", stringsAsFactors = FALSE, check.names = FALSE)
  
  # 转置代谢物数据
  metabolite_names <- metabolite_data[[1]]
  metabolite_matrix <- as.matrix(metabolite_data[, -1])
  rownames(metabolite_matrix) <- metabolite_names
  metabolite_matrix_t <- t(metabolite_matrix)
  
  # 样本匹配
  common_samples <- intersect(rownames(metabolite_matrix_t), clinical_data$Sample)
  metabolite_final <- metabolite_matrix_t[common_samples, ]
  clinical_final <- clinical_data[clinical_data$Sample %in% common_samples, ]
  clinical_final <- clinical_final[match(common_samples, clinical_final$Sample), ]
  
  # 创建阶段数值映射
  stage_levels <- c("Benign", "0", "I", "II", "III")
  clinical_final$Stage_numeric <- match(clinical_final$Stage_processed, stage_levels)
  
  cat("=== 数据加载完成 ===\n")
  cat("- 样本数:", nrow(metabolite_final), "\n")
  cat("- 代谢物数:", ncol(metabolite_final), "\n")
  cat("- 临床分期分布:\n")
  print(table(clinical_final$Stage_processed))
  
  return(list(
    metabolite_data = metabolite_final,
    clinical_data = clinical_final,
    stage_levels = stage_levels
  ))
}

# 3. 提取聚类代谢物并计算表达轨迹
# ==========================================
extract_cluster_trajectories <- function(data_list) {
  
  metabolite_data <- data_list$metabolite_data
  clinical_data <- data_list$clinical_data
  stage_levels <- data_list$stage_levels
  
  # 创建聚类代谢物映射
  cluster_mapping <- data.frame(
    metabolite = character(),
    cluster = character(),
    stringsAsFactors = FALSE
  )
  
  for(i in 1:length(cluster_metabolites)) {
    cluster_name <- names(cluster_metabolites)[i]
    metabolites <- cluster_metabolites[[i]]
    temp_df <- data.frame(
      metabolite = metabolites,
      cluster = cluster_name,
      stringsAsFactors = FALSE
    )
    cluster_mapping <- rbind(cluster_mapping, temp_df)
  }
  
  # 找到可用的聚类代谢物
  available_metabolites <- intersect(cluster_mapping$metabolite, colnames(metabolite_data))
  cluster_mapping_available <- cluster_mapping[cluster_mapping$metabolite %in% available_metabolites, ]
  
  cat("\n=== 聚类代谢物匹配情况 ===\n")
  cluster_match_summary <- list()
  
  for(cluster_name in names(cluster_metabolites)) {
    cluster_mets <- cluster_metabolites[[cluster_name]]
    available_in_cluster <- intersect(cluster_mets, colnames(metabolite_data))
    cluster_match_summary[[cluster_name]] <- available_in_cluster
    
    cat(sprintf("- %s: %d/%d (%.1f%%)\n", 
                cluster_name, length(available_in_cluster), length(cluster_mets),
                length(available_in_cluster)/length(cluster_mets)*100))
  }
  
  # 计算每个阶段的平均表达
  stage_expression <- list()
  
  for(i in 1:length(stage_levels)) {
    stage <- stage_levels[i]
    stage_samples <- clinical_data$Sample[clinical_data$Stage_processed == stage]
    stage_samples <- intersect(stage_samples, rownames(metabolite_data))
    
    if(length(stage_samples) > 0) {
      if(length(stage_samples) == 1) {
        stage_expr <- metabolite_data[stage_samples, available_metabolites]
      } else {
        stage_expr <- colMeans(metabolite_data[stage_samples, available_metabolites], na.rm = TRUE)
      }
      stage_expression[[i]] <- stage_expr
    } else {
      stage_expression[[i]] <- rep(NA, length(available_metabolites))
    }
  }
  
  # 创建表达矩阵（代谢物 x 阶段）
  expr_matrix <- matrix(NA, nrow = length(available_metabolites), ncol = length(stage_levels))
  rownames(expr_matrix) <- available_metabolites
  colnames(expr_matrix) <- stage_levels
  
  for(i in 1:length(stage_levels)) {
    if(!all(is.na(stage_expression[[i]]))) {
      expr_matrix[, i] <- stage_expression[[i]]
    }
  }
  
  # 标准化表达数据
  expr_matrix_scaled <- t(scale(t(expr_matrix)))
  
  return(list(
    expr_matrix = expr_matrix_scaled,
    cluster_mapping = cluster_mapping_available,
    cluster_match_summary = cluster_match_summary,
    stage_levels = stage_levels,
    available_metabolites = available_metabolites
  ))
}

# 4. 增强版可视化函数
# ==========================================
create_enhanced_cluster_visualization <- function(cluster_results) {
  expr_matrix <- cluster_results$expr_matrix
  cluster_mapping <- cluster_results$cluster_mapping
  stage_levels <- cluster_results$stage_levels
  cluster_match_summary <- cluster_results$cluster_match_summary
  
  # 定义聚类颜色
  cluster_colors <- c("#00BCD4", "#6A3D9A", "#1F78B4", "#33A02C", "#FF7F00", "#B15928", "#FFD700", "#E91E63","#E31A1C",  "#607D8B")
  
  
  # 创建每个聚类的可视化
  cluster_plots <- list()
  
  for(cluster_name in names(cluster_match_summary)) {
    cluster_mets <- cluster_match_summary[[cluster_name]]
    
    if(length(cluster_mets) == 0) next
    
    # 获取聚类索引
    cluster_idx <- as.numeric(gsub("Cluster_", "", cluster_name))
    cluster_color <- cluster_colors[cluster_idx]
    
    # 提取该聚类的表达数据
    cluster_expr <- expr_matrix[cluster_mets, , drop = FALSE]
    
    # 计算聚类中心
    if(nrow(cluster_expr) == 1) {
      cluster_center <- cluster_expr[1, ]
    } else {
      cluster_center <- colMeans(cluster_expr, na.rm = TRUE)
    }
    
    # 准备绘图数据
    plot_data <- data.frame(
      Stage = rep(1:length(stage_levels), nrow(cluster_expr)),
      Expression = as.vector(t(cluster_expr)),
      Metabolite = rep(rownames(cluster_expr), each = length(stage_levels)),
      Stage_Label = rep(stage_levels, nrow(cluster_expr)),
      stringsAsFactors = FALSE
    )
    
    # 计算代谢物重要性评分（基于表达变化幅度）
    metabolite_importance <- sapply(cluster_mets, function(met) {
      met_expr <- cluster_expr[met, ]
      if(all(is.na(met_expr))) return(0)
      return(max(met_expr, na.rm = TRUE) - min(met_expr, na.rm = TRUE))
    })
    
    # 标准化重要性评分到0.3-1范围
    if(max(metabolite_importance) > min(metabolite_importance)) {
      metabolite_importance_norm <- 0.3 + 0.7 * (metabolite_importance - min(metabolite_importance)) / 
        (max(metabolite_importance) - min(metabolite_importance))
    } else {
      metabolite_importance_norm <- rep(0.8, length(metabolite_importance))
    }
    names(metabolite_importance_norm) <- cluster_mets
    
    # 添加重要性评分到数据
    plot_data$Importance <- rep(metabolite_importance_norm[cluster_mets], each = length(stage_levels))
    
    # 聚类中心数据
    center_data <- data.frame(
      Stage = 1:length(stage_levels),
      Expression = cluster_center,
      Stage_Label = stage_levels,
      stringsAsFactors = FALSE
    )
    
    # 创建图表
    p <- ggplot() +
      theme_classic() +
      theme(
        panel.grid.major.y = element_line(color = "grey95", size = 0.3),
        panel.grid.minor.y = element_line(color = "grey98", size = 0.2),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5, 
                                  color = cluster_color),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = "grey80"),
        plot.margin = margin(15, 15, 15, 15)  # 增加右边距以容纳完整代谢物名称
      )
    
    # 添加个体代谢物轨迹
    p <- p + geom_line(data = plot_data, 
                       aes(x = Stage, y = Expression, group = Metabolite,
                           color = Importance, alpha = Importance),
                       size = 1) +
      scale_color_gradient(
        low = alpha(cluster_color, 0.4), 
        high = cluster_color,
        name = "Importance\nScore",
        limits = c(0.3, 1),
        breaks = seq(0.3, 1, by = 0.2),
        labels = function(x) sprintf("%.1f", x)
      ) +
      scale_alpha_continuous(
        range = c(0.4, 0.8),
        limits = c(0.3, 1),
        guide = "none"
      )
    
    # 添加聚类中心线
    p <- p + geom_line(data = center_data,
                       aes(x = Stage, y = Expression),
                       color = "black", size = 3, alpha = 0.8) +
      geom_point(data = center_data,
                 aes(x = Stage, y = Expression),
                 color = "black", size = 4, shape = 19) +
      geom_point(data = center_data,
                 aes(x = Stage, y = Expression),
                 color = cluster_color, size = 3, shape = 19)
    
    # 设置坐标轴
    p <- p + scale_x_continuous(
      breaks = 1:length(stage_levels), 
      labels = stage_levels,
      expand = expansion(mult = c(0.05, 0.35))  # 增加右侧空间给完整代谢物名称
    ) +
      coord_cartesian(xlim = c(0.8, length(stage_levels) + 3)) +  # 增加右侧空间
      labs(
        title = paste(cluster_name, sprintf("(n = %d)", length(cluster_mets))),
        x = "Clinical Stage",
        y = "Standardized Expression"
      )
    
    # 添加代谢物名称标注
    if(length(cluster_mets) > 0) {
      # 选择最重要的代谢物进行标注 (减少数量以避免重叠)
      top_metabolites <- names(sort(metabolite_importance, decreasing = TRUE))
      max_labels <- min(8, length(top_metabolites))  # 减少到8个以避免长名称重叠
      top_metabolites <- top_metabolites[1:max_labels]
      
      # 获取最后阶段的表达值用于标注位置
      last_stage_expr <- cluster_expr[top_metabolites, length(stage_levels)]
      label_data <- data.frame(
        metabolite = top_metabolites,
        y_pos = last_stage_expr,
        importance = metabolite_importance_norm[top_metabolites]
      )
      
      # 排序并分配y位置避免重叠 (增加间距以适应完整名称)
      label_data <- label_data[order(label_data$y_pos, decreasing = TRUE), ]
      y_range <- range(plot_data$Expression, na.rm = TRUE)
      y_spacing <- diff(y_range) / (max_labels + 2)  # 增加间距避免长名称重叠
      y_positions <- seq(y_range[2] - y_spacing/2, y_range[1] + y_spacing/2, length.out = max_labels)
      label_data$y_final <- y_positions
      
      # 添加标注
      for(i in 1:nrow(label_data)) {
        # 计算文本颜色
        text_alpha <- label_data$importance[i]
        text_color <- alpha(cluster_color, text_alpha)
        
        # 保留完整的代谢物名称
        met_name <- label_data$metabolite[i]
        
        p <- p + annotate("text", 
                          x = length(stage_levels) + 0.1, 
                          y = label_data$y_final[i],
                          label = met_name,
                          hjust = 0, vjust = 0.5, 
                          size = 5,  # 稍微减小字体以适应长名称
                          color = text_color,
                          fontface = ifelse(text_alpha > 0.7, "bold", "plain"))
      }
      
      # 如果有更多代谢物未显示
      if(length(cluster_mets) > max_labels) {
        p <- p + annotate("text", 
                          x = length(stage_levels) + 0.1, 
                          y = y_range[1] - 0.1 * diff(y_range),
                          label = paste("... and", length(cluster_mets) - max_labels, "more"),
                          hjust = 0, vjust = 0.5, 
                          size = 3,
                          color = "grey50",
                          fontface = "italic")
      }
    }
    
    cluster_plots[[cluster_name]] <- p
  }
  
  return(cluster_plots)
}

# 5. 主分析函数
# ==========================================
main_analysis <- function() {
  
  cat("=== 增强版聚类代谢物验证分析 ===\n")
  
  # 1. 加载数据
  data_list <- load_and_process_data()
  
  # 2. 提取聚类轨迹
  cluster_results <- extract_cluster_trajectories(data_list)
  
  # 3. 创建增强版可视化
  cat("\n=== 生成增强版可视化 ===\n")
  cluster_plots <- create_enhanced_cluster_visualization(cluster_results)
  
  # 4. 显示所有图表并保存
  cat("\n=== 保存单个聚类图表 ===\n")
  for(cluster_name in names(cluster_plots)) {
    cat(sprintf("显示并保存 %s 图表...\n", cluster_name))
    
    # 显示图表
    print(cluster_plots[[cluster_name]])
    
    # 保存单个图表 (增加宽度以适应完整代谢物名称)
    filename <- paste0("./mFuzz/cluster_validation_", cluster_name, ".png")
    ggsave(filename, plot = cluster_plots[[cluster_name]], 
           width = 10, height = 5, dpi = 300, bg = "white")  # 增加宽度
    cat(sprintf("  已保存: %s\n", filename))
  }
  
  # 5. 生成并保存组合图
  if(length(cluster_plots) >= 2) {
    cat("\n=== 生成并保存组合图表 ===\n")
    
    # 计算网格布局
    n_plots <- length(cluster_plots)
    n_cols <- ifelse(n_plots <= 4, 2, 2)
    n_rows <- ceiling(n_plots / n_cols)
    
    # 为组合图调整每个子图的边距（减少右边距）
    cluster_plots_combined <- lapply(cluster_plots, function(p) {
      p + theme(plot.margin = margin(10, 30, 10, 10))  # 减少右边距
    })
    
    # 创建组合图
    combined_plot <- grid.arrange(
      grobs = cluster_plots, 
      ncol = n_cols,
      top = grid::textGrob("mFuzz-Style Cluster Validation Across Clinical Stages", 
                           gp = grid::gpar(fontsize = 16, fontface = "bold")),
      widths = rep(1, n_cols),  # 统一子图宽度
      heights = rep(1, n_rows),  # 统一子图高度
      padding = unit(0.5, "cm")  # 减少子图间距
    )
    
    print(combined_plot)
    
    # 保存组合图 (增加宽度以适应完整代谢物名称)
    combined_filename <- "FigureS8BC.png"
    ggsave(combined_filename, plot = combined_plot, 
           width = n_cols * 10, height = n_rows * 5, dpi = 300, bg = "white")  # 增加宽度
    cat(sprintf("已保存组合图: %s\n", combined_filename))
    
    # 同时保存PDF版本（用于出版）
    # pdf_filename <- "cluster_validation_combined.pdf"
    # pdf(pdf_filename, width = n_cols * 10, height = n_rows * 5)  # 增加宽度
    # print(combined_plot)
    # dev.off()
    # cat(sprintf("已保存PDF版本: %s\n", pdf_filename))
  }
  
  # 6. 总结
  cat("\n=== 验证总结 ===\n")
  total_metabolites <- length(cluster_results$available_metabolites)
  cat("- 共验证了", total_metabolites, "个聚类代谢物\n")
  
  for(cluster_name in names(cluster_results$cluster_match_summary)) {
    count <- length(cluster_results$cluster_match_summary[[cluster_name]])
    cat(sprintf("  %s: %d个代谢物\n", cluster_name, count))
  }
  
  cat("- 临床分期: Benign → 0 → I → II\n")
  cat("- 可视化样式: mFuzz增强版，显示表达轨迹和聚类中心\n")
  
  cat("\n=== 保存的文件 ===\n")
  for(cluster_name in names(cluster_plots)) {
    cat(sprintf("- %s.png\n", paste0("./mFuzz/cluster_validation_", cluster_name)))
  }
  if(length(cluster_plots) >= 2) {
    cat("- cluster_validation_combined.png (高分辨率PNG)\n")
    cat("- cluster_validation_combined.pdf (出版质量PDF)\n")
  }
  cat("所有图片保存在当前工作目录中\n")
  
  return(list(
    data = data_list,
    cluster_results = cluster_results,
    plots = cluster_plots,
    saved_files = c(
      paste0("./mFuzz/cluster_validation_", names(cluster_plots), ".png"),
      if(length(cluster_plots) >= 2) c("cluster_validation_combined.png", "cluster_validation_combined.pdf")
    )
  ))
}

# 6. 运行分析
# ==========================================
cat("增强版代码加载完成！\n")
cat("运行分析请执行: results <- main_analysis()\n")

results <- main_analysis()