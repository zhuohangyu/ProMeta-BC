# 清理工作空间
rm(list = ls())
gc()


# 加载必要的包
library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# 基于分析结果创建数据框，按频率降序排列以匹配图例
metabolite_data <- data.frame(
  Category1 = c("Lipid", "Amino Acid", "Xenobiotics", "Unclassified", "Peptide", 
                "Carbohydrate", "Nucleotide", "Cofactors and Vitamins", "Energy", "Secondary metabolites"),
  Freq = c(268, 134, 90, 34, 33, 27, 23, 19, 1, 1),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Freq)) # 确保按频率降序排列

# 创建对应的图例标签（按数据顺序）
legend_labels <- c("lipid (43.4%)", "amino acid (19.6%)", "xenobiotics (14.6%)", 
                   "unclassified (5.5%)", "peptide (5.3%)", "carbohydrate (4.4%)",
                   "nucleotide (3.7%)", "cofactors and vitamins (3.1%)", 
                   "energy (0.2%)", "secondary metabolites (0.2%)")
# 计算百分比
metabolite_data <- metabolite_data %>%
  mutate(
    percentage = round((Freq / sum(Freq)) * 100, 1)
  )

# 总代谢物数量
total_metabolites <- sum(metabolite_data$Freq)

# 排序
metabolite_data <- metabolite_data %>% arrange(desc(Freq))

# 设置颜色并命名
colors <- c( "#DDA0DD", "#87CEEB", "#F0E68C", "#98FB98", "#FFB6C1",
            "#B0C4DE", "#20B2AA", "#90EE90", "#F5DEB3", "#D3D3D3")
names(colors) <- metabolite_data$Category1

legend_labels <- paste0(tolower(metabolite_data$Category1), " (", metabolite_data$percentage, "%)")
# 使用基础R创建带有彩色图例的饼状图
create_pie_with_colored_legend <- function() {
  par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
  
  # 饼图
  pie(metabolite_data$Freq, 
      labels = "", 
      col = colors,  # 确保和全局定义的 colors 对应
      border = "black",
      main = "Tissue samples",
      cex.main = 1.2)
  
  # 图例
  par(mar = c(1, 0, 3, 1))
  plot.new()
  title(main = paste0("Metabolite class (n = ", total_metabolites, ")"), 
        cex.main = 1.0, font.main = 2)
  
  y_positions <- seq(0.9, 0.1, length.out = length(legend_labels))
  
  for(i in 1:length(legend_labels)) {
    rect(0.05, y_positions[i] - 0.02, 0.15, y_positions[i] + 0.02, 
         col = colors[i], border = "black")
    text(0.2, y_positions[i], legend_labels[i], 
         col = colors[i], adj = 0, cex = 0.9, font = 1)
  }
  
  rect(0.02, 0.08, 0.98, 0.95, col = NA, border = "black", lwd = 2)
}
# 创建饼状图
pie_chart <- ggplot(metabolite_data, aes(x = "", y = Freq, fill = Category1)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +
  coord_polar("y", start = 0) +
  scale_fill_manual(
    values = colors,
    labels = legend_labels
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 24),
    legend.background = element_rect(color = "black", fill = "white", size = 1),
    legend.margin = margin(30, 30, 30, 30),
    legend.box.margin = margin(10, 10, 10, 20),
    plot.margin = margin(30, 30, 30, 30)
  ) +
  labs(
    title = "Plasma samples",
    fill = paste0("Metabolite class (n = ", total_metabolites, ")")
  ) +
  guides(fill = guide_legend(
    override.aes = list(
      color = colors,  # 图例项边框颜色
      fill = colors    # 图例项填充颜色
    )
  )) +
  # 手动设置图例文字颜色
  theme(
    legend.text = element_text(size = 22, color = "black")  # 先设为黑色，然后用自定义方法
  )
# 保存标准ggplot版本
ggsave("Figure4A.png", pie_chart, width = 10, height = 6, dpi = 300)