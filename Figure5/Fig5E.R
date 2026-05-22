library(tidyverse)
library(igraph)
library(viridis)
library(patchwork)

## =============================
## 1. 读取数据
## =============================
prot  <- read.csv("Protein.csv", check.names = FALSE)
metab <- read.csv("Metabolism.csv", check.names = FALSE)
cross <- read.csv("pathway_correlation_results_FDR_BH.csv", check.names = FALSE)


cat("📊 Loading data...\n")

## =============================
## 2. 构建网络
## =============================

proteins      <- unique(prot$Gene)
proteins_pref <- paste0("P_", proteins)

prot_pathways     <- unique(cross$Protein)
prot_path_pref    <- paste0("PP_", prot_pathways)


meta_pathways  <- unique(cross$Metabolite)
meta_path_pref    <- paste0("MP_", meta_pathways)

metabolites       <- unique(metab$Metabolite)
metabolites_pref  <- paste0("M_", metabolites)

nodes <- tibble(
  name = c(proteins_pref, prot_path_pref, meta_path_pref, metabolites_pref),
  type = c(
    rep("Protein",            length(proteins_pref)),
    rep("Protein_Pathway",    length(prot_path_pref)),
    rep("Metabolite_Pathway", length(meta_path_pref)),
    rep("Metabolite",         length(metabolites_pref))
  )
)

# 绑定原始名字 + Module / Cluster
nodes <- nodes %>%
  mutate(raw_name = gsub("^(P_|PP_|MP_|M_)", "", name))

# 蛋白的 Module
prot_module <- prot %>%
  distinct(Gene, Module)

# 代谢物的 Cluster
metab_cluster <- metab %>%
  distinct(Metabolite, Cluster)

nodes <- nodes %>%
  left_join(prot_module,  by = c("raw_name" = "Gene")) %>%
  left_join(metab_cluster, by = c("raw_name" = "Metabolite"))

# 构建边
edges_p_to_pp <- prot %>%
  filter(Terms %in% prot_pathways) %>%
  transmute(from = paste0("P_", Gene), to = paste0("PP_", Terms)) %>%
  distinct()

edges_pp_to_mp <- cross %>%
  filter( FDR_BH < 0.001) %>%   # ✅ 只保留显著相关
  transmute(
    from = paste0("PP_", Protein),
    to   = paste0("MP_", Metabolite),
    corr = Correlation,
    pval = FDR_BH
  ) %>%
  distinct(from, to, .keep_all = TRUE)


edges_mp_to_m <- metab %>%
  filter(Terms %in% meta_pathways) %>%
  transmute(from = paste0("MP_", Terms), to = paste0("M_", Metabolite)) %>%
  distinct()

edges_all <- bind_rows(
  edges_p_to_pp %>% mutate(corr = NA_real_, pval = NA_real_),
  edges_pp_to_mp,  # ✅ 自带 corr / pval
  edges_mp_to_m %>% mutate(corr = NA_real_, pval = NA_real_)
) %>%
  filter(from %in% nodes$name & to %in% nodes$name)


g <- graph_from_data_frame(edges_all, vertices = nodes, directed = TRUE)


cat(sprintf("✓ Network: %d nodes, %d edges\n", vcount(g), ecount(g)))

## =============================
## 3. 节点属性
## =============================

degree_df <- tibble(name = V(g)$name, degree = degree(g, mode = "all"))
nodes <- nodes %>% left_join(degree_df, by = "name")

## =============================
## 4. 等高布局（按 Module / Cluster 排序）
## =============================

LAYER_X <- c(
  "Protein" = 2,
  "Protein_Pathway" = 8,
  "Metabolite_Pathway" = 14,
  "Metabolite" = 20
)

TOTAL_HEIGHT <- 20

# 先得到蛋白 Module 和代谢物 Cluster 的排序
prot_module_order <- nodes %>%
  filter(type == "Protein", !is.na(Module)) %>%
  pull(Module) %>%
  unique() %>%
  sort()

metab_cluster_order <- nodes %>%
  filter(type == "Metabolite", !is.na(Cluster)) %>%
  pull(Cluster) %>%
  unique() %>%
  sort()

nodes_layout <- nodes %>%
  mutate(
    x = LAYER_X[type],
    layer_num = match(type, names(LAYER_X)),
    
    # 为 Protein / Metabolite 生成数值排序键
    sort_key = dplyr::case_when(
      type == "Protein"    ~ as.numeric(match(Module,  prot_module_order)),
      type == "Metabolite" ~ as.numeric(match(Cluster, metab_cluster_order)),
      TRUE                 ~ NA_real_
    )
  ) %>%
  group_by(layer_num, type) %>%
  arrange(sort_key, raw_name, .by_group = TRUE) %>%  # 先按分类，再按名字
  mutate(
    n_nodes = n(),
    row_idx = row_number()
  ) %>%
  ungroup() %>%
  mutate(
    spacing = ifelse(n_nodes == 1, 0, TOTAL_HEIGHT / (n_nodes - 1)),
    y = ifelse(
      n_nodes == 1,
      0,
      -(row_idx - 1) * spacing + TOTAL_HEIGHT / 2
    )
  )

nodes_layout <- nodes_layout %>%
  mutate(
    raw_name_wrapped = case_when(
      # 只对包含这个通路名的节点换行
      grepl("Parathyroid hormone synthesis, secretion and action", raw_name) ~
        gsub(",", ",\n", raw_name),
      
      # 如果你希望所有包含逗号的名称都换行，则改成：
      # grepl(",", raw_name) ~ gsub(",", ",\n", raw_name),
      
      TRUE ~ raw_name
    )
  )


## =============================
## 4.1 组信息:Module / Cluster → group
## =============================
wrap_string <- function(x, width = 22){
  sapply(x, function(s){
    paste(strwrap(s, width = width), collapse = "\n")
  })
}

nodes_layout <- nodes_layout %>%
  mutate(
    raw_name_wrapped = wrap_string(raw_name, width = 22)
  )



nodes_layout <- nodes_layout %>%
  mutate(
    group = dplyr::case_when(
      type == "Protein"    & !is.na(Module)  ~ paste0("Module ", Module),
      type == "Metabolite" & !is.na(Cluster) ~ paste0("Cluster ", Cluster),
      type == "Protein_Pathway"              ~ "Protein Pathway",
      type == "Metabolite_Pathway"           ~ "Metabolite Pathway",
      TRUE ~ as.character(type)
    )
  )

## =============================
## 5. 边数据(添加曲线参数)
## =============================

# 这里用 from = Protein Pathway 节点名，to = Metabolite Pathway 节点名
# 注意：要写的是带前缀的名字，比如 "PP_xxx" / "MP_xxx"
key_edges <- tribble(
  ~from,                            ~to,
  "PP_IgSF CAM signaling",         "MP_Phenylalanine metabolism",
  "PP_Sphingolipid signaling pathway", "MP_Sphingolipid signaling pathway"
)

# 
# edge_data <- edges %>%
#   left_join(nodes_layout %>% select(name, x_from = x, y_from = y), 
#             by = c("from" = "name")) %>%
#   left_join(nodes_layout %>% select(name, x_to = x, y_to = y), 
#             by = c("to" = "name")) %>%
#   mutate(
#     edge_type = case_when(
#       x_to == LAYER_X["Protein_Pathway"]    ~ "Type1", # Protein → Protein Pathway
#       x_to == LAYER_X["Metabolite_Pathway"] ~ "Type2", # Protein Pathway → Metabolite Pathway
#       x_to == LAYER_X["Metabolite"]         ~ "Type3"  # Metabolite Pathway → Metabolite
#     ),
#     curvature = (y_to - y_from) / 20 * 0.3
#   ) %>%
#   # 标记关键通路-通路边：只在 Type2（PP→MP）里检测
#   left_join(
#     key_edges %>% mutate(is_key = TRUE),
#     by = c("from", "to")
#   ) %>%
#   mutate(
#     is_key = ifelse(is.na(is_key), FALSE, is_key),
#     edge_group = dplyr::case_when(
#       is_key ~ "Key_Pathway_Cross",   # ⭐ 关键边
#       TRUE   ~ edge_type              # 其他边用原来的 Type1/2/3
#     )
#   )
# 
# edge_data <- edge_data %>%
#   mutate(
#     edge_size = dplyr::case_when(
#       edge_group == "Key_Pathway_Cross" ~ 2.0,   # ⭐ 关键边：更粗
#       TRUE                               ~ 1.2    # 普通边
#     )
#   )
edge_data <- edges_all %>%
  left_join(nodes_layout %>% select(name, x_from = x, y_from = y),
            by = c("from" = "name")) %>%
  left_join(nodes_layout %>% select(name, x_to = x, y_to = y),
            by = c("to" = "name")) %>%
  mutate(
    edge_type = case_when(
      x_to == LAYER_X["Protein_Pathway"]    ~ "Type1", # Protein → Protein Pathway
      x_to == LAYER_X["Metabolite_Pathway"] ~ "Type2", # Protein Pathway → Metabolite Pathway
      x_to == LAYER_X["Metabolite"]         ~ "Type3"  # Metabolite Pathway → Metabolite
    ),
    curvature = (y_to - y_from) / 20 * 0.3
  )

# ✅ 只对 Type2 用相关性决定线宽（P<0.05 已在 edges_pp_to_mp 过滤）
#    线宽映射范围你可以自行调大/调小：to = c(0.6, 3.2)
edge_data <- edge_data %>%
  mutate(
    edge_group = edge_type,
    edge_size = case_when(
      edge_type == "Type2" & !is.na(corr) ~ scales::rescale(abs(corr), to = c(0.6, 3.2)),
      TRUE                                ~ 1.2
    )
  )

## =============================
## 6. 精致配色方案
## =============================

# Module颜色 - 使用科研渐变色(蓝色系)
module_levels  <- nodes_layout %>%
  filter(type == "Protein" & !is.na(Module)) %>%
  pull(Module) %>%
  unique() %>%
  sort()

# Cluster颜色 - 使用暖色系渐变
cluster_levels <- nodes_layout %>%
  filter(type == "Metabolite" & !is.na(Cluster)) %>%
  pull(Cluster) %>%
  unique() %>%
  sort()

# 使用专业配色
module_pal <- if (length(module_levels) > 0) {
  # 蓝色渐变:从浅到深
  cols <- c( "#42A5F5", "#2f44c3","#44c5c3")
  stats::setNames(cols, paste0("Module ", module_levels))
} else {
  character(0)
}

cluster_pal <- if (length(cluster_levels) > 0) {
  # 暖色渐变:橙红色系
  cols <- c("#00BCD4", "#E31A1C")
  stats::setNames(cols, paste0("Cluster ", cluster_levels))
} else {
  character(0)
}

# 通路节点使用优雅的紫色和青色
pathway_cols <- c(
  "Protein Pathway"    = "#7E57C2",    # 深紫色
  "Metabolite Pathway" = "#26A69A"     # 蒂芙尼蓝
)

# 综合颜色
GROUP_COLORS <- c(module_pal, pathway_cols, cluster_pal)

# 背景颜色 - 更柔和
BG_COLORS <- c(
  "Protein" = "#F5F9FF",
  "Protein_Pathway" = "#F3EEF9",
  "Metabolite_Pathway" = "#E8F5F3",
  "Metabolite" = "#FFF3F0"
)

EDGE_COLORS <- c(
  "Type1"            = "#42A5F5",  # Protein → Protein Pathway
  "Type2"            = "#7E57C2",  # Protein Pathway → Metabolite Pathway（普通）
  "Type3"            = "#26A69A",  # Metabolite Pathway → Metabolite
  "Key_Pathway_Cross" = "#e11b1b"  # ⭐ 关键通路-通路交叉
)

## =============================
## 7. 背景矩形(带圆角效果)
## =============================

bg_rects <- nodes_layout %>%
  group_by(type, x) %>%
  summarise(
    xmin = unique(x) - 2.0,
    xmax = unique(x) + 2.0,
    ymin = -TOTAL_HEIGHT / 2 - 1.2,
    ymax = TOTAL_HEIGHT / 2 + 1.2,
    .groups = "drop"
  ) %>%
  mutate(fill_color = BG_COLORS[type])

## =============================
## 8. 层标题
## =============================
layer_titles <- tibble(
  x = LAYER_X,
  y = TOTAL_HEIGHT / 2 + 2.8,
  label = c(
    "Proteins",
    "Protein Pathways",
    "Metabolite Pathways",
    "Metabolites"
  ),
  label_color = c(
    "#42A5F5",                      # Proteins（蓝色）
    pathway_cols["Protein Pathway"],        # Protein Pathway 的颜色（紫色）
    pathway_cols["Metabolite Pathway"],     # Metabolite Pathway 的颜色（青绿）
    "#e11b1b"                       # Metabolites（橙红）
  )
)





## =============================
## 9. 增强绘图
## =============================

# 先算一下节点范围
x_range <- range(nodes_layout$x)
y_range <- range(nodes_layout$y)

p <- ggplot() +
  
  # 1. 背景矩形 - 添加边框
  geom_rect(
    data = bg_rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = bg_rects$fill_color,
    color = "white",
    alpha = 0.4,
    linewidth = 0.8
  ) +
  
  # 2. 边 - 使用更柔和的曲线
  geom_curve(
    data = edge_data,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to, 
        color = edge_group,
        linewidth = edge_size),      # ⭐ 用分组粗细
    curvature = 0.2,
    alpha = 0.6,                     # 可以稍微增加一点透明度，显得更干净
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    lineend = "round"
  ) +
  scale_linewidth_identity() +        # ⭐ 让 linewidth 不被 scale_size 覆盖

  
  # 3. 节点外圈阴影效果
  geom_point(
    data = nodes_layout,
    aes(x = x, y = y, size = degree),
    shape = 21,
    color = "#00000020",
    fill = "#00000010",
    stroke = 0,
    alpha = 0.3
  ) +
  
  # 4. 节点主体
  geom_point(
    data = nodes_layout,
    aes(x = x, y = y, fill = group, size = degree),
    shape = 21,
    color = "white",
    stroke = 3,
    alpha = 0.95
  ) +
  
  # 4a. 非通路节点：圆形 (Protein / Metabolite)
  geom_point(
    data = nodes_layout %>% 
      dplyr::filter(type %in% c("Protein", "Metabolite")),
    aes(x = x, y = y, fill = group, size = degree),
    shape = 21,                # 圆形
    color = "white",
    stroke = 2.2,
    alpha = 0.95
  ) +
  
  # 4b. 通路节点：正方形 (Protein_Pathway / Metabolite_Pathway)
  geom_point(
    data = nodes_layout %>% 
      dplyr::filter(type %in% c("Protein_Pathway", "Metabolite_Pathway")),
    aes(x = x, y = y, fill = group, size = degree),
    shape = 22,                # 正方形
    color = "white",
    stroke = 2.2,
    alpha = 0.95
  ) +
  # 5. 节点标签
  geom_text(
    data = nodes_layout,
    aes(x = x + 0.6, y = y, label = raw_name_wrapped),
    hjust = 0,
    size = 7,
    fontface = "bold",
    color = "#37474F",
    lineheight = 0.9
  ) +
  
  # 6. 层标题 - 使用灰色框
  geom_label(
    data = layer_titles,
    aes(
      x = x,
      y = y,
      label = label,
      colour = I(label_color)   # ⭐ 每个标题使用自己的颜色，不走 scale
    ),
    size = 9,
    fontface = "bold",
    fill = "#FFFFFF",
    alpha = 0.9,
    label.padding = unit(0.35, "lines"),
    label.size = 0.5,
    lineheight = 0.85,
    show.legend = FALSE
  ) +
  

  
  # 7. 主标题
  annotate("text",
           x = mean(LAYER_X),
           y = TOTAL_HEIGHT / 2 + 6.0,
           label = "Protomic-Metabolomics Integrative Network Analysis",
           size = 9,
           fontface = "bold",
           color = "#263238") +
  
  # annotate("text",
  #          x = mean(LAYER_X),
  #          y = TOTAL_HEIGHT / 2 + 5.0,
  #          label = "Protein-Pathway-Metabolite Cross-talk Visualization",
  #          size = 4.5,
  #          color = "#546E7A",
  #          fontface = "italic") +
  
  # # 8. 添加统计信息框
  # annotate("label",
  #          x = 23,
  #          y = -TOTAL_HEIGHT/2 - 0.5,
  #          label = sprintf("Network Statistics\nNodes: %d | Edges: %d", 
  #                          vcount(g), ecount(g)),
  #          size = 3,
  #          hjust = 0,
  #          color = "#37474F",
  #          fill = "#FFFFFF",
  #          alpha = 0.9,
  #          label.size = 0.5) +
  
  # 9. 颜色映射
  scale_fill_manual(
    values = GROUP_COLORS,
    name   = "Node Groups"
  ) +
  
  scale_color_manual(
    values = EDGE_COLORS,
    name = "Connection Types",
    labels = c(
      "Type1" = "Protein → Protein Pathway",
      "Type2" = "Protein Pathway → Metabolite Pathway",
      "Type3" = "Metabolite Pathway → Metabolite",
      "Key_Pathway_Cross" = "Key cross-talk between pathways"
    )
  ) +
  
  scale_size_continuous(
    range = c(8, 15), 
    name = "Node Degree",
    breaks = c(2, 5, 10, 15)
  ) +
  
  # 10. 坐标系统
  coord_cartesian(
    xlim = x_range + c(-2, 5),   # 左右各留一点空白即可
    ylim = y_range + c(-2, 4),   # 上面多留一点给标题
    clip = "off"
  ) +
  # 11. 精致主题
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    
    legend.position = "right",
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.spacing.y = unit(1.4, "cm"),
    
    # 🔺全局字号
    text = element_text(size = 20, family = "sans"),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text  = element_text(size = 20),
    
    legend.key.width  = unit(1.6, "cm"),
    legend.key.height = unit(1.3, "cm"),
    
    plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt")
  ) +
  
  # ========== 图例样式 ==========
guides(
  color = guide_legend(
    title = "Connection Types",
    override.aes = list(linewidth = 2, alpha = 0.9),
    keyheight = unit(1.2, "cm"),  # ⭐ 只给 color legend 增加内部间距
    ncol = 1,
    order = 1          # 第一块：边
  ),
  fill = guide_legend(
    title = "Node Groups",
    override.aes = list(size = 7, stroke = 2),
    keyheight = unit(1.2, "cm"),  # ⭐ 只给 color legend 增加内部间距
    ncol = 1,
    order = 2          # 第二块：节点
  ),
  size = guide_legend(
    title = "Node Degree",
    override.aes = list(fill = "#90CAF9", color = "white", stroke = 2),
    keyheight = unit(1.2, "cm"),  # ⭐ 只给 color legend 增加内部间距
    ncol = 1,
    order = 3          # 第三块：度数
  )
)


## =============================
## 10. 保存高质量图片
## =============================

ggsave("Network_Enhanced_Beautiful.png",
       p, width = 24, height = 14, dpi = 600, bg = "white")

ggsave("Network_Enhanced_Beautiful.pdf",
       p, width = 24, height = 14, device = cairo_pdf, bg = "white")

ggsave("Network_Enhanced_Beautiful.tiff",
       p, width = 24, height = 14, dpi = 600, 
       compression = "lzw", bg = "white")

