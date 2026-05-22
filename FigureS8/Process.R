## =========================================================
## 0. 清空环境 & 加载包
## =========================================================
rm(list = ls())

library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)  # 为了避免版本问题，这里也加载上

## =========================================================
## 1. 读入节点表和边表
## =========================================================
# 节点：Mfuzz_node.csv
nodes <- read.csv("Mfuzz_node.csv", stringsAsFactors = FALSE)

# 边：mFuzz_edge.csv
edges <- read.csv("mFuzz_edge.csv", stringsAsFactors = FALSE)

# 确保关键列为字符型
nodes <- nodes %>%
  mutate(
    ID   = as.character(ID),
    Type1 = as.character(Type1)        # Metabolites / Pathway
  )

edges <- edges %>%
  mutate(
    Source = as.character(Source),
    Target = as.character(Target)
  )

## 检查边表中的节点 ID 是否都在节点表中
missing_in_nodes <- setdiff(unique(c(edges$Source, edges$Target)), nodes$ID)
cat("边表中找不到的节点数:", length(missing_in_nodes), "\n")
if (length(missing_in_nodes) > 0) {
  cat("缺失的节点示例:\n")
  print(head(missing_in_nodes))
}

## =========================================================
## 2. 构建图对象，并把节点属性贴上去
##    这里用 tidygraph::as_tbl_graph + left_join，避免之前的报错
## =========================================================
g_tbl <- as_tbl_graph(edges, directed = FALSE) %>%
  # name 是图中自动生成的节点 ID，与 edges 的 Source/Target 一致
  left_join(nodes, by = c("name" = "ID")) %>%
  mutate(
    # 把 Type1 变成 factor，控制形状次序
    Type1 = factor(Type1, levels = c("Metabolites", "Pathway")),
    # 颜色分组：用 Significance 决定边框颜色
    color_group = case_when(
      Significance == "Up-regulated"   ~ "up",
      Significance == "Down-regulated" ~ "down",
      TRUE                             ~ "ns"
    ),
    color_group = factor(color_group, levels = c("up", "down", "ns")),
    # 先对原始 Log2FC 做裁剪，防止极值
    Log2FC_clip = pmax(pmin(Log2FC, 1), -1)
  )

## =========================================================
## 3. 让「通路节点颜色」由其连接代谢物的 Log2FC 决定
## =========================================================
# 转回 igraph，方便按邻居计算
g_igraph <- as.igraph(g_tbl)

# 取当前每个节点的剪裁后 Log2FC
log2fc_clip <- V(g_igraph)$Log2FC_clip
node_type   <- V(g_igraph)$Type1

# 计算每个节点邻居的平均 Log2FC
avg_neighbor_log2fc <- sapply(V(g_igraph), function(v) {
  nb <- neighbors(g_igraph, v)
  if (length(nb) == 0) return(NA_real_)
  vals <- log2fc_clip[nb]
  if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
})

# 构造用于上色的 Log2FC_pathway：
# - 对 Metabolites：用自己的 Log2FC_clip
# - 对 Pathway：用相邻节点 Log2FC 的平均值，如果没有邻居则置 0
log2fc_pathway <- ifelse(
  node_type == "Pathway",
  ifelse(is.na(avg_neighbor_log2fc), 0, avg_neighbor_log2fc),
  log2fc_clip
)

# 再次裁剪到 [-1, 1]
log2fc_pathway <- pmax(pmin(log2fc_pathway, 1), -1)

# 写回 igraph 对象
V(g_igraph)$Log2FC_pathway <- log2fc_pathway

# 再转回 tidygraph
g_tbl2 <- as_tbl_graph(g_igraph)

## =========================================================
## 4. 绘图：网络 + 标签
##    - 形状：Metabolites = 圆, Pathway = 菱形
##    - 大小：degree
##    - 填充颜色：Log2FC_pathway（Metabolites 自己，Pathway 来自邻居）
##    - 边框颜色： Significance (up/down/ns)
##    - 标签：Label，自动避让
## =========================================================
p <- ggraph(g_tbl2, layout = "fr") +
  # 边
  geom_edge_link(
    colour = "grey80",
    width  = 0.4
  ) +
  # 节点
  geom_node_point(
    aes(
      shape  = Type1,             # 代谢物 / 通路
      size   = centrality_degree(),  # 节点度数
      fill   = Log2FC_pathway,    # 关键：用于上色的数值
      colour = color_group        # 边框颜色：up / down / ns
    ),
    stroke = 0.7
  ) +
  # 标签（使用 geom_node_text + repel=TRUE 自动避让）
  geom_node_text(
    aes(label = Label),
    repel = TRUE,    # ggraph >= 2.0 支持，内部调用 ggrepel
    size  = 3.8
  ) +
  # 形状映射
  scale_shape_manual(
    values = c(
      "Metabolites" = 21,   # 圆形
      "Pathway"     = 23    # 菱形
    ),
    name = "shape"
  ) +
  # 大小映射
  scale_size_continuous(
    range  = c(3, 12),
    breaks = c(2, 4, 6),   # 少量刻度
    name   = "size"
  ) +
  # 填充渐变（可按你喜欢调整颜色）
  scale_fill_gradient2(
    low      = "#009E73",   # 负向（整体 down）
    mid      = "white",
    high     = "#B2182B",   # 正向（整体 up）
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Log2FC"
  ) +
  # 边框颜色：显著性
  scale_colour_manual(
    values = c(
      "down" = "#009E73",
      "up"   = "#B2182B",
      "ns"   = "grey60"
    ),
    name   = "color",
    labels = c("down" = "down", "up" = "up", "ns" = "ns")
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 9)
  )

# 显示
print(p)

# 保存为高分辨率图片
ggsave("FigureS8A.png",
       p, width = 10, height = 8, dpi = 300)
