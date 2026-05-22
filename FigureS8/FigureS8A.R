## ================================
## 1. 读入通路注释与定义背景集
## ================================

# 读入 All_pathway.csv
pathway_map <- read.csv("All_pathway.csv",
                        stringsAsFactors = FALSE,
                        check.names = FALSE)

# 看一下列名，确认代谢物列名字
colnames(pathway_map)
# 假定代谢物列名为 "Metabolite name"
# 通路列为 "Terms"，以及 "Pathway ID", "Pathway level 1", "Pathway level 2" 等

## =========================================
## 2. 定义函数：给某个 cluster 做 代谢物–通路 映射
## =========================================
library(dplyr)

map_cluster_metabolites <- function(cluster_idx,
                                    pathway_map,
                                    file_pattern = "cluster_%d_metabolites_for_pathway_benign.txt") {
  # 文件名
  fn <- sprintf(file_pattern, cluster_idx)
  if (!file.exists(fn)) {
    warning("文件不存在: ", fn, "，跳过。")
    return(NULL)
  }
  
  # 读入该 cluster 的代谢物列表
  mets <- readLines(fn)
  mets <- trimws(mets)
  mets <- mets[mets != ""]
  
  if (length(mets) == 0) {
    warning("Cluster ", cluster_idx, " 没有读到任何代谢物，跳过。")
    return(NULL)
  }
  
  cat("Cluster", cluster_idx, "代谢物数:", length(mets), "\n")
  
  # 做一个 data.frame：Cluster + Metabolite
  df_cluster <- data.frame(
    Cluster   = paste0("Cluster_", cluster_idx, "_benign"),
    Metabolite = mets,
    stringsAsFactors = FALSE
  )
  
  # 与 All_pathway 做左连接：
  # Metabolite 对应 All_pathway$`Metabolite name`
  merged <- df_cluster %>%
    left_join(pathway_map,
              by = c("Metabolite" = "Metabolite name"))
  
  # 你可以选择：
  # 1）保留所有代谢物，即使没有通路（Terms 为 NA）
  # 2）只保留有通路注释的行：
  # merged <- merged %>% filter(!is.na(Terms))
  
  return(merged)
}

## =========================================
## 3. 对 10 个 cluster 循环映射
## =========================================
all_list <- list()

for (i in 1:10) {
  res_i <- map_cluster_metabolites(
    cluster_idx = i,
    pathway_map = pathway_map
  )
  if (!is.null(res_i)) {
    all_list[[paste0("Cluster_", i)]] <- res_i
  }
}

# 合并所有 cluster
met_path_all <- dplyr::bind_rows(all_list)

# 看一下结果
head(met_path_all)

## =========================================
## 4. 导出为一个 CSV 文件
## =========================================
write.csv(
  met_path_all,
  file = "Metabolite_Pathway_mapping_by_cluster.csv",
  row.names = FALSE
)

cat("已导出：Metabolite_Pathway_mapping_by_cluster.csv\n")
