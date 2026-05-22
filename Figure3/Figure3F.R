# 清理工作空间
rm(list = ls())
gc()

# Mfuzz Analysis - Clustering of Differential Proteins
# 完整修改版本：包含Benign阶段处理的蛋白组学分析 (benign-0-I-II-III-IV)
# ================================================================

# Set global random seed for reproducibility
RANDOM_SEED <- 42  # Change this number to get different but reproducible results
set.seed(RANDOM_SEED)
cat("=== Proteomics Mfuzz Analysis with Random Seed:", RANDOM_SEED, "===\n")

# 1. Load Required Packages
# ================================================================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("Mfuzz", "Biobase", "dplyr", "ggplot2", 
                       "limma", "pheatmap", "RColorBrewer", "reshape2", "gridExtra")

for (pkg in required_packages) { 
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% c("Mfuzz", "Biobase", "limma")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}
library(viridis)
library(scales)
library(grid)

# 2. Data Loading and Preprocessing (修改版：蛋白组学数据处理)
# ================================================================
# Read data files
clinical_data <- read.csv("Clinical.csv", stringsAsFactors = FALSE)
proteomics_data <- read.csv("Combined_data_filtered.csv", stringsAsFactors = FALSE)

# Check the actual column names
cat("Clinical data columns:\n")
print(colnames(clinical_data))
cat("\nProteomics data columns (first 10):\n")
print(colnames(proteomics_data)[1:10])

# Rename clinical stage column (it's actually column 3, "Clinical stage")
colnames(clinical_data)[3] <- "Clinical_Stage"

# Set up proteomics matrix - use Gene name as identifier, or Protein accession if Gene name is empty
# First, create a unique identifier for each protein
protein_ids <- ifelse(is.na(proteomics_data$Gene.name) | proteomics_data$Gene.name == "", 
                      proteomics_data$Protein.accession, 
                      proteomics_data$Gene.name)

# Handle duplicates by adding suffix
protein_ids <- make.unique(protein_ids)

# Set row names and extract expression matrix (excluding first 3 annotation columns)
rownames(proteomics_data) <- protein_ids
proteomics_matrix <- proteomics_data[, -c(1:3)]  # Remove first 3 columns (annotation)

# Convert to numeric matrix
proteomics_matrix <- apply(proteomics_matrix, 2, function(x) as.numeric(as.character(x)))
rownames(proteomics_matrix) <- protein_ids

# Check original clinical stages BEFORE cleaning
cat("\nOriginal clinical stages found:", paste(unique(clinical_data$Clinical_Stage), collapse = ", "), "\n")
cat("Number of NA values in Clinical_Stage:", sum(is.na(clinical_data$Clinical_Stage)), "\n")
cat("Number of empty strings in Clinical_Stage:", sum(clinical_data$Clinical_Stage == "", na.rm = TRUE), "\n")

# ================================================================
# NEW: Handle NA and empty values by replacing with "benign"
# ================================================================
cat("\n=== BENIGN STAGE PROCESSING ===\n")

# Replace NA values with "benign"
clinical_data$Clinical_Stage[is.na(clinical_data$Clinical_Stage)] <- "benign"

# Replace empty strings with "benign"
clinical_data$Clinical_Stage[clinical_data$Clinical_Stage == ""] <- "benign"

# Replace any variation of "normal", "control", "healthy" with "benign" (case insensitive)
clinical_data$Clinical_Stage[grepl("^normal$|^control$|^healthy$|^benign$", 
                                   clinical_data$Clinical_Stage, ignore.case = TRUE)] <- "benign"

cat("After replacing NA/empty with 'benign':", paste(unique(clinical_data$Clinical_Stage), collapse = ", "), "\n")

# ================================================================
# Clean and standardize clinical stages: remove A, B, C suffixes and standardize to benign, 0, I, II, III, IV
# ================================================================
clinical_data$Clinical_Stage_Clean <- clinical_data$Clinical_Stage

# Remove any A, B, C suffixes (case insensitive) - but preserve "benign"
clinical_data$Clinical_Stage_Clean[clinical_data$Clinical_Stage_Clean != "benign"] <- 
  gsub("[ABC]$", "", clinical_data$Clinical_Stage_Clean[clinical_data$Clinical_Stage_Clean != "benign"], ignore.case = TRUE)

# Standardize stage names (keep "benign" as is)
clinical_data$Clinical_Stage_Clean <- gsub("^0$", "0", clinical_data$Clinical_Stage_Clean)
clinical_data$Clinical_Stage_Clean <- gsub("^1$|^I$", "I", clinical_data$Clinical_Stage_Clean)
clinical_data$Clinical_Stage_Clean <- gsub("^2$|^II$", "II", clinical_data$Clinical_Stage_Clean)
clinical_data$Clinical_Stage_Clean <- gsub("^3$|^III$", "III", clinical_data$Clinical_Stage_Clean)
clinical_data$Clinical_Stage_Clean <- gsub("^4$|^IV$", "IV", clinical_data$Clinical_Stage_Clean)

# Update the main Clinical_Stage column
clinical_data$Clinical_Stage <- clinical_data$Clinical_Stage_Clean

# Check cleaned stages
unique_stages_clean <- unique(clinical_data$Clinical_Stage)
cat("\nCleaned clinical stages:", paste(unique_stages_clean, collapse = ", "), "\n")

# ================================================================
# NEW: Create ordered factor with benign as first level, followed by cancer stages
# ================================================================
stage_levels <- c("benign", "0", "I", "II", "III", "IV")
# Only use stages that exist in the data
existing_stages <- intersect(stage_levels, unique_stages_clean)
cat("\nUsing stage levels in order:", paste(existing_stages, collapse = " → "), "\n")

clinical_data$Clinical_Stage <- factor(clinical_data$Clinical_Stage, 
                                       levels = existing_stages,
                                       ordered = TRUE)

# Match samples between datasets
common_samples <- intersect(colnames(proteomics_matrix), clinical_data$Sample)
cat("\nNumber of common samples:", length(common_samples), "\n")
cat("Total samples in proteomics:", ncol(proteomics_matrix), "\n")
cat("Total samples in clinical:", nrow(clinical_data), "\n")

# Filter to common samples only - ensure exact matching
proteomics_matrix <- proteomics_matrix[, common_samples]
clinical_data_matched <- clinical_data[match(common_samples, clinical_data$Sample), ]

# Remove any NA rows and ensure no missing data
clinical_data_matched <- clinical_data_matched[!is.na(clinical_data_matched$Sample), ]
valid_samples <- clinical_data_matched$Sample

# Final filtering to ensure perfect alignment
proteomics_matrix <- proteomics_matrix[, valid_samples]
clinical_data <- clinical_data_matched[clinical_data_matched$Sample %in% valid_samples, ]

# Ensure the order is exactly the same
clinical_data <- clinical_data[match(colnames(proteomics_matrix), clinical_data$Sample), ]

# Final verification
cat("\nAfter matching:")
cat("\n  Proteomics samples:", ncol(proteomics_matrix))
cat("\n  Clinical samples:", nrow(clinical_data))
cat("\n  Sample names match:", all(colnames(proteomics_matrix) == clinical_data$Sample))

# Remove any samples with missing stage information after cleaning
valid_stage_samples <- !is.na(clinical_data$Clinical_Stage)
clinical_data <- clinical_data[valid_stage_samples, ]

# Re-align proteomics data after stage cleaning
common_samples_final <- intersect(colnames(proteomics_matrix), clinical_data$Sample)
proteomics_matrix <- proteomics_matrix[, common_samples_final]
clinical_data <- clinical_data[clinical_data$Sample %in% common_samples_final, ]

# Ensure exact order matching
clinical_data <- clinical_data[match(colnames(proteomics_matrix), clinical_data$Sample), ]

cat("\nAfter stage cleaning:")
cat("\n  Final common samples:", length(common_samples_final))
cat("\n  Proteomics samples:", ncol(proteomics_matrix))
cat("\n  Clinical samples:", nrow(clinical_data))

cat("Data dimensions:")
cat("\nProteomics matrix:", dim(proteomics_matrix))
cat("\nClinical data:", dim(clinical_data))

# ================================================================
# Display stage distribution with BENIGN stage
# ================================================================
cat("\n=== FINAL STAGE DISTRIBUTION (including BENIGN) ===\n")
stage_table <- table(clinical_data$Clinical_Stage)
print(stage_table)

# Create a more detailed summary
stage_summary <- data.frame(
  Stage = names(stage_table),
  Count = as.numeric(stage_table),
  Percentage = round(as.numeric(stage_table) / sum(stage_table) * 100, 1),
  stringsAsFactors = FALSE
)
stage_summary$Description <- c(
  "benign" = "Benign/Normal/Control samples",
  "0" = "Stage 0 (Carcinoma in situ)",
  "I" = "Stage I (Early-stage localized cancer)", 
  "II" = "Stage II (Locally advanced cancer)",
  "III" = "Stage III (Regional lymph node involvement)",
  "IV" = "Stage IV (Distant metastasis)"
)[stage_summary$Stage]

cat("\nDetailed stage distribution:\n")
print(stage_summary)

# Check if we have enough samples per stage
min_samples_per_stage <- 3
stages_with_few_samples <- names(stage_table)[stage_table < min_samples_per_stage]
if (length(stages_with_few_samples) > 0) {
  cat("\nWarning: The following stages have fewer than", min_samples_per_stage, "samples:\n")
  print(stage_table[stages_with_few_samples])
  cat("Consider combining stages or using a different grouping strategy.\n")
}

# Show progression from benign to advanced cancer
cat("\n=== CANCER PROGRESSION ANALYSIS ===\n")
cat("Sample progression from benign to advanced stages:\n")
for (i in 1:length(existing_stages)) {
  stage <- existing_stages[i]
  count <- stage_table[stage]
  cat(sprintf("  %s: %d samples (%.1f%%)\n", 
              stage, count, count/sum(stage_table)*100))
}

# 3. Differential Analysis to Identify Significant Proteins
# ================================================================
# Convert proteomics matrix to numeric (in case there are character values)
# Also handle any non-numeric values
proteomics_numeric <- apply(proteomics_matrix, 2, function(x) as.numeric(as.character(x)))
rownames(proteomics_numeric) <- rownames(proteomics_matrix)

# Ensure sample order matches between proteomics and clinical data
sample_order <- colnames(proteomics_numeric)
clinical_data <- clinical_data[match(sample_order, clinical_data$Sample), ]

# Verify perfect alignment before proceeding
if (!all(colnames(proteomics_numeric) == clinical_data$Sample)) {
  stop("Sample names do not match between proteomics and clinical data!")
}

cat("\nFinal data verification:")
cat("\n  Proteomics dimensions:", dim(proteomics_numeric))
cat("\n  Clinical data dimensions:", dim(clinical_data))
cat("\n  Perfect sample alignment:", all(colnames(proteomics_numeric) == clinical_data$Sample))

# Remove proteins with too many missing values
missing_threshold <- 0.3  # Remove if >30% missing
valid_proteins <- rowSums(is.na(proteomics_numeric)) < (ncol(proteomics_numeric) * missing_threshold)
proteomics_filtered <- proteomics_numeric[valid_proteins, ]

cat("\nProteins after missing value filtering:", nrow(proteomics_filtered))

# Create design matrix for limma - now includes benign as reference
design <- model.matrix(~ Clinical_Stage, data = clinical_data)
cat("\nDesign matrix dimensions:", dim(design))
cat("\nDesign matrix columns:", colnames(design), "\n")

# Note: "benign" will be the reference level (intercept) in the design matrix
# Other stages will be compared against benign

# Double-check dimensions before fitting
cat("\nDimension check before lmFit:")
cat("\n  proteomics_filtered:", dim(proteomics_filtered))
cat("\n  design matrix:", dim(design))
cat("\n  Samples match:", ncol(proteomics_filtered) == nrow(design))

# Fit linear model with error checking
if (ncol(proteomics_filtered) != nrow(design)) {
  cat("\nERROR: Dimension mismatch!")
  cat("\n  proteomics_filtered samples:", colnames(proteomics_filtered)[1:min(5, ncol(proteomics_filtered))])
  cat("\n  clinical_data samples:", clinical_data$Sample[1:min(5, nrow(clinical_data))])
  cat("\n  design matrix rows:", nrow(design))
  stop("Cannot proceed with lmFit due to dimension mismatch.")
}

fit <- lmFit(proteomics_filtered, design)
fit <- eBayes(fit)

# Extract differential proteins for each stage comparison (vs benign)
coef_names <- colnames(design)[-1]  # Remove intercept (benign reference)
cat("\nCoefficients for comparison (vs benign):", coef_names, "\n")

# Collect significant proteins from all comparisons
significant_proteins <- c()
p_threshold <- 0.05
logfc_threshold <- 0.5

for (coef in coef_names) {
  if (coef %in% colnames(fit$coefficients)) {
    results <- topTable(fit, coef = coef, number = Inf, sort.by = "P")
    sig_proteins <- rownames(results)[results$P.Value < p_threshold & 
                                        abs(results$logFC) > logfc_threshold]
    significant_proteins <- c(significant_proteins, sig_proteins)
    cat("Significant proteins for", coef, "vs benign:", length(sig_proteins), "\n")
  }
}

# Remove duplicates
significant_proteins <- unique(significant_proteins)
cat("\nTotal unique significant proteins found:", length(significant_proteins))

# If too few significant proteins, use top variable proteins instead
if (length(significant_proteins) < 20) {
  cat("\nToo few significant proteins, using top variable proteins instead\n")
  protein_vars <- apply(proteomics_filtered, 1, var, na.rm = TRUE)
  top_var_proteins <- names(sort(protein_vars, decreasing = TRUE)[1:min(100, nrow(proteomics_filtered))])
  significant_proteins <- top_var_proteins
  cat("Using top", length(significant_proteins), "variable proteins\n")
}

# 4. Create Stage-wise Mean Expression for Significant Proteins (INCLUDING BENIGN)
# ================================================================
# Filter matrix to significant proteins only
sig_proteomics <- proteomics_filtered[significant_proteins, ]

# Calculate mean expression for each stage (including benign)
stage_means <- data.frame(row.names = rownames(sig_proteomics))

for (stage in levels(clinical_data$Clinical_Stage)) {
  stage_samples <- clinical_data$Sample[clinical_data$Clinical_Stage == stage]
  # Only use samples that exist in both datasets
  stage_samples <- intersect(stage_samples, colnames(sig_proteomics))
  if (length(stage_samples) > 0) {
    stage_data <- sig_proteomics[, stage_samples, drop = FALSE]
    # Use standardized column names: Stage_benign, Stage_0, Stage_I, etc.
    stage_means[[paste0("Stage_", stage)]] <- rowMeans(stage_data, na.rm = TRUE)
    cat("Stage", stage, "- samples:", length(stage_samples), "\n")
  } else {
    cat("Warning: No samples found for stage", stage, "\n")
  }
}

# Remove any remaining NA values
stage_means_clean <- stage_means[complete.cases(stage_means), ]
cat("\nFinal proteins for clustering:", nrow(stage_means_clean))

# Show the stage means structure
cat("\nStage means data structure (including benign):\n")
print(head(stage_means_clean))

# 5. Create ExpressionSet for Mfuzz (INCLUDING BENIGN)
# ================================================================
# Create phenotype data matching the stage means columns
stage_levels <- levels(clinical_data$Clinical_Stage)
phenoData <- data.frame(
  Stage = stage_levels,
  Time = seq_along(stage_levels) - 1,  # 0-based numbering (benign=0, 0=1, I=2, etc.)
  row.names = paste0("Stage_", stage_levels)
)

cat("\nPhenotype data (including benign as baseline):\n")
print(phenoData)

# Create ExpressionSet
eset <- new("ExpressionSet",
            exprs = as.matrix(stage_means_clean),
            phenoData = new("AnnotatedDataFrame", data = phenoData))

cat("\nExpressionSet created with dimensions:", dim(exprs(eset)))

# Print progression information
cat("\n=== PROGRESSION MODELING ===\n")
cat("Cancer progression modeled as:\n")
for (i in 1:length(stage_levels)) {
  cat(sprintf("  Time %d: Stage %s\n", i-1, stage_levels[i]))
}

# 6. Mfuzz Analysis
# ================================================================
# Set random seed for reproducible results
set.seed(RANDOM_SEED)  # Use the global seed variable

# Fill missing values and standardize
eset_filled <- fill.NA(eset)
eset_std <- standardise(eset_filled)

# Estimate fuzzification parameter
m_est <- mestimate(eset_std)
cat("\nEstimated fuzzification parameter m:", m_est)

# Set seed again before cluster validation
set.seed(RANDOM_SEED)

# Determine optimal cluster number using different criteria
cluster_range <- 2:8
cluster_validity <- numeric(length(cluster_range))

for (i in seq_along(cluster_range)) {
  k <- cluster_range[i]
  set.seed(RANDOM_SEED)  # Use the global seed for each clustering attempt
  cl_temp <- mfuzz(eset_std, c = k, m = m_est)
  cluster_validity[i] <- cl_temp$withinerror
}

# Plot cluster validation
pdf("./mFuzz/proteomics_cluster_validation.pdf", width = 6, height = 6)
plot(cluster_range, cluster_validity, type = "b", pch = 19,
     xlab = "Number of clusters", ylab = "Within-cluster error",
     main = "Protein Cluster Number Optimization")
grid()
dev.off()

# Select optimal cluster number
optimal_k <- 2  # You can change this or use: cluster_range[which.min(cluster_validity)]
cat("\nUsing", optimal_k, "clusters")

# 7. Final Clustering
# ================================================================
# Set seed again before final clustering for reproducible results
set.seed(RANDOM_SEED)

cl <- mfuzz(eset_std, c = optimal_k, m = m_est)

cat("\nClustering completed with seed", RANDOM_SEED, "for reproducible results")

# 8. Create Publication-Style Visualization
# ================================================================
create_publication_mfuzz_plot <- function(eset_std, cl, membership_threshold = 0.2) {
  
  expr_data <- exprs(eset_std)
  n_clusters <- ncol(cl$membership)
  
  # Set up colors for clusters
  cluster_colors <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928", "#FFD700", "#E91E63", "#00BCD4", "#607D8B")
  
  # Create PDF with proper layout
  pdf("./mFuzz/proteomics_clustering_differential_proteins.pdf", width = 14, height = 16)
  
  # Set up multi-panel layout
  par(mfrow = c(3, 2), mar = c(4, 4, 3, 10), oma = c(2, 2, 4, 2))
  
  for (k in 1:n_clusters) {
    # Get proteins assigned to this cluster with sufficient membership
    cluster_assignment <- apply(cl$membership, 1, which.max)
    max_membership <- apply(cl$membership, 1, max)
    
    cluster_members <- which(cluster_assignment == k & max_membership >= membership_threshold)
    
    if (length(cluster_members) == 0) next
    
    # Set y-axis limits based on data range
    y_range <- range(expr_data, na.rm = TRUE)
    
    # Create base plot
    n_stages <- ncol(expr_data)
    plot(1:n_stages, cl$centers[k, ], type = "n", 
         xlim = c(0.8, n_stages + 1.4), ylim = y_range,
         xlab = "Clinical Stage", ylab = "Relative abundance",
         main = paste("Cluster", k),
         xaxt = "n", cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1)
    
    # Add x-axis labels using actual stage names
    stage_labels <- levels(clinical_data$Clinical_Stage)
    axis(1, at = 1:n_stages, labels = stage_labels, cex.axis = 1.1)
    
    # Add light grid
    abline(h = seq(from = floor(y_range[1]*2)/2, to = ceiling(y_range[2]*2)/2, by = 0.5), 
           col = "lightgray", lty = 1, lwd = 0.5)
    abline(v = 1:n_stages, col = "lightgray", lty = 1, lwd = 0.5)
    
    # Plot individual protein trajectories
    colors_alpha <- adjustcolor(cluster_colors[k], alpha.f = 0.6)
    
    for (i in seq_along(cluster_members)) {
      protein_idx <- cluster_members[i]
      lines(1:n_stages, expr_data[protein_idx, ], 
            col = colors_alpha, lwd = 1.2)
    }
    
    # Plot cluster center as thick black line
    lines(1:n_stages, cl$centers[k, ], col = "black", lwd = 3)
    
    # Add protein names on the right side
    protein_names <- rownames(expr_data)[cluster_members]
    
    # Limit displayed names to prevent overcrowding
    max_names <- min(15, length(protein_names))
    selected_names <- protein_names[1:max_names]
    
    # Calculate positions for text
    text_y_start <- y_range[2] - 0.1 * diff(y_range)
    text_y_step <- diff(y_range) * 0.8 / max_names
    
    for (i in seq_along(selected_names)) {
      text_y <- text_y_start - (i-1) * text_y_step
      text(n_stages + 0.15, text_y, selected_names[i], 
           pos = 4, cex = 0.8, col = "black")
    }
    
    # Add count information
    mtext(paste0("n = ", length(cluster_members)), side = 3, line = 0.5, 
          cex = 0.9, col = "darkgray")
  }
  
  # Add main title
  mtext("Clustering of differential proteins (Benign as baseline)", outer = TRUE, cex = 1.6, font = 2, line = 1)
  
  dev.off()
}

# Generate the publication plot
create_publication_mfuzz_plot(eset_std, cl, membership_threshold = 0.2)

# 9. Enhanced ggplot2 Visualization
# ================================================================
create_enhanced_ggplot_version <- function(eset_std, cl, membership_threshold = 0.2, color_max = NULL) {
  
  expr_data <- exprs(eset_std)
  n_clusters <- ncol(cl$membership)
  
  # Get cluster assignments and membership values
  cluster_assignment <- apply(cl$membership, 1, which.max)
  max_membership <- apply(cl$membership, 1, max)
  
  # Calculate dynamic color range if not provided
  if (is.null(color_max)) {
    overall_max_membership <- max(max_membership, na.rm = TRUE)
    color_max <- ceiling(overall_max_membership * 10) / 10
  }
  
  cat("Dynamic color range: ", membership_threshold, " to ", color_max, "\n")
  
  # Prepare data for ggplot
  plot_data_list <- list()
  
  for (k in 1:n_clusters) {
    cluster_members <- which(cluster_assignment == k & max_membership >= membership_threshold)
    
    if (length(cluster_members) == 0) next
    
    # Individual protein data with membership values
    n_stages <- ncol(expr_data)
    protein_data <- data.frame(
      Stage = rep(1:n_stages, length(cluster_members)),
      Expression = as.vector(t(expr_data[cluster_members, ])),
      Protein = rep(rownames(expr_data)[cluster_members], each = n_stages),
      Membership = rep(max_membership[cluster_members], each = n_stages),
      Cluster = k,
      Type = "Individual"
    )
    
    # Cluster center data
    center_data <- data.frame(
      Stage = 1:n_stages,
      Expression = cl$centers[k, ],
      Protein = "Center",
      Membership = color_max,
      Cluster = k,
      Type = "Center"
    )
    
    plot_data_list[[k]] <- rbind(protein_data, center_data)
  }
  
  # Combine all data
  plot_data <- do.call(rbind, plot_data_list)
  plot_data$Stage_Label <- factor(plot_data$Stage, levels = 1:ncol(expr_data), 
                                  labels = levels(clinical_data$Clinical_Stage))
  
  # Define color palette for clusters
  cluster_colors <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928", "#FFD700", "#E91E63", "#00BCD4", "#607D8B")
  
  # Create plots for each cluster
  cluster_plots <- list()
  
  for (k in 1:n_clusters) {
    cluster_data <- plot_data[plot_data$Cluster == k, ]
    individual_data <- cluster_data[cluster_data$Type == "Individual", ]
    center_data <- cluster_data[cluster_data$Type == "Center", ]
    
    if (nrow(individual_data) == 0) next
    
    # Get protein count
    n_proteins <- length(unique(individual_data$Protein))
    
    # Create the plot
    p <- ggplot() +
      theme_classic() +
      theme(
        panel.grid.major.y = element_line(color = "grey95", size = 0.3),
        panel.grid.minor.y = element_line(color = "grey98", size = 0.2),
        plot.title = element_text(size = 36, face = "bold", hjust = 0.5, 
                                  color = cluster_colors[k]),
        axis.title = element_text(size = 36, face = "bold"),
        axis.text = element_text(size = 32),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 32, face = "bold", margin = margin(b = 15)),
        legend.text = element_text(size = 32),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = "grey80"),
        strip.background = element_rect(fill = "grey95"),
        plot.margin = margin(10, 10, 10, 10),
        legend.key.height = unit(1.4, "cm")
      )
    
    # Add individual protein lines
    p <- p + geom_line(data = individual_data, 
                       aes(x = Stage, y = Expression, group = Protein,
                           color = Membership, alpha = Membership),
                       size = 1.2) +
      scale_color_gradient(
        low = alpha(cluster_colors[k], 0.3), 
        high = cluster_colors[k],
        name = "Membership\nScore",
        limits = c(membership_threshold, color_max),
        breaks = seq(membership_threshold, color_max, by = 0.1),
        labels = function(x) sprintf("%.1f", x)
      ) +
      scale_alpha_continuous(
        range = c(0.4, 0.9),
        limits = c(membership_threshold, color_max),
        guide = "none"
      )
    
    # Add cluster center as thick line
    p <- p + geom_line(data = center_data,
                       aes(x = Stage, y = Expression, group = 1),
                       color = "black", size = 5, alpha = 0.9) +
      geom_point(data = center_data,
                 aes(x = Stage, y = Expression),
                 color = "black", size = 7, shape = 19) +
      geom_point(data = center_data,
                 aes(x = Stage, y = Expression),
                 color = cluster_colors[k], size = 5.5, shape = 19)
    
    # Set scales and labels
    p <- p + scale_x_continuous(
      breaks = 1:ncol(expr_data), 
      labels = levels(clinical_data$Clinical_Stage),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
      coord_cartesian(xlim = c(0.8, ncol(expr_data) + 1.5)) +
      labs(
        title = paste("Cluster", k, sprintf("(n = %d)", n_proteins)),
        x = "Clinical Stage",
        y = "Relative Abundance"
      )
    
    # Add protein names
    # Add protein names (Top 10) + "Top 10" label
    if (nrow(individual_data) > 0) {
      protein_info <- individual_data %>%
        dplyr::filter(Stage == ncol(expr_data)) %>%
        dplyr::arrange(desc(Membership)) %>%
        dplyr::select(Protein, Expression, Membership) %>%
        dplyr::distinct()
      
      max_names <- min(15, nrow(protein_info))   # ✅ Top15
      protein_info <- protein_info[1:max_names, , drop = FALSE]
      
      # ✅ y_range 必须定义
      y_range <- range(c(individual_data$Expression, center_data$Expression), na.rm = TRUE)
      
      # ===== 留白参数（可调）=====
      top_reserved    <- 0.12 * diff(y_range)  # 顶部留白：给“Top 10”标题
      bottom_reserved <- 0.18 * diff(y_range)  # ✅底部留白：给“…and more”等文字
      
      name_top    <- y_range[2] - top_reserved
      name_bottom <- y_range[1] + bottom_reserved
      
      # Top10 名字只排布在 [name_bottom, name_top] 之间
      y_positions <- seq(from = name_top, to = name_bottom, length.out = max_names)
      
      
      # ✅ 标题：Top 10（放在图内，避免看不见）
      p <- p + annotate(
        "text",
        x = ncol(expr_data) + 0.2,
        y = y_range[2] - 0.02 * diff(y_range),
        label = "Top 15",
        hjust = 0, vjust = 1,
        size = 9,
        color = "grey20",
        fontface = "bold"
      )
      
      
      # 右侧蛋白名
      for (i in 1:nrow(protein_info)) {
        text_color <- grDevices::colorRamp(c(alpha(cluster_colors[k], 0.5), cluster_colors[k]))(
          (protein_info$Membership[i] - membership_threshold) / (color_max - membership_threshold)
        )
        text_color_hex <- grDevices::rgb(text_color[1], text_color[2], text_color[3], maxColorValue = 255)
        
        p <- p + annotate(
          "text",
          x = ncol(expr_data) + 0.2,
          y = y_positions[i],
          label = protein_info$Protein[i],
          hjust = 0, vjust = 0.5,
          size = 9,
          color = text_color_hex,
          fontface = ifelse(
            protein_info$Membership[i] > (membership_threshold + color_max) / 2,
            "bold", "plain"
          )
        )
      }
      
      # 可选：更多提示
      n_proteins <- length(unique(individual_data$Protein))
      if (n_proteins > max_names) {
        bottom_text_y <- y_range[1] + 0.06 * diff(y_range)  # 在底部保留区里（偏上）
        p <- p + annotate(
          "text",
          x = ncol(expr_data) - 0.5,
          y = bottom_text_y,
          label = paste("... and", n_proteins - max_names, "more"),
          hjust = 0, vjust = 0,
          size = 6.8,
          color = "grey50",
          fontface = "italic"
        )
        
      }
      
      # ✅ 防止右侧文字被 panel 裁剪
      p <- p + coord_cartesian(
        xlim = c(0.8, ncol(expr_data) + 1.5),
        clip = "off"
      )
    }
    
    
    cluster_plots[[k]] <- p
  }
  
  # Remove NULL plots
  cluster_plots <- cluster_plots[!sapply(cluster_plots, is.null)]
  
  # Create enhanced combined plot
  png("Figure3F.png",
      width = 25, height = 10, units = "in", res = 300)
  
  combined_plot <- grid.arrange(
    grobs = cluster_plots,
    ncol = 2
  )
  
  print(combined_plot)
  dev.off()
  
  
  return(list(plots = cluster_plots, color_max = color_max))
}

# Calculate color range
cluster_assignment <- apply(cl$membership, 1, which.max)
max_membership <- apply(cl$membership, 1, max)
overall_max_membership <- max(max_membership, na.rm = TRUE)
color_max <- ceiling(overall_max_membership * 10) / 10

# Generate enhanced visualization
enhanced_result <- create_enhanced_ggplot_version(eset_std, cl, membership_threshold = 0.2, color_max = color_max)

# 10. Results Export and Analysis
# ================================================================
# Create detailed results
membership_threshold <- 0.2
cluster_assignment <- apply(cl$membership, 1, which.max)
max_membership <- apply(cl$membership, 1, max)

results_summary <- data.frame(
  Protein = rownames(cl$membership),
  Cluster = cluster_assignment,
  Max_Membership = max_membership,
  High_Confidence = max_membership >= membership_threshold,
  stringsAsFactors = FALSE
)

# Add expression pattern description
for (i in 1:nrow(results_summary)) {
  cluster_id <- results_summary$Cluster[i]
  center_pattern <- cl$centers[cluster_id, ]
  
  # Classify pattern based on trajectory (now relative to benign)
  if (center_pattern[ncol(cl$centers)] > center_pattern[1] + 0.2) {
    if (length(center_pattern) > 2 && center_pattern[2] > center_pattern[1] + 0.1) {
      pattern <- "Early and sustained upregulation"
    } else {
      pattern <- "Late cancer-stage upregulation"
    }
  } else if (center_pattern[ncol(cl$centers)] < center_pattern[1] - 0.2) {
    if (length(center_pattern) > 2 && center_pattern[2] < center_pattern[1] - 0.1) {
      pattern <- "Early and sustained downregulation"
    } else {
      pattern <- "Late cancer-stage downregulation"
    }
  } else {
    # Check for intermediate patterns
    max_diff <- max(center_pattern) - min(center_pattern)
    if (max_diff > 0.3) {
      mid_point <- length(center_pattern) %/% 2
      if (center_pattern[mid_point] > center_pattern[1] && center_pattern[mid_point] > center_pattern[length(center_pattern)]) {
        pattern <- "Mid-stage peak expression"
      } else {
        pattern <- "Variable across cancer stages"
      }
    } else {
      pattern <- "Stable relative to benign"
    }
  }
  
  results_summary$Pattern[i] <- pattern
}

# Export results
write.csv(results_summary, "./mFuzz/mfuzz_differential_proteins_results_with_benign.csv", row.names = FALSE)

# Export cluster-specific results
for (k in 1:optimal_k) {
  cluster_proteins <- results_summary[results_summary$Cluster == k & 
                                        results_summary$High_Confidence, ]
  
  if (nrow(cluster_proteins) > 0) {
    write.csv(cluster_proteins, 
              paste0("./mFuzz/cluster_", k, "_differential_proteins_benign_baseline.csv"), 
              row.names = FALSE)
    
    writeLines(cluster_proteins$Protein, 
               paste0("./mFuzz/cluster_", k, "_proteins_for_pathway_benign.txt"))
  }
}

# 11. Reference Image Style Individual Cluster Plots
# ================================================================
# Simplified cluster plotting function - Reference image style
plot_cluster_proteins_simple <- function(cluster_num, mfuzz_result, standardized_data, 
                                         min_membership = 0.2, max_proteins = 20, 
                                         stage_labels = NULL, show_center = FALSE,
                                         plot_width = 10, plot_height = 8,
                                         base_font_size = 16) {
  
  # Get core proteins for this cluster
  cluster_assignment <- apply(mfuzz_result$membership, 1, which.max)
  max_membership <- apply(mfuzz_result$membership, 1, max)
  
  # Get proteins assigned to this cluster with sufficient membership
  cluster_members <- which(cluster_assignment == cluster_num & max_membership >= min_membership)
  
  if (length(cluster_members) == 0) {
    warning(paste("Cluster", cluster_num, "has no proteins with membership >=", min_membership))
    return(NULL)
  }
  
  # Get protein names
  core_proteins <- rownames(standardized_data)[cluster_members]
  
  # If too many proteins, select top ones by membership
  if (length(core_proteins) > max_proteins) {
    membership_values <- max_membership[cluster_members]
    top_indices <- order(membership_values, decreasing = TRUE)[1:max_proteins]
    core_proteins <- core_proteins[top_indices]
    cluster_members <- cluster_members[top_indices]
  }
  
  # Extract standardized expression data for these proteins
  expr_data <- standardized_data[core_proteins, , drop = FALSE]
  
  # Prepare plotting data
  plot_data <- data.frame()
  
  for (i in 1:nrow(expr_data)) {
    protein <- rownames(expr_data)[i]
    
    protein_data <- data.frame(
      Protein = protein,
      Stage = colnames(expr_data),
      Expression = as.numeric(expr_data[i, ]),
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, protein_data)
  }
  
  # Set stage labels
  if (!is.null(stage_labels)) {
    plot_data$Stage <- factor(plot_data$Stage, 
                              levels = colnames(expr_data), 
                              labels = stage_labels)
  } else {
    plot_data$Stage <- factor(plot_data$Stage, levels = colnames(expr_data))
  }
  
  # Generate colors similar to reference image
  n_proteins <- length(core_proteins)
  # Use professional color palette similar to the reference image
  reference_colors <- c(
    "#FF6B9D", "#C44569", "#F8B500", "#FF7043", "#FFA726", "#FFCC02", 
    "#8BC34A", "#66BB6A", "#26A69A", "#4DB6AC", "#26C6DA", "#29B6F6",
    "#42A5F5", "#5C6BC0", "#7986CB", "#9575CD", "#BA68C8", "#E91E63",
    "#FF5722", "#795548", "#78909C", "#90A4AE", "#A1887F", "#BCAAA4"
  )
  
  # Select colors based on number of proteins
  if (n_proteins <= length(reference_colors)) {
    selected_colors <- reference_colors[1:n_proteins]
  } else {
    # If more proteins than predefined colors, interpolate additional colors
    selected_colors <- colorRampPalette(reference_colors)(n_proteins)
  }
  names(selected_colors) <- core_proteins
  
  # Create plot similar to reference image
  p <- ggplot(plot_data, aes(x = Stage, y = Expression, group = Protein, color = Protein)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 2.5, alpha = 0.9) +
    
    scale_color_manual(values = selected_colors, name = paste("Cluster", cluster_num)) +
    
    labs(
      y = "Log2 (Relative abundance)",  # Using standardized data
      x = "Clinical Stage (Benign → Cancer Progression)"
    ) +
    
    # Clean theme similar to reference image with larger fonts throughout
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_line(color = "grey95", size = 0.3),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      
      # All axis elements scaled up
      axis.title.x = element_text(size = base_font_size + 2, color = "black", face = "bold"),
      axis.title.y = element_text(size = base_font_size + 2, color = "black", face = "bold"),
      axis.text = element_text(size = base_font_size, color = "black", face = "bold"),
      axis.text.x = element_text(size = base_font_size, color = "black", face = "bold"),
      axis.text.y = element_text(size = base_font_size, color = "black", face = "bold"),
      
      # Legend elements scaled up
      legend.title = element_text(size = base_font_size + 4, face = "bold"),
      legend.text = element_text(size = base_font_size + 2, face = "bold"),  # Large protein names
      legend.position = "right",
      legend.key.size = unit(1.0, "cm"),  # Larger legend keys
      legend.key.width = unit(1.5, "cm"), # Wider legend keys for better line visibility
      legend.background = element_rect(fill = "white", color = "black", size = 1),
      legend.margin = margin(l = 15, r = 15, t = 15, b = 15),  # More margin around legend
      
      # Plot elements
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Adjust plot margins to accommodate larger fonts
      plot.margin = margin(t = 15, r = 20, b = 15, l = 15)
    )
  
  # Optionally add cluster center line
  if (show_center) {
    center_data <- data.frame(
      Stage = if(!is.null(stage_labels)) stage_labels else levels(plot_data$Stage),
      Expression = mfuzz_result$centers[cluster_num, ]
    )
    
    p <- p + geom_line(data = center_data, 
                       aes(x = Stage, y = Expression, group = 1), 
                       color = "black", size = 3, linetype = "dashed", 
                       inherit.aes = FALSE) +
      geom_point(data = center_data,
                 aes(x = Stage, y = Expression),
                 color = "black", size = 4, shape = 19,
                 inherit.aes = FALSE)
  }
  
  return(p)
}

# Generate reference image style individual cluster plots
# ================================================================
# Create directory for cluster plots
dir.create("proteomics_cluster_plots_benign_baseline", showWarnings = FALSE)

# Define stage labels using the actual clinical stages from your data
stage_labels_simple <- levels(clinical_data$Clinical_Stage)

# Generate simple style plots for each cluster
cluster_plots_simple <- list()
for (i in 1:optimal_k) {
  cat(paste("\nGenerating benign baseline plot for cluster", i, "..."))
  
  p <- plot_cluster_proteins_simple(
    cluster_num = i,
    mfuzz_result = cl,
    standardized_data = exprs(eset_std),  # Use standardized data
    min_membership = 0.2,
    max_proteins = 20,
    stage_labels = stage_labels_simple,
    show_center = FALSE,  # Set to TRUE to show cluster center line
    plot_width = 10,      # Narrower width
    plot_height = 8,      # Standard height
    base_font_size = 14   # Large font size for ALL text elements
  )
  
  if (!is.null(p)) {
    cluster_plots_simple[[i]] <- p
    
    # Save individual plots (narrower width with larger legend text)
    output_file_pdf <- paste0("proteomics_cluster_plots_benign_baseline/cluster_", i, "_proteins_benign.pdf")
    output_file_png <- paste0("proteomics_cluster_plots_benign_baseline/cluster_", i, "_proteins_benign.png")
    
    # Use the plot dimensions from the function parameters
    plot_width <- 12   # Can be customized in function call above
    plot_height <- 10   # Can be customized in function call above
    
    ggsave(output_file_pdf, p, width = plot_width, height = plot_height)
    ggsave(output_file_png, p, width = plot_width, height = plot_height, dpi = 300)
    
    cat(" - Saved")
  }
}

# Create combined plot
# ================================================================
# Remove any NULL plots
cluster_plots_simple <- cluster_plots_simple[!sapply(cluster_plots_simple, is.null)]

if (length(cluster_plots_simple) > 0) {
  # Create combined plot with narrower width
  pdf("./mFuzz/combined_proteomics_cluster_plots_benign_baseline.pdf", width = 18, height = 16)  # Reduced width from 20 to 18
  
  combined_plot_simple <- grid.arrange(
    grobs = cluster_plots_simple, 
    ncol = 2,
    top = textGrob("Protein Clustering: Benign → Cancer Progression Patterns", 
                   gp = gpar(fontsize = 18, fontface = "bold")),
    bottom = textGrob(paste("Membership threshold: 0.2 | Fuzzification parameter m:", round(m_est, 3),
                            "| Clinical stages:", paste(levels(clinical_data$Clinical_Stage), collapse = " → "),
                            "| Random seed:", RANDOM_SEED, "| Baseline: Benign"),
                      gp = gpar(fontsize = 12, col = "grey40"))
  )
  
  print(combined_plot_simple)
  dev.off()
}

cat("\nStage sample counts:\n")
stage_counts <- table(clinical_data$Clinical_Stage)
for (i in 1:length(stage_counts)) {
  cat(sprintf("  Stage %s: %d samples\n", names(stage_counts)[i], stage_counts[i]))
}

cat("\nCluster Summary (patterns relative to benign):\n")
for (k in 1:optimal_k) {
  cluster_proteins <- sum(results_summary$Cluster == k & results_summary$High_Confidence)
  if (cluster_proteins > 0) {
    avg_membership <- mean(results_summary$Max_Membership[results_summary$Cluster == k])
    patterns <- unique(results_summary$Pattern[results_summary$Cluster == k])
    pattern <- paste(patterns, collapse = ", ")
    
    cat(sprintf("Cluster %d: %d proteins (avg membership: %.3f) - %s\n", 
                k, cluster_proteins, avg_membership, pattern))
  }
}

# 13. Save Analysis Results (Updated for Benign Integration)
# ================================================================
# 保存原有变量 + 新增变量
save(cl, eset_std, m_est, optimal_k, 
     clinical_data, proteomics_filtered, significant_proteins, stage_means_clean, RANDOM_SEED,
     file = "./mFuzz/proteomics_mfuzz_analysis_complete_with_benign_baseline.RData")

cat("✓ 已保存完整的蛋白组学Mfuzz分析结果（Benign基线）到 proteomics_mfuzz_analysis_complete_with_benign_baseline.RData\n")

# 创建整合分析变量 (针对Benign基线)
# ================================================================
# 1. 创建 mfuzz_clusters - Mfuzz聚类结果
mfuzz_clusters <- data.frame(
  Protein = rownames(cl$membership),
  Cluster = apply(cl$membership, 1, which.max),
  Max_Membership = apply(cl$membership, 1, max),
  stringsAsFactors = FALSE
)

# 添加聚类标签（例如：Cluster_1, Cluster_2等）
mfuzz_clusters$Cluster_Label <- paste0("Cluster_", mfuzz_clusters$Cluster)

# 添加高置信度标记
membership_threshold <- 0.2
mfuzz_clusters$High_Confidence <- mfuzz_clusters$Max_Membership >= membership_threshold

# 添加所有聚类的成员概率
membership_df <- as.data.frame(cl$membership)
colnames(membership_df) <- paste0("Membership_Cluster_", 1:ncol(membership_df))
mfuzz_clusters <- cbind(mfuzz_clusters, membership_df)

# 2. 创建 significant_results - 差异分析结果
significant_results <- data.frame(
  Feature = rownames(fit$coefficients),
  stringsAsFactors = FALSE
)

# 添加各个比较的统计结果 (vs benign)
coef_names <- colnames(design)[-1]  # 移除截距项
for (coef in coef_names) {
  if (coef %in% colnames(fit$coefficients)) {
    results <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
    
    # 添加统计结果
    significant_results[[paste0(coef, "_vs_benign_logFC")]] <- results$logFC
    significant_results[[paste0(coef, "_vs_benign_pvalue")]] <- results$P.Value
    significant_results[[paste0(coef, "_vs_benign_padj")]] <- results$adj.P.Val
    significant_results[[paste0(coef, "_vs_benign_significant")]] <- results$adj.P.Val < 0.05 & abs(results$logFC) > 0.1
  }
}

# 创建综合显著性标记
p_cols <- grep("_vs_benign_pvalue$", colnames(significant_results), value = TRUE)
significant_results$Any_Significant_vs_Benign <- rowSums(significant_results[, gsub("_pvalue", "_significant", p_cols)], na.rm = TRUE) > 0

# 3. 创建聚类间相关性矩阵
cluster_centers <- cl$centers
correlation_matrix <- cor(t(cluster_centers), use = "complete.obs")
rownames(correlation_matrix) <- paste0("Cluster_", 1:nrow(correlation_matrix))
colnames(correlation_matrix) <- paste0("Cluster_", 1:ncol(correlation_matrix))

# 4. 创建调整后的p值矩阵（基于聚类特征，相对于benign）
p_adjusted <- matrix(NA, nrow = optimal_k, ncol = length(levels(clinical_data$Clinical_Stage)) - 1)
rownames(p_adjusted) <- paste0("Cluster_", 1:optimal_k)
colnames(p_adjusted) <- paste0("Stage_", levels(clinical_data$Clinical_Stage)[-1], "_vs_benign")

# 5. 创建目标配对信息（相对于benign的关联）
target_pair <- data.frame(
  Source_Type = character(0),
  Source_ID = character(0),
  Target_Type = character(0),
  Target_ID = character(0),
  Relationship = character(0),
  Strength = numeric(0),
  P_value = numeric(0),
  stringsAsFactors = FALSE
)

# 添加聚类间的显著关联
for (i in 1:(nrow(correlation_matrix)-1)) {
  for (j in (i+1):ncol(correlation_matrix)) {
    cor_val <- correlation_matrix[i, j]
    if (abs(cor_val) > 0.5) {  # 显著相关阈值
      target_pair <- rbind(target_pair, data.frame(
        Source_Type = "Mfuzz_Cluster",
        Source_ID = rownames(correlation_matrix)[i],
        Target_Type = "Mfuzz_Cluster", 
        Target_ID = colnames(correlation_matrix)[j],
        Relationship = ifelse(cor_val > 0, "Positive_Correlation", "Negative_Correlation"),
        Strength = abs(cor_val),
        P_value = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 保存整合分析所需的核心变量（针对Benign基线）
save(mfuzz_clusters, significant_results, correlation_matrix, p_adjusted, target_pair,
     file = "./mFuzz/proteomics_mfuzz_integration_variables_benign_baseline.RData")


# 14. Optional: 合并高度相似的两个聚类，仅用于可视化展示
# ================================================================
cat("\n\n=== Step 14: 检查并合并相似聚类（仅用于可视化） ===\n")

# 仅在 Mfuzz 的聚类数为 2 时考虑合并
if (ncol(cl$membership) == 2) {
  center_cor <- cor(cl$centers[1, ], cl$centers[2, ], use = "complete.obs")
  cat("Cluster 1 与 Cluster 2 中心轨迹的相关系数:", round(center_cor, 3), "\n")
  
  # 相关性阈值可以自行调整，例如 0.9 / 0.95
  merge_threshold <- 0.9
  
  if (center_cor > merge_threshold) {
    cat(">>> 两个聚类中心轨迹高度相似，将合并为一个 MergedCluster 进行额外可视化展示。\n")
    
    # 构建一个“仅用于画图”的合并版本 cl_plot
    cl_plot <- cl
    
    # 合并中心轨迹（取两条中心的平均）
    merged_center <- colMeans(cl$centers)
    cl_plot$centers <- matrix(merged_center, nrow = 1)
    rownames(cl_plot$centers) <- "MergedCluster"
    
    # 对于 c = 2，membership_1 + membership_2 = 1，所以合并后 membership 可以设为 1
    merged_membership <- matrix(1, nrow = nrow(cl$membership), ncol = 1)
    rownames(merged_membership) <- rownames(cl$membership)
    colnames(merged_membership) <- "MergedCluster"
    cl_plot$membership <- merged_membership
    
    # -----------------------------
    # 14.1 合并后的 Mfuzz 轨迹图（Publication-style 版本）
    #    使用新的文件名，不覆盖原始结果
    # -----------------------------
    
    create_publication_mfuzz_plot_merged <- function(eset_std, cl, membership_threshold = 0.2) {
      expr_data <- exprs(eset_std)
      n_clusters <- ncol(cl$membership)
      
      # 颜色（与前面风格一致）
      cluster_colors <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", 
                          "#B15928", "#FFD700", "#E91E63", "#00BCD4", "#607D8B")
      
      # 新文件名：_merged 版本
      pdf("./mFuzz/proteomics_clustering_differential_proteins_merged.pdf", 
          width = 10, height = 8)
      
      par(mfrow = c(1, n_clusters), mar = c(4, 4, 3, 10), oma = c(2, 2, 4, 2))
      
      for (k in 1:n_clusters) {
        cluster_assignment <- apply(cl$membership, 1, which.max)
        max_membership <- apply(cl$membership, 1, max)
        cluster_members <- which(cluster_assignment == k & max_membership >= membership_threshold)
        if (length(cluster_members) == 0) next
        
        y_range <- range(expr_data, na.rm = TRUE)
        n_stages <- ncol(expr_data)
        
        plot(1:n_stages, cl$centers[k, ], type = "n",
             xlim = c(0.8, n_stages + 1.4), ylim = y_range,
             xlab = "Clinical Stage", ylab = "Relative abundance",
             main = paste("Merged Cluster", k),
             xaxt = "n", cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1)
        
        stage_labels <- levels(clinical_data$Clinical_Stage)
        axis(1, at = 1:n_stages, labels = stage_labels, cex.axis = 1.1)
        
        abline(h = seq(from = floor(y_range[1]*2)/2, 
                       to = ceiling(y_range[2]*2)/2, by = 0.5),
               col = "lightgray", lty = 1, lwd = 0.5)
        abline(v = 1:n_stages, col = "lightgray", lty = 1, lwd = 0.5)
        
        colors_alpha <- adjustcolor(cluster_colors[k], alpha.f = 0.6)
        for (i in seq_along(cluster_members)) {
          protein_idx <- cluster_members[i]
          lines(1:n_stages, expr_data[protein_idx, ],
                col = colors_alpha, lwd = 1.2)
        }
        
        lines(1:n_stages, cl$centers[k, ], col = "black", lwd = 3)
        
        protein_names <- rownames(expr_data)[cluster_members]
        max_names <- min(15, length(protein_names))
        selected_names <- protein_names[1:max_names]
        
        text_y_start <- y_range[2] - 0.1 * diff(y_range)
        text_y_step <- diff(y_range) * 0.8 / max_names
        
        for (i in seq_along(selected_names)) {
          text_y <- text_y_start - (i-1) * text_y_step
          text(n_stages + 0.15, text_y, selected_names[i],
               pos = 4, cex = 0.8, col = "black")
        }
        
        mtext(paste0("n = ", length(cluster_members)), side = 3, line = 0.5,
              cex = 0.9, col = "darkgray")
      }
      
      mtext("Merged clustering of differential proteins (Benign as baseline)",
            outer = TRUE, cex = 1.6, font = 2, line = 1)
      dev.off()
    }
    
    # 调用合并版绘图函数
    create_publication_mfuzz_plot_merged(eset_std, cl_plot, membership_threshold = 0.2)
    cat(">>> 已生成合并后的 Mfuzz 轨迹图: ./mFuzz/proteomics_clustering_differential_proteins_merged.pdf\n")
    
    
    # -----------------------------
    # 14.2 合并后的 ggplot 加强版可视化
    # -----------------------------
    create_enhanced_ggplot_version_merged <- function(eset_std, cl, membership_threshold = 0.2,
                                                      color_max = NULL) {
      expr_data <- exprs(eset_std)
      n_clusters <- ncol(cl$membership)
      
      cluster_assignment <- apply(cl$membership, 1, which.max)
      max_membership <- apply(cl$membership, 1, max)
      
      if (is.null(color_max)) {
        overall_max_membership <- max(max_membership, na.rm = TRUE)
        color_max <- ceiling(overall_max_membership * 10) / 10
      }
      
      cat("Merged version color range: ", membership_threshold, " to ", color_max, "\n")
      
      plot_data_list <- list()
      for (k in 1:n_clusters) {
        cluster_members <- which(cluster_assignment == k & max_membership >= membership_threshold)
        if (length(cluster_members) == 0) next
        
        n_stages <- ncol(expr_data)
        protein_data <- data.frame(
          Stage = rep(1:n_stages, length(cluster_members)),
          Expression = as.vector(t(expr_data[cluster_members, ])),
          Protein = rep(rownames(expr_data)[cluster_members], each = n_stages),
          Membership = rep(max_membership[cluster_members], each = n_stages),
          Cluster = k,
          Type = "Individual"
        )
        center_data <- data.frame(
          Stage = 1:n_stages,
          Expression = cl$centers[k, ],
          Protein = "Center",
          Membership = color_max,
          Cluster = k,
          Type = "Center"
        )
        plot_data_list[[k]] <- rbind(protein_data, center_data)
      }
      
      plot_data <- do.call(rbind, plot_data_list)
      plot_data$Stage_Label <- factor(
        plot_data$Stage,
        levels = 1:ncol(expr_data),
        labels = levels(clinical_data$Clinical_Stage)
      )
      
      cluster_colors <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A",
                          "#B15928", "#FFD700", "#E91E63", "#00BCD4", "#607D8B")
      
      cluster_plots <- list()
      for (k in 1:n_clusters) {
        cluster_data <- plot_data[plot_data$Cluster == k, ]
        individual_data <- cluster_data[cluster_data$Type == "Individual", ]
        center_data <- cluster_data[cluster_data$Type == "Center", ]
        if (nrow(individual_data) == 0) next
        
        n_proteins <- length(unique(individual_data$Protein))
        
        p <- ggplot() +
          theme_classic() +
          theme(
            panel.grid.major.y = element_line(color = "grey95", size = 0.3),
            panel.grid.minor.y = element_line(color = "grey98", size = 0.2),
            plot.title = element_text(size = 30, face = "bold", hjust = 0.5,
                                      color = cluster_colors[k]),
            axis.title = element_text(size = 30, face = "bold"),
            axis.text = element_text(size = 28),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            legend.title = element_text(size = 28, face = "bold", margin = margin(b = 15)),
            legend.text = element_text(size = 28),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = "grey80"),
            strip.background = element_rect(fill = "grey95"),
            plot.margin = margin(10, 10, 10, 10),
            legend.key.height = unit(1.4, "cm")
          ) +
          geom_line(data = individual_data,
                    aes(x = Stage, y = Expression, group = Protein,
                        color = Membership, alpha = Membership),
                    size = 1.2) +
          scale_color_gradient(
            low = alpha(cluster_colors[k], 0.3),
            high = cluster_colors[k],
            name = "Membership\nScore",
            limits = c(membership_threshold, color_max),
            breaks = seq(membership_threshold, color_max, by = 0.1),
            labels = function(x) sprintf("%.1f", x)
          ) +
          scale_alpha_continuous(
            range = c(0.4, 0.9),
            limits = c(membership_threshold, color_max),
            guide = "none"
          ) +
          geom_line(data = center_data,
                    aes(x = Stage, y = Expression, group = 1),
                    color = "black", size = 5, alpha = 0.9) +
          geom_point(data = center_data,
                     aes(x = Stage, y = Expression),
                     color = "black", size = 6, shape = 19) +
          geom_point(data = center_data,
                     aes(x = Stage, y = Expression),
                     color = cluster_colors[k], size = 4.5, shape = 19) +
          scale_x_continuous(
            breaks = 1:ncol(expr_data),
            labels = levels(clinical_data$Clinical_Stage),
            expand = expansion(mult = c(0.05, 0.05))
          ) +
          coord_cartesian(xlim = c(0.8, ncol(expr_data) + 1.5)) +
          labs(
            title = paste("Merged Cluster", sprintf("(n = %d)", n_proteins)),
            x = "Clinical Stage",
            y = "Relative Abundance"
          )
        # ===== Right-side labels: Top 10 proteins =====
        if (nrow(individual_data) > 0) {
          protein_info <- individual_data %>%
            dplyr::filter(Stage == ncol(expr_data)) %>%
            dplyr::arrange(desc(Membership)) %>%
            dplyr::select(Protein, Expression, Membership) %>%
            dplyr::distinct()
          
          max_names <- min(15, nrow(protein_info))   # ✅ Top10
          protein_info <- protein_info[1:max_names, ]
          
          # ✅ 必须定义 y_range
          y_range <- range(c(individual_data$Expression, center_data$Expression), na.rm = TRUE)
          
          # 给标题+名字留间距（都在图内）
          top_pad <- 0.12 * diff(y_range)
          name_top <- y_range[2] - top_pad
          name_bottom <- y_range[1] + 0.06 * diff(y_range)
          
          y_positions <- seq(from = name_top, to = name_bottom, length.out = max_names)
          
          # ✅ 标题放在图内，不会被裁剪
          p <- p + annotate(
            "text",
            x = ncol(expr_data) + 0.2,
            y = y_range[2] - 0.02 * diff(y_range),
            label = "Top 15",
            hjust = 0, vjust = 1,
            size = 8,
            color = "grey20",
            fontface = "bold"
          )
          
          # 右侧蛋白名
          for (i in 1:nrow(protein_info)) {
            text_color <- grDevices::colorRamp(c(alpha(cluster_colors[k], 0.5), cluster_colors[k]))(
              (protein_info$Membership[i] - membership_threshold) / (color_max - membership_threshold)
            )
            text_color_hex <- grDevices::rgb(text_color[1], text_color[2], text_color[3], maxColorValue = 255)
            
            p <- p + annotate(
              "text",
              x = ncol(expr_data) + 0.2,
              y = y_positions[i],
              label = protein_info$Protein[i],
              hjust = 0, vjust = 0.5,
              size = 7.2,
              color = text_color_hex,
              fontface = ifelse(
                protein_info$Membership[i] > (membership_threshold + color_max) / 2,
                "bold", "plain"
              )
            )
          }
          
          # 可选：更多提示
          n_proteins <- length(unique(individual_data$Protein))
          if (n_proteins > max_names) {
            p <- p + annotate(
              "text",
              x = ncol(expr_data) + 0.2,
              y = y_range[1] + 0.02 * diff(y_range),
              label = paste("... and", n_proteins - max_names, "more"),
              hjust = 0, vjust = 0,
              size = 6.5,
              color = "grey50",
              fontface = "italic"
            )
          }
        }
        
        
        
        cluster_plots[[k]] <- p
      }
      
      cluster_plots <- cluster_plots[!sapply(cluster_plots, is.null)]
      
      # 保存为 PNG（推荐）
      png("Figure4F.png",
          width = 8, height = 8, units = "in", res = 300)
      
      combined_plot <- grid.arrange(grobs = cluster_plots, ncol = n_clusters)
      print(combined_plot)
      
      dev.off()
      
    }
    
    # 调用合并版 ggplot 可视化
    overall_max_membership_merged <- 1  # 合并后 membership = 1
    color_max_merged <- overall_max_membership_merged
    create_enhanced_ggplot_version_merged(eset_std, cl_plot,
                                          membership_threshold = 0.2,
                                          color_max = color_max_merged)
    cat(">>> 已生成合并后的 ggplot 轨迹图: ./mFuzz/enhanced_ggplot_proteomics_clustering_differential_proteins_merged.pdf\n")
    
    
    # -----------------------------
    # 14.3 参考图风格的单一合并 Cluster 折线图
    # -----------------------------
    dir.create("proteomics_cluster_plots_benign_baseline_merged", showWarnings = FALSE)
    
    stage_labels_simple <- levels(clinical_data$Clinical_Stage)
    
    p_merged_simple <- plot_cluster_proteins_simple(
      cluster_num = 1,
      mfuzz_result = cl_plot,
      standardized_data = exprs(eset_std),
      min_membership = 0.2,
      max_proteins = 20,
      stage_labels = stage_labels_simple,
      show_center = FALSE,
      plot_width = 10,
      plot_height = 8,
      base_font_size = 14
    )
    
    if (!is.null(p_merged_simple)) {
      ggsave("proteomics_cluster_plots_benign_baseline_merged/merged_cluster_proteins_benign.pdf",
             p_merged_simple, width = 12, height = 10)
      ggsave("proteomics_cluster_plots_benign_baseline_merged/merged_cluster_proteins_benign.png",
             p_merged_simple, width = 12, height = 10, dpi = 300)
      cat(">>> 已生成合并后的参考风格折线图: proteomics_cluster_plots_benign_baseline_merged/merged_cluster_proteins_benign.*\n")
    }
    
  } else {
    cat(">>> 两个聚类的中心轨迹相关性不足（≤ 阈值），保留原始两个聚类，不做合并可视化。\n")
  }
  
} else {
  cat("当前聚类数 != 2，不执行“两个聚类合并展示”的逻辑。\n")
}
