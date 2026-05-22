# 清理工作空间
rm(list = ls())
gc()

# Load required libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(gridExtra)
library(grid)
library(reshape2)
library(scales)

# Custom Nature journal color palette
nature_colors <- c(
  primary = "#42a5f5",
  secondary = "#A23B72", 
  accent1 = "#F18F01",
  accent2 = "#DC3545",
  benign = "#28A745",
  malignant = "#DC3545"
)

# UNIFIED molecular colors (consistent with previous charts)
molecular_colors_unified <- c(
  "HR+/HER2-" = "#42a5f5",  # Light blue
  "HER2+" = "#DC3545",      # Pink
  "TNBC" = "#F18F01"        # Orange
)

# Additional bright color schemes for consistency
molecular_colors <- c(
  "HR+/HER2-" = "#42a5f5",  # Light blue
  "HER2+" = "#DC3545",      # Pink
  "TNBC" = "#F18F01"        # Orange
)

molecular_colors_gradient <- c(
  "HR+/HER2-" = "#81D4FA",  # Light blue gradient
  "HER2+" = "#F48FB1",      # Light pink gradient
  "TNBC" = "#FFCC80"        # Light orange gradient
)

# Set ggplot2 theme with LARGER FONTS
theme_nature <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "sans", size = 24),  # Increased from 11 to 16
      plot.title = element_text(size = 26, face = "bold", hjust = 0.5),  # Increased from 14 to 20
      axis.title = element_text(size = 24, face = "bold"),  # Increased from 12 to 18
      axis.text = element_text(size = 24),  # Increased from 10 to 16
      legend.title = element_text(size = 24, face = "bold"),  # Increased from 11 to 16
      legend.text = element_text(size = 22),  # Increased from 10 to 14
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      strip.text = element_text(size = 22, face = "bold"),  # Increased from 11 to 16
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
}

# ===================================================
# 1. Data Loading and Preprocessing
# ===================================================

# Load data
data <- read_csv("Clinical.csv")

# Data preprocessing
data_clean <- data %>%
  mutate(
    # Calculate BMI
    BMI = BMI,
    
    # BMI grouping
    BMI_Group = case_when(
      BMI < 24 ~ "Normal (<24)",
      BMI >= 24 & BMI < 28 ~ "Overweight (24-28)",
      BMI >= 28 ~ "Obese (>28)",
      TRUE ~ NA_character_
    ),
    
    # Age grouping
    Age_Group = case_when(
      Age < 40 ~ "<40 years",
      Age >= 40 & Age <= 60 ~ "40-60 years",
      Age > 60 ~ ">60 years",
      TRUE ~ NA_character_
    ),
    
    # Clinical stage grouping
    Stage_Group = case_when(
      `Clinical stage` == "0" ~ "Stage 0",
      str_detect(`Clinical stage`, "^I[AB]?$") ~ "Stage I",
      str_detect(`Clinical stage`, "^II[AB]?$") ~ "Stage II",
      str_detect(`Clinical stage`, "^III[ABC]?$") ~ "Stage III",
      str_detect(`Clinical stage`, "^IV") ~ "Stage IV",
      TRUE ~ NA_character_
    ),
    
    # Early vs Late stage
    Early_Late = case_when(
      Stage_Group %in% c("Stage 0", "Stage I", "Stage II") ~ "Early (0-II)",
      Stage_Group %in% c("Stage III", "Stage IV") ~ "Late (III-IV)",
      TRUE ~ NA_character_
    ),
    
    # Molecular classification
    Molecular_Type = case_when(
      `Molecular Classification` == 1 ~ "HR+/HER2-",
      `Molecular Classification` == 2 ~ "HER2+",
      `Molecular Classification` == 3 ~ "TNBC",
      TRUE ~ NA_character_
    ),
    
    # T stage cleaning
    T_stage = case_when(
      T == "is" ~ "Tis",
      T %in% c(0, 1, 2, 3, 4) ~ paste0("T", T),
      TRUE ~ NA_character_
    ),
    
    # N stage cleaning
    N_stage = case_when(
      N %in% c(0, 1, 2, 3) ~ paste0("N", N),
      TRUE ~ NA_character_
    ),
    
    # M stage cleaning
    M_stage = case_when(
      M %in% c(0, 1) ~ paste0("M", M),
      TRUE ~ NA_character_
    )
  )

# Separate benign and malignant data
benign_data <- data_clean %>% filter(Label == "Benign")
malignant_data <- data_clean %>% filter(Label == "Malignant")


# ===================================================
# Panel 1: Clinical Features Distribution in three cohorts
# ===================================================

# 1.1: Age distribution by cohort (additional analysis using real data)
age_by_cohort <- data_clean %>%
  filter(!is.na(Age)) %>%
  ggplot(aes(x = Cohort, y = Age, fill = Cohort, color = Cohort)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21, outlier.stroke = 1.2, 
               outlier.size = 3, width = 0.6, size = 1, outlier.alpha = 0.7) +  # Increased outlier size
  geom_jitter(width = 0.15, alpha = 0.6, size = 2, show.legend = FALSE) +  # Increased point size
  scale_fill_manual(values = c("Cohort_1" = "#2196F3", 
                               "Cohort_2" = "#FF9800", 
                               "Cohort_3" = "#4CAF50"),
                    name = "") +
  scale_color_manual(values = c("Cohort_1" = "#1565C0", 
                                "Cohort_2" = "#F57C00", 
                                "Cohort_3" = "#388E3C"),
                     name = "") +
  labs(
    title = "Age Distribution by Cohort",
    x = "",
    y = "Age (years)"
  ) +
  theme_nature() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 26, hjust = 0.5),  # Increased from 13 to 18
    axis.title = element_text(size = 24, face = "bold"),  # Increased from 11 to 16
    legend.title = element_text(size = 22, face = "bold"),  # Increased from 10 to 14
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave("FigureS1B1.png", age_by_cohort, width = 7, height = 8, dpi = 300)

# 1.2: BMI distribution by cohort (additional analysis using real data)
bmi_by_cohort <- data_clean %>%
  filter(!is.na(BMI)) %>%
  ggplot(aes(x = Cohort, y = BMI, fill = Cohort, color = Cohort)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21, outlier.stroke = 1.2, 
               outlier.size = 3, width = 0.6, size = 1, outlier.alpha = 0.7) +  # Increased outlier size
  geom_jitter(width = 0.15, alpha = 0.6, size = 2, show.legend = FALSE) +  # Increased point size
  scale_fill_manual(values = c("Cohort_1" = "#2196F3", 
                               "Cohort_2" = "#FF9800", 
                               "Cohort_3" = "#4CAF50"),
                    name = "") +
  scale_color_manual(values = c("Cohort_1" = "#1565C0", 
                                "Cohort_2" = "#F57C00", 
                                "Cohort_3" = "#388E3C"),
                     name = "") +
  labs(
    title = "BMI Distribution by Cohort",
    x = "",
    y = "BMI (kg/m²)"
  ) +
  theme_nature() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 26, hjust = 0.5),  # Increased from 13 to 18
    axis.title = element_text(size = 24, face = "bold"),  # Increased from 11 to 16
    legend.title = element_text(size = 22, face = "bold"),  # Increased from 10 to 14
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )


ggsave("FigureS1B2.png", bmi_by_cohort, width = 7, height = 8, dpi = 300)



# 1.3: Clinical stage distribution by cohort (using real data)
stage_by_cohort <- malignant_data %>%
  filter(!is.na(Stage_Group)) %>%
  count(Cohort, Stage_Group) %>%
  group_by(Cohort) %>%
  mutate(percentage = n / sum(n) * 100,
         total = sum(n)) %>%
  ungroup()

# Clinical stage colors (red gradient)
stage_colors <- c("Stage 0" = "#ffebee",
                  "Stage I" = "#ffcdd2", 
                  "Stage II" = "#ef5350", 
                  "Stage III" = "#f44336",
                  "Stage IV" = "#c62828")

p_stage_by_cohort <- stage_by_cohort %>%
  ggplot(aes(x = Cohort, y = percentage, fill = Stage_Group)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.9, width = 0.6) +
  geom_text(aes(label = ifelse(percentage > 5, paste0(round(percentage, 1), "%"), "")), 
            position = position_stack(vjust = 0.5), size = 5, color = "white", fontface = "bold") +
  scale_fill_manual(values = stage_colors, name = "") +
  labs(
    title = "Clinical Stages Distribution by Cohort",
    x = "",
    y = "Frequency"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_nature() +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 26, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 24, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_legend(nrow = 1, reverse = FALSE))  # 改为 reverse = FALSE


ggsave("FigureS1B3.png", p_stage_by_cohort, width = 9, height = 8, dpi = 300)


# 1.4: Molecular subtype distribution by cohort (using real data)
molecular_by_cohort <- malignant_data %>%
  filter(!is.na(Molecular_Type)) %>%
  count(Cohort, Molecular_Type) %>%
  group_by(Cohort) %>%
  mutate(percentage = n / sum(n) * 100,
         total = sum(n)) %>%
  ungroup()

# Molecular subtype colors (blue gradient)
molecular_colors_gradient <- c("HR+/HER2-" = "#90caf9",
                               "HER2+" = "#42a5f5", 
                               "TNBC" = "#1976d2")

p_molecular_by_cohort <- molecular_by_cohort %>%
  ggplot(aes(x = Cohort, y = percentage, fill = Molecular_Type)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.9, width = 0.6) +
  geom_text(aes(label = ifelse(percentage > 8, paste0(round(percentage, 1), "%"), "")), 
            position = position_stack(vjust = 0.5), size = 5, color = "white", fontface = "bold") +  # Increased from 3.5 to 5
  scale_fill_manual(values = molecular_colors_gradient, name = "") +
  labs(
    title = "Molecular Subtype Distribution by Cohort",
    x = "",
    y = "Frequency"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_nature() +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 26, hjust = 0.5, face = "bold"),  # Increased from 14 to 20
    axis.title = element_text(size = 24, face = "bold"),  # Increased from 12 to 18
    legend.title = element_text(size = 22, face = "bold"),  # Increased from 11 to 16
    panel.grid.major.x = element_blank(),
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_legend(nrow = 1, reverse = FALSE))

ggsave("FigureS1B4.png", p_molecular_by_cohort, width = 9, height = 8, dpi = 300)


# ===================================================
# Panel 2: Clinical Features Distribution in benign and BC
# ===================================================

# Enhanced Age Distribution Boxplot with vibrant colors
p_age_box_enhanced <- data_clean %>%
  filter(!is.na(Age)) %>%
  ggplot(aes(x = Label, y = Age, fill = Label, color = Label)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 1, outlier.stroke = 1, 
               outlier.size = 2, width = 0.5, size = 0.8) +  # Increased outlier size
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, show.legend = FALSE) +  # Increased point size
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     label.x = 1.5, label.y = max(data_clean$Age, na.rm = TRUE) * 0.92,
                     size = 8) +  # Increased statistical test label size
  scale_fill_manual(values = c("Benign" = "#28A745", 
                               "Malignant" = "#DC3545")) +
  scale_color_manual(values = c("Benign" = "#28A745", 
                                "Malignant" = "#DC3545")) +
  labs(
    title = "Age Distribution \n Between Benign and Malignant",
    x = "Pathological Type",
    y = "Age (years)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 20)),  # Increased from 14 to 20
    axis.title = element_text(size = 24, face = "bold"),  # Increased from 12 to 18
    axis.text = element_text(size = 22),  # Increased from 11 to 16
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Enhanced BMI Distribution Boxplot with vibrant colors
p_bmi_box_enhanced <- data_clean %>%
  filter(!is.na(BMI)) %>%
  ggplot(aes(x = Label, y = BMI, fill = Label, color = Label)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 1, outlier.stroke = 1, 
               outlier.size = 2, width = 0.5, size = 0.8) +  # Increased outlier size
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, show.legend = FALSE) +  # Increased point size
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     label.x = 1.5, label.y = max(data_clean$BMI, na.rm = TRUE) * 0.92,
                     size = 8) +  # Increased statistical test label size
  scale_fill_manual(values = c("Benign" = "#28A745", 
                               "Malignant" = "#DC3545")) +
  scale_color_manual(values = c("Benign" = "#28A745", 
                                "Malignant" = "#DC3545")) +
  labs(
    
    title = "BMI Distribution \n Between Benign and Malignant",
    x = "Pathological Type",
    y = "BMI (kg/m²)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5, margin = margin(b = 20)),  # Increased from 14 to 20
    axis.title = element_text(size = 24, face = "bold"),  # Increased from 12 to 18
    axis.text = element_text(size = 22),  # Increased from 11 to 16
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Also save as PNG for immediate viewing
ggsave("FigureS1C1.png", p_age_box_enhanced, width = 7, height = 8, dpi = 300)
ggsave("FigureS1C2.png", p_bmi_box_enhanced, width = 7, height = 8, dpi = 300)



# ===================================================
# Panel D2: Molecular Subtype and Clinical Feature Association
# ===================================================

# D1.2: Molecular subtype and early/late stage association heatmap
molecular_stage_data <- malignant_data %>%
  filter(!is.na(Molecular_Type) & !is.na(Early_Late)) %>%
  count(Molecular_Type, Early_Late) %>%
  group_by(Molecular_Type) %>%
  mutate(percentage = n / sum(n) * 100)

p_molecular_stage_heatmap <- molecular_stage_data %>%
  ggplot(aes(x = Early_Late, y = Molecular_Type, fill = percentage)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")), 
            color = "white", fontface = "bold", size = 5) +  # Increased text size
  scale_fill_gradient(low = nature_colors["primary"], high = nature_colors["accent2"],
                      name = "Percentage (%)") +
  labs(
    x = "Stage Group",
    y = "Molecular Subtype"
  ) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # Increased axis text size
        axis.text.y = element_text(size = 16))  # Increased axis text size
# Also save as PNG for immediate viewing
ggsave("FigureS1D1.png", p_molecular_stage_heatmap, width = 8, height = 8, dpi = 300)


# D3.1: Molecular subtype and T stage cross-tabulation visualization (with proper order)

# D3.1: Molecular subtype and T stage (updated colors)
molecular_t_data <- malignant_data %>%
  filter(!is.na(Molecular_Type) & !is.na(T_stage)) %>%
  count(Molecular_Type, T_stage) %>%
  group_by(T_stage) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(T_stage_ordered = factor(T_stage, 
                                  levels = c("T0", "Tis", "T1", "T2", "T3", "T4"),
                                  labels = c("T0", "Tis", "T1", "T2", "T3", "T4")))

p_molecular_t_updated <- molecular_t_data %>%
  ggplot(aes(x = T_stage_ordered, y = percentage, fill = Molecular_Type)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.9, width = 0.7) +
  geom_text(aes(label = ifelse(percentage > 5, paste0(n), "")), 
            position = position_stack(vjust = 0.5), size = 4.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = molecular_colors_unified, name = "Molecular Subtype") +
  labs(
    title = "Molecular Subtype and T Stage Association",
    x = "T Stage",
    y = "Percentage (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 21),
    legend.position = "bottom",
    legend.title = element_text(size = 21, face = "bold"),
    legend.text = element_text(size = 20),
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", size = 0.5)
  ) +
  guides(fill = guide_legend(nrow = 1))
ggsave("FigureS1D2.png", p_molecular_t_updated, width = 8, height = 8, dpi = 300)



# D3.2: Molecular subtype and N stage cross-tabulation visualization
molecular_n_data <- malignant_data %>%
  filter(!is.na(Molecular_Type) & !is.na(N_stage)) %>%
  count(Molecular_Type, N_stage) %>%
  group_by(N_stage) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(N_stage_ordered = factor(N_stage, 
                                  levels = c("N0", "N1", "N2", "N3"),
                                  labels = c("N0", "N1", "N2", "N3")))

p_molecular_n_updated <- molecular_n_data %>%
  ggplot(aes(x = N_stage_ordered, y = percentage, fill = Molecular_Type)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.9, width = 0.7) +
  geom_text(aes(label = ifelse(percentage > 5, paste0(n), "")), 
            position = position_stack(vjust = 0.5), size = 4.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = molecular_colors_unified, name = "Molecular Subtype") +
  labs(
    title = "Molecular Subtype and N Stage Association", 
    x = "N Stage",
    y = "Percentage (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 21),
    legend.position = "bottom",
    legend.title = element_text(size = 21, face = "bold"),
    legend.text = element_text(size = 20),
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", size = 0.5)
  ) +
  guides(fill = guide_legend(nrow = 1))
ggsave("FigureS1D3.png", p_molecular_n_updated, width = 8, height = 8, dpi = 300)

