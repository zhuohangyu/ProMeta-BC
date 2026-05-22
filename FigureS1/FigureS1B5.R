# 清理工作空间
rm(list = ls())
gc()


# Updated TNM Stage Distribution with Separated Legends
# ===================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(viridis)
library(gridExtra)
library(grid)
library(scales)

# Load and preprocess data (same as before)
data <- read_csv("Clinical.csv")

data_clean <- data %>%
  mutate(
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

malignant_data <- data_clean %>% filter(Label == "Malignant")

# ===================================================
# Updated TNM Stage Distribution with Separated Colors
# ===================================================

# Prepare TNM data with proper stage ordering
tnm_enhanced_data <- malignant_data %>%
  select(T_stage, N_stage, M_stage) %>%
  pivot_longer(cols = everything(), names_to = "TNM_Type", values_to = "Stage") %>%
  filter(!is.na(Stage)) %>%
  mutate(
    TNM_Type = case_when(
      TNM_Type == "T_stage" ~ "T Stage",
      TNM_Type == "N_stage" ~ "N Stage", 
      TNM_Type == "M_stage" ~ "M Stage"
    ),
    # Create ordered factors for proper sorting
    Stage_ordered = case_when(
      TNM_Type == "T Stage" ~ factor(Stage, levels = c("T0", "Tis", "T1", "T2", "T3", "T4")),
      TNM_Type == "N Stage" ~ factor(Stage, levels = c("N0", "N1", "N2", "N3")),
      TNM_Type == "M Stage" ~ factor(Stage, levels = c("M0", "M1")),
      TRUE ~ factor(Stage)
    )
  ) %>%
  count(TNM_Type, Stage_ordered) %>%
  group_by(TNM_Type) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup() %>%
  # Create a combined factor for coloring
  mutate(
    Stage_Type_Combined = paste0(TNM_Type, "_", Stage_ordered),
    Stage_Label = as.character(Stage_ordered)
  )

# Define separate color palettes for each TNM type with correct T stage order
t_stage_colors <- c(
  "T0" = "#fff3e0",     # Very light orange
  "Tis" = "#ffe0b2",    # Light orange  
  "T1" = "#ffcc80",     # Medium light orange
  "T2" = "#ffb74d",     # Medium orange
  "T3" = "#ffa726",     # Medium dark orange
  "T4" = "#ff9800"      # Dark orange
)

n_stage_colors <- c(
  "N0" = "#f3e5f5",     # Very light purple
  "N1" = "#ce93d8",     # Light purple
  "N2" = "#ba68c8",     # Medium purple
  "N3" = "#9c27b0"      # Dark purple
)

m_stage_colors <- c(
  "M0" = "#e8eaf6",     # Very light blue
  "M1" = "#3f51b5"      # Dark blue
)

# Combine all colors into one palette
all_tnm_colors <- c(t_stage_colors, n_stage_colors, m_stage_colors)

# Updated TNM data with correct T stage ordering (T0, Tis, T1, T2, T3, T4)
tnm_enhanced_data <- malignant_data %>%
  select(T_stage, N_stage, M_stage) %>%
  pivot_longer(cols = everything(), names_to = "TNM_Type", values_to = "Stage") %>%
  filter(!is.na(Stage)) %>%
  mutate(
    TNM_Type = case_when(
      TNM_Type == "T_stage" ~ "T Stage",
      TNM_Type == "N_stage" ~ "N Stage", 
      TNM_Type == "M_stage" ~ "M Stage"
    ),
    # Create ordered factors with Tis between T0 and T1
    Stage_ordered = case_when(
      TNM_Type == "T Stage" ~ factor(Stage, levels = c("T0", "Tis", "T1", "T2", "T3", "T4")),
      TNM_Type == "N Stage" ~ factor(Stage, levels = c("N0", "N1", "N2", "N3")),
      TNM_Type == "M Stage" ~ factor(Stage, levels = c("M0", "M1")),
      TRUE ~ factor(Stage)
    )
  ) %>%
  count(TNM_Type, Stage_ordered) %>%
  group_by(TNM_Type) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup() %>%
  # Create a combined factor for coloring
  mutate(
    Stage_Type_Combined = paste0(TNM_Type, "_", Stage_ordered),
    Stage_Label = as.character(Stage_ordered)
  )

# Create a version with simple square+label legends
create_tnm_with_simple_legends <- function() {
  # Create the main plot without legend
  main_plot <- tnm_enhanced_data %>%
    ggplot(aes(x = TNM_Type, y = percentage, fill = Stage_Label)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.9, width = 0.7) +
    geom_text(aes(label = ifelse(percentage > 3, 
                                 paste0(Stage_Label, " ", round(percentage, 1), "%"), "")), 
              position = position_stack(vjust = 0.5), 
              size = 3.5, color = "black", fontface = "bold") +
    scale_fill_manual(values = all_tnm_colors) +
    labs(
      title = "TNM Stage Distribution",
      x = "Stage Type",
      y = "Frequency (%)"
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1), 
                       limits = c(0, 105)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 12, face = "bold"),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black", size = 0.5)
    )
  
  # Create simple square+label legend function with perfect squares
  create_simple_legend <- function(stages, colors, title) {
    # Create data frame for simple legend
    legend_df <- data.frame(
      y = length(stages):1,
      stage = stages,
      color = colors,
      stringsAsFactors = FALSE
    )
    
    ggplot(legend_df, aes(x = 1, y = y)) +
      # Create perfect square tiles
      geom_tile(aes(fill = stage), width = 0.7, height = 0.7, 
                alpha = 0.9, color = "grey60", size = 0.3) +
      geom_text(aes(x = 1.9, label = stage), size = 4.2, fontface = "bold", 
                color = "black", hjust = 0) +
      scale_fill_manual(values = setNames(colors, stages)) +
      labs(title = title) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", hjust = 0, margin = margin(b = 10)),
        plot.margin = margin(10, 30, 10, 10),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      coord_fixed(ratio = 1) +  # This ensures perfect squares
      scale_x_continuous(limits = c(0.5, 4)) +
      scale_y_continuous(limits = c(0.5, length(stages) + 0.5))
  }
  
  # Create individual simple legends with correct T stage order
  t_legend <- create_simple_legend(
    stages = c("T0", "Tis", "T1", "T2", "T3", "T4"),  # Correct order: Tis between T0 and T1
    colors = t_stage_colors,
    title = "T Stage"
  )
  
  n_legend <- create_simple_legend(
    stages = c("N0", "N1", "N2", "N3"),
    colors = n_stage_colors,
    title = "N Stage"
  )
  
  m_legend <- create_simple_legend(
    stages = c("M0", "M1"),
    colors = m_stage_colors,
    title = "M Stage"
  )
  
  # Combine legends vertically with proper spacing
  combined_legends <- grid.arrange(
    t_legend, 
    n_legend, 
    m_legend, 
    ncol = 1, 
    heights = c(2, 1.3, 0.7)
  )
  
  # Combine main plot with legends
  final_plot <- grid.arrange(main_plot, combined_legends, ncol = 2, widths = c(3.5, 1))
  
  return(final_plot)
}

# Create an even cleaner version with better spacing
create_tnm_clean_simple <- function() {
  # Enhanced main plot with "Stage XX%" format
  main_plot <- tnm_enhanced_data %>%
    ggplot(aes(x = TNM_Type, y = percentage, fill = Stage_Label)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.9, width = 0.7) +
    geom_text(aes(label = ifelse(percentage > 4, 
                                 paste0(Stage_Label, " ", round(percentage, 1), "%"), 
                                 ifelse(percentage > 1.5, Stage_Label, ""))), 
              position = position_stack(vjust = 0.5), 
              size = 3.2, color = "black", fontface = "bold") +
    scale_fill_manual(values = all_tnm_colors) +
    labs(
      title = "TNM Stage Distribution",
      x = "Stage Type",
      y = "Frequency (%)"
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1), 
                       limits = c(0, 105)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 12, face = "bold"),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black", size = 0.5)
    )
  
  # Enhanced legend function with perfect squares
  create_clean_legend <- function(stages, colors, title) {
    legend_df <- data.frame(
      x = 1,
      y = length(stages):1,
      stage = stages,
      color = colors,
      stringsAsFactors = FALSE
    )
    
    ggplot(legend_df, aes(x = x, y = y)) +
      # Create perfect square tiles
      geom_tile(aes(fill = stage), width = 0.8, height = 0.8, 
                alpha = 0.9, color = "grey70", size = 0.4) +
      geom_text(aes(x = 2.1, label = stage), size = 4.5, fontface = "bold", 
                color = "black", hjust = 0, vjust = 0.5) +
      scale_fill_manual(values = setNames(colors, stages)) +
      labs(title = title) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 13, face = "bold", hjust = 0, 
                                  margin = margin(b = 8, l = 10)),
        plot.margin = margin(8, 25, 15, 5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      coord_fixed(ratio = 1) +  # This ensures perfect squares
      scale_x_continuous(limits = c(0.5, 4.5), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0.3, length(stages) + 0.7), expand = c(0, 0))
  }
  
  # Create clean legends with correct T stage order
  t_legend_clean <- create_clean_legend(
    stages = c("T0", "Tis", "T1", "T2", "T3", "T4"),  # Correct order: Tis between T0 and T1
    colors = t_stage_colors,
    title = "T Stage"
  )
  
  n_legend_clean <- create_clean_legend(
    stages = c("N0", "N1", "N2", "N3"),
    colors = n_stage_colors,
    title = "N Stage"
  )
  
  m_legend_clean <- create_clean_legend(
    stages = c("M0", "M1"),
    colors = m_stage_colors,
    title = "M Stage"
  )
  
  # Combine legends with optimized spacing
  combined_legends_clean <- grid.arrange(
    t_legend_clean,
    n_legend_clean, 
    m_legend_clean,
    ncol = 1,
    heights = c(1.8, 1.2, 0.7)
  )
  
  # Final combination with better proportions
  final_clean <- grid.arrange(main_plot, combined_legends_clean, ncol = 2, widths = c(4, 1.2))
  
  return(final_clean)
}

p_tnm_clean <- create_tnm_clean_simple()

ggsave("FigureS1B5.png", p_tnm_clean, width = 14, height = 8, dpi = 300)
