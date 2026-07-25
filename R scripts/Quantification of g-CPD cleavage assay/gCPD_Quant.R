# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(Cairo)
# --- 1. Data Processing Function (Robust for CSV header issues) ---
process_wb_data <- function(file_path) {
  df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  
  # Clean column names to prevent the "NA or "" names" error
  colnames(df) <- make.names(colnames(df), unique = TRUE)
  
  # Flexible column selection
  df <- df %>% 
    select(
      Group_Raw = matches("Group|Group_Raw"),
      Label = matches("Label"),
      Signal = matches("Signal")
    )
  
  # Process and Normalize
  df_processed <- df %>%
    mutate(Group = case_when(
      Group_Raw == "1486WT" ~ "WT",
      Group_Raw == "IP6" ~ "IP6",
      Group_Raw == "1486IP6" ~ "IP6mut",
      Group_Raw == "1486IP6+IP6" ~ "IP6mut+IP6",
      TRUE ~ as.character(Group_Raw)
    )) %>%
    tidyr::extract(Label, into = c("Protein", "Time"), 
                   regex = "([a-zA-Z0-9]+)-(\\d+)", convert = TRUE) %>%
    group_by(Group, Protein, Time) %>%
    mutate(Replicate = row_number()) %>%
    ungroup() %>%
    group_by(Group, Protein, Replicate) %>%
    mutate(Fold_Change = Signal / Signal[Time == 0]) %>%
    ungroup()
  
  return(df_processed)
}

# --- 2. Load and Combine Data ---
# Ensure these filenames match your local files
data1 <- process_wb_data("gCPD_RoiSet.csv")
data2 <- process_wb_data("DNAJC13_1_1486_RoiSet.csv")

full_data <- bind_rows(data1, data2)

# Calculate Mean and SD
summary_stats <- full_data %>%
  group_by(Group, Protein, Time) %>%
  summarise(
    mean_fc = mean(Fold_Change, na.rm = TRUE),
    sd_fc = sd(Fold_Change, na.rm = TRUE),
    .groups = 'drop'
  )

# --- 3. Shaded Band Plotting & Saving Function ---
plot_shaded_bands <- function(data, protein_name) {
  
  # Set Group factor levels for legend order
  data$Group <- factor(data$Group, levels = c("WT", "IP6", "IP6mut", "IP6mut+IP6"))
  
  plot_data <- data %>% filter(Protein == protein_name)
  
  p <- ggplot(plot_data, aes(x = Time, y = mean_fc, color = Group, fill = Group)) +
    # 1. Add baseline reference
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", size = 0.5) +
    
    # 2. Add Shaded Bands (Ribbon)
    geom_ribbon(aes(ymin = mean_fc - sd_fc, ymax = mean_fc + sd_fc), 
                color = NA, alpha = 0.15) + 
    
    # 3. Add Mean Lines and Points
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    
    # 4. Professional Colors
    scale_color_manual(values = c("WT" = "#000000", "IP6" = "#2166ac", 
                                  "IP6mut" = "#b2182b", "IP6mut+IP6" = "#1b7837")) +
    scale_fill_manual(values = c("WT" = "#000000", "IP6" = "#2166ac", 
                                 "IP6mut" = "#b2182b", "IP6mut+IP6" = "#1b7837")) +
    
    # 5. Theme and Labels
    labs(title = paste("Kinetic Profile:", protein_name),
         subtitle = "Shaded area represents ± Standard Deviation (n=3)",
         x = "Time (min)", 
         y = "Relative Expression (Fold Change)") +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      axis.line = element_line(size = 0.8),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # --- SAVE BLOCK ---
  # Saves each plot in 3 formats at a standard 7x5 inch size
  ggsave(paste0(protein_name, "_plot.png"), plot = p, width = 7, height = 5, dpi = 300)
  ggsave(paste0(protein_name, "_plot.pdf"), plot = p, width = 7, height = 5)
  ggsave(paste0(protein_name, "_plot.eps"), plot = p, width = 7, height = 5)
  
  return(p)
}

# --- 4. Generate and Save Plots ---
# This will display the plots and save the files to your working directory
p_d13  <- plot_shaded_bands(summary_stats, "D13")
p_gcpd <- plot_shaded_bands(summary_stats, "gCPD")
p_cpd  <- plot_shaded_bands(summary_stats, "CPD")

# View them in RStudio
print(p_d13)
print(p_gcpd)
print(p_cpd)
