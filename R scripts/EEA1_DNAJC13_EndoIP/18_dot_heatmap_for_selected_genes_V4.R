# 📦 Load Necessary Libraries
rm(list = ls())

library(readr)     # For efficient data loading
library(dplyr)     # For data manipulation
library(tidyr)     # For data pivoting
library(ggplot2)   # For plotting
library(pheatmap)  # Keep, in case you use it elsewhere

# 📌 File and Directory Setup
# NOTE: Using 'here' is a best practice to avoid setwd() issues.
# You might need to install.packages("here")
# If you don't use 'here', manually ensure the file path is correct.
# library(here)
# data_file <- here("R_analysis", "aDNAJC13IP_aEEA1IP_sameplex", 
#                   "ea20112ea20113all_stats_results_exp1analysis1.csv")
setwd("~/HMS Dropbox/Fu Tao/tao/3_MassSpec/20241021EEA1_DNAJC13_ALFAendoIP/R_analysis/aDNAJC13IP_aEEA1IP_sameplex/1_final_figures")

# Temporary path based on your original setwd()
data_file <- "ea20112ea20113all_stats_results_exp1analysis1.csv" 

# 💾 Load and Prepare Data
df <- read_csv(data_file)

# 🧹 Data Cleaning and Row Naming
# Set Gene.Symbol as row names, making them unique first
df <- df %>%
  mutate(Gene.Symbol = make.names(Gene.Symbol, unique = TRUE)) %>%
  tibble::column_to_rownames("Gene.Symbol")

# 📋 Define Columns and Genes of Interest
log2fc_col_1 <- "log2FC_EEA1_rescue-EEA1KO"
pval_col_1 <- "q.val_EEA1_rescue-EEA1KO" 
log2fc_col_2 <- "log2FC_DNAJC13rescue-DNAJC13KO"
pval_col_2 <- "q.val_DNAJC13rescue-DNAJC13KO" 

# List of target genes (using the original order for plotting)
genes_to_plot <- c(
  "DNAJC13","EEA1","WASH2P","WASHC2A","WASHC2C","WASHC4","WASHC5","CAPZA1","CAPZB","TBC1D23","RAB11B","RAB11FIP1","RAB11FIP2","RAB6A","BICD2","RAB5A","RAB5B","RAB5C","RAB7A","ATP6V1A","ATP6V1B2","ATP6V1C1","ATP6V1D","ATP6V1E1","ATP6V1G1","ATP6V1H","ATP6AP1","ATP6AP2","ATP6V0A1","ATP6V0A2","ATP6V0C","ATP6V0D1","RUSC1","RUSC2","AP4E1","WDR11","C17orf75","FAM91A1")

# 🛠️ Subset and Pivot Data for Dot Heatmap
df_dot_plot <- df %>%
  # 1. Filter to include only the genes of interest
  filter(rownames(.) %in% genes_to_plot) %>%
  # 2. Select the necessary log2FC and p-value columns
  select(all_of(c(log2fc_col_1, pval_col_1, log2fc_col_2, pval_col_2))) %>%
  # 3. Convert rownames to a proper column
  tibble::rownames_to_column("Gene") %>%
  # 4. Pivot to long format for ggplot2
  pivot_longer(
    cols = -Gene,
    names_to = c(".value", "Group"),
    # The new columns created are named 'log2FC' and 'p.val'
    names_pattern = "(log2FC|q.val)_(.*)" 
  ) %>%
  # 5. Calculate -log10(p-value) for size aesthetic (better scaling)
  mutate(
    # NOTE: The resulting column from the pivot is named 'p.val'
    `Significance (-log10 p)`= -log10(q.val) 
  )

# 6. Format Gene order for the Y-axis (vertical)
df_dot_plot$Gene <- factor(
  df_dot_plot$Gene, 
  levels = rev(genes_to_plot)
)

# 📊 Generate the Dot Heatmap
# log2FC = color (log-fold change)
# -log10(p.val) = size (statistical significance)
p <- ggplot(df_dot_plot, aes(x = Group, y = Gene)) +
  # Use geom_point with shape 21 (filled circle) for better control over border/fill
  geom_point(
    aes(size = `Significance (-log10 p)`, fill = log2FC),
    shape = 21, 
    color = "black", # Border color of the dots
    stroke = 0.5     # Thickness of the border
  ) +
  
  # Color Scale for log2FC: Diverging scale centered at 0
  scale_fill_gradient2(
    #low = "#2166AC", mid = "white", high = "#B2182B", # Blue-White-Red Palette
    low = "#04a3ff", mid = "white", high = "#ff349c", # Blue-White-Red Palette#04a3ff
    midpoint = 0, 
    name = expression(log[2]~FC), # Use LaTeX for log2FC label
    limit = max(abs(df_dot_plot$log2FC)) * c(-0.5, 1) # Ensure scale is symmetrical
  ) +
  
  # Size Scale for Significance
  scale_size_continuous(
    range = c(2, 10), 
    name = expression(-log[10](italic(p) - value))
  ) +
  
  # 🎨 Theming for Publication Quality
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "lightgrey", linetype = "solid",size=0.5), # Adds grid for alignment
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  
  labs(
    x = "Comparison Group", 
    y = "Gene",
    title = "Dot Heatmap of Differential Protein Abundance"
  )

print(p)
# 

# 💾 Save Output
ggsave("18_dot_heatmap_for_selected_genes_V4.pdf", p, width = 4.2, height = 14)
ggsave("18_dot_heatmap_for_selected_genes_V4.eps", p, width = 4.2, height = 14, device = "eps")

