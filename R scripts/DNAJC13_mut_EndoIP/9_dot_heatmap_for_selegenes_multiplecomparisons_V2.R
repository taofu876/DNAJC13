# 📦 Load Necessary Libraries
rm(list = ls())
options(stringsAsFactors = FALSE)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# 📌 1. File & Gene Setup
# Set your working directory if needed:
# setwd("/Users/taofu/HMS Dropbox/Fu Tao/tao/3_MassSpec/20250506_ALFAendoIP_D13mutants/R_plot/")

data_file <- "all_stats_results_sa02558_02559_sa02561_02562_combine_annotate_organelles.csv" 

genes_to_plot <- c("DNAJC13", "WASHC1", "WASHC2A","WASHC2C","WASHC3","WASHC4","WASHC5", "CAPZA1","CAPZB","FKBP15","VPS35","VPS29","VPS26A","VPS26B","TBC1D5","RAB5A","RAB5B","RAB7A","RAB11A","RAB11B","RAB11FIP1","RAB11FIP5","ATP6V1A","ATP6V1B2","ATP6V1C1","ATP6V1D","ATP6V1E1","ATP6V1G1","ATP6V1H","ATP6AP1","ATP6AP2","ATP6V0A1","ATP6V0A2","ATP6V0C","ATP6V0D1")

# 📋 2. DEFINE YOUR COMPARISONS
comparison_list <- list(
  "157resWT.157KO"     = c(fc = "log2FC_157resWT.157KO",     p = "q.val_157resWT.157KO"),
  "157resmut2.157KO"   = c(fc = "log2FC_157resmut2.157KO",   p = "q.val_157resmut2.157KO"),
  "157resmut1.157KO"   = c(fc = "log2FC_157resmut1.157KO",   p = "q.val_157resmut1.157KO"),
  "157resmut1V2.157KO" = c(fc = "log2FC_157resmut1V2.157KO", p = "q.val_157resmut1V2.157KO")
)

# 💾 Load and Prepare Data
df <- read_csv(data_file) %>%
  mutate(Gene.Symbol = make.names(Gene.Symbol, unique = TRUE))

# 🛠️ Data Wrangling
df_dot_plot <- df %>%
  filter(Gene.Symbol %in% genes_to_plot) %>%
  select(Gene.Symbol, all_of(unname(unlist(comparison_list)))) %>%
  pivot_longer(-Gene.Symbol, names_to = "Full_Col", values_to = "Value") %>%
  mutate(
    Type = ifelse(grepl("log2FC", Full_Col), "log2FC", "q.val"),
    # FIXED: Extract the group name by removing the prefix 'log2FC_' or 'q.val_'
    Group = gsub("log2FC_|q.val_", "", Full_Col)
  ) %>%
  select(-Full_Col) %>%
  pivot_wider(names_from = Type, values_from = Value) %>%
  mutate(
    `Significance (-log10 p)` = -log10(q.val)
  )

# Set Factor Levels for Plotting
df_dot_plot$Gene.Symbol <- factor(df_dot_plot$Gene.Symbol, levels = rev(genes_to_plot))
df_dot_plot$Group <- factor(df_dot_plot$Group, levels = names(comparison_list))

# 📊 Generate the Dot Heatmap
p <- ggplot(df_dot_plot, aes(x = Group, y = Gene.Symbol)) +
  geom_point(aes(size = `Significance (-log10 p)`, fill = log2FC), 
             shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient2(
    low = "#04a3ff", mid = "white", high = "#ff349c", midpoint = 0, 
    name = expression(log[2]~FC),
    limit = max(abs(df_dot_plot$log2FC), na.rm = TRUE) * c(-1, 1)
  ) +
  scale_size_continuous(range = c(2, 10), 
                        name = expression(-log[10](italic(q) - value))) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "lightgrey", size = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(x = "Comparison Group", y = "Gene", title = "Differential Protein Abundance")

print(p)

# 💾 Save Output
ggsave("9_dot_heatmap_for_selegenes_multiplecomparisons_V2.pdf", p, width = 6, height = 12)








