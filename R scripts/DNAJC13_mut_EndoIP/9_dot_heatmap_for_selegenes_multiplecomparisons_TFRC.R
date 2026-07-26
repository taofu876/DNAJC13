rm(list = ls())
options(stringsAsFactors = FALSE)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

data_file <- "all_stats_results_sa02558_02559_sa02561_02562_combine_annotate_organelles.csv"

genes_to_plot <- c(
  "TFRC", "TF", "M6PR", "LRP1", "LRP11", "LRP8", "IGF2R", "PROCR", "SORT1", "VLDLR"
)

comparison_list <- list(
  "157resWT.157KO"     = c(fc = "log2FC_157resWT.157KO",     p = "q.val_157resWT.157KO"),
  "157resmut2.157KO"   = c(fc = "log2FC_157resmut2.157KO",   p = "q.val_157resmut2.157KO"),
  "157resmut1.157KO"   = c(fc = "log2FC_157resmut1.157KO",   p = "q.val_157resmut1.157KO"),
  "157resmut1V2.157KO" = c(fc = "log2FC_157resmut1V2.157KO", p = "q.val_157resmut1V2.157KO")
)

df <- read_csv(data_file)

df_dot_plot <- bind_rows(lapply(names(comparison_list), function(group_name) {
  fc_col <- comparison_list[[group_name]]["fc"]
  p_col  <- comparison_list[[group_name]]["p"]
  
  df %>%
    filter(Gene.Symbol %in% genes_to_plot) %>%
    transmute(
      Gene.Symbol = Gene.Symbol,
      Group = group_name,
      log2FC = as.numeric(.data[[fc_col]]),
      q.val = as.numeric(.data[[p_col]])
    )
})) %>%
  mutate(
    q.val = ifelse(q.val <= 0, NA, q.val),
    `Significance (-log10 q)` = -log10(q.val),
    Gene.Symbol = factor(Gene.Symbol, levels = rev(genes_to_plot)),
    Group = factor(Group, levels = names(comparison_list))
  )

# Autoscale log2FC color range symmetrically around 0
fc_limit <- max(abs(df_dot_plot$log2FC), na.rm = TRUE)

# Autoscale -log10(q) dot size breaks
sig_min <- floor(min(df_dot_plot$`Significance (-log10 q)`, na.rm = TRUE))
sig_max <- ceiling(max(df_dot_plot$`Significance (-log10 q)`, na.rm = TRUE))
sig_breaks <- pretty(c(sig_min, sig_max), n = 4)

p <- ggplot(df_dot_plot, aes(x = Group, y = Gene.Symbol)) +
  geom_point(
    aes(size = `Significance (-log10 q)`, fill = log2FC),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  scale_fill_gradient2(
    low = "#04a3ff",
    mid = "white",
    high = "#ff349c",
    midpoint = 0,
    name = expression(log[2]~FC),
    limits = c(-fc_limit, fc_limit),
    breaks = pretty(c(-fc_limit, fc_limit), n = 5)
  ) +
  scale_size_continuous(
    range = c(2, 10),
    limits = c(sig_min, sig_max),
    breaks = sig_breaks,
    name = expression(-log[10](italic(q) - value))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "lightgrey", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Comparison Group",
    y = "Gene",
    title = "Differential Protein Abundance"
  )

print(p)

ggsave(
  "9_dot_heatmap_for_selegenes_multiplecomparisons_TFRC.pdf",
  p,
  width = 5,
  height = 6
)