# Load necessary packages

rm(list = ls())
library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(gghalves)
library(ggsignif)
library(ggbeeswarm)

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

setwd("/Volumes/SeagateHub/NIC_imaging_data/Tao/20260628/JoCAP/JoCAP_R/JoCAP_Otsu_C2C3_V2_figures/forfigureV2")

df <- read_csv(
  "/Volumes/SeagateHub/NIC_imaging_data/Tao/20260628/JoCAP/JoCAP_R/JoCAP_Otsu_C2C3_V2_figures/forfigureV2/consolidate_colocalisation_intensity_by_ROI_Ch2_1800_to_5500_filtered_corrected.csv",
  show_col_types = FALSE
)


plot_metric <- "Pearson's Coefficient"
#plot_metric <- "Ch2_Mean"
#plot_metric <- "Ch3_Mean"

# If Group does not exist or is empty, extract it from SampleID
df <- df %>%
  mutate(
    Group = ifelse(
      is.na(Group) | Group == "",
      str_extract(SampleID, "^[A-Z]\\d+_"),
      Group
    )
  )

# ------------------------------------------------------------
# Define custom group order
# ------------------------------------------------------------

custom_group_order <- c("B2_", "A2_")

# Filter and factor
df <- df %>%
  filter(Group %in% custom_group_order) %>%
  mutate(Group = factor(Group, levels = custom_group_order))

# ------------------------------------------------------------
# Choose metric to plot
# ------------------------------------------------------------

df_plot <- df %>%
  filter(!is.na(.data[[plot_metric]]))

# ------------------------------------------------------------
# Calculate sample sizes for x-axis labels
# ------------------------------------------------------------

sample_sizes <- df_plot %>%
  group_by(Group) %>%
  tally() %>%
  ungroup() %>%
  mutate(label = paste0(Group, "\n(n = ", n, ")")) %>%
  select(Group, label) %>%
  tibble::deframe()

# ------------------------------------------------------------
# Define custom color palette
# ------------------------------------------------------------

my_color <- colorRampPalette(c(
  "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
))(length(custom_group_order))

names(my_color) <- custom_group_order

# ------------------------------------------------------------
# Create the plot
# ------------------------------------------------------------

p <- ggplot(df_plot, aes(x = Group, y = .data[[plot_metric]], fill = Group)) +
  geom_violin(alpha = 1, trim = FALSE, width = 0.8, color = NA) +
  geom_quasirandom(
    shape = 21,
    size = 2.5,
    stroke = 0.5,
    width = 0.2,
    color = "black"
  ) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    alpha = 1,
    fill = "white",
    color = "black",
    size = 0.8
  ) +
  geom_signif(
    comparisons = list(
      c("B2_", "A2_")
    ),
    test = "t.test",
    map_signif_level = function(p) sprintf("p = %.2g", p),
    y_position = 0.95,
    step_increase = 0.15,
    tip_length = 0.02,
    size = 0.8,
    textsize = 5
  ) +
  theme_classic(base_size = 16) +
  labs(
    title = "DNAJC13mutants",
    x = NULL,
    y = plot_metric
  ) +
  scale_fill_manual(values = my_color) +
  scale_x_discrete(labels = sample_sizes) +
  theme(
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold.italic", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  ylim(-0.2, 1.2)

# Show plot
print(p)

# Save outputs
ggsave("JoCAP_with_n_C2C3_merged_by_ROI.png", plot = p, width = 3, height = 5, dpi = 300)
ggsave("JoCAP_with_n_C2C3_merged_by_ROI.eps", plot = p, width = 3, height = 5, device = "eps")
ggsave("JoCAP_with_n_C2C3_merged_by_ROI.pdf", plot = p, width = 3, height = 5)

