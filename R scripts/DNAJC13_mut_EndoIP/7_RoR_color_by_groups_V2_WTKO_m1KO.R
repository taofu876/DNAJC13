rm(list = ls())
library(MASS)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(readr)
library(broom)
library(dplyr)
library(RColorBrewer) 

#########################################################################
# 1. 配置部分 (Configuration)
#########################################################################
# Set your working directory
setwd("/Users/taofu/HMS Dropbox/Fu Tao/tao/3_MassSpec/20250506_ALFAendoIP_D13mutants/R_plot/")

main_data_file  <- "all_stats_results_sa02558_02559_sa02561_02562_combine_annotate_organelles.csv"
color_data_file <- "7_color_by_groups_refined.csv" 

col_x     <- "log2FC_157resWT.157KO"
col_y     <- "log2FC_157resmut1.157KO"
col_gene  <- "Gene.Symbol"

#########################################################################
# 2. 数据处理与颜色逻辑 (Data Processing)
#########################################################################

data_raw <- read.csv(main_data_file)
term_data <- read.csv(color_data_file)

# Join the main data with your classification file
data <- data.frame(
  x = data_raw[[col_x]], 
  y = data_raw[[col_y]], 
  Gene.Symbol = data_raw[[col_gene]]
) %>%
  left_join(term_data, by = c("Gene.Symbol" = "Gene.names"))

# Handle missing terms
data$term[is.na(data$term)] <- "Other"

# --- FIX 1: Explicitly define factor levels so your groups come first ---
# This ensures Group 1 gets Color 1 (Red), Group 2 gets Color 2 (Blue), etc.
unique_terms <- setdiff(unique(data$term), "Other")
data$term <- factor(data$term, levels = c(unique_terms, "Other"))

# Dynamic Color Palette Generation
all_levels <- levels(data$term)
num_groups <- length(all_levels)

# --- FIX 2: Corrected base_colors (removed trailing comma) ---
base_colors <- c("#c0392b", "#085990", "#8e44ad", "#27ae60", "#e67e22")

if(num_groups > length(base_colors)) {
  # If you have more groups than colors, generate more
  extra_colors <- colorRampPalette(brewer.pal(8, "Set3"))(num_groups - length(base_colors))
  final_palette <- c(base_colors, extra_colors)
} else {
  final_palette <- base_colors[1:num_groups]
}

# Assign names to colors to map them correctly to terms
names(final_palette) <- all_levels

# --- FIX 3: Ensure "Other" is grey ---
if("Other" %in% names(final_palette)) {
  final_palette["Other"] <- "grey75" # Light grey for background
}

# Sort data so colored points (your groups) are plotted last (on top)
data <- data %>% arrange(desc(term == "Other"))

# Identify labels (only for non-"Other" genes)
label_data <- data %>% 
  filter(term != "Other") %>% 
  distinct(Gene.Symbol, .keep_all = TRUE)

#########################################################################
# 3. 统计信息 (Statistics)
#########################################################################

lm_model <- lm(y ~ x, data = data)
lm_summary <- summary(lm_model)
r_squared  <- lm_summary$r.squared
p_value    <- tidy(lm_model)$p.value[2]
intercept  <- coef(lm_model)[1]
slope      <- coef(lm_model)[2]

# Axis helper function
draw_axis_line_with_edges <- function(length_x, length_y, tick_step = 1, lab_step = 2) {
  tick_x_frame <- data.frame(ticks = seq(-length_x, length_x, by = tick_step))
  tick_y_frame <- data.frame(ticks = seq(-length_y, length_y, by = tick_step))
  lab_x_frame  <- subset(data.frame(lab = seq(-length_x, length_x, by = lab_step), zero = 0), lab != 0)
  lab_y_frame  <- subset(data.frame(lab = seq(-length_y, length_y, by = lab_step), zero = 0), lab != 0)
  
  ggplot() +
    geom_segment(aes(y = 0, yend = 0, x = -length_x, xend = length_x), linewidth = 0.5) +
    geom_segment(aes(x = 0, xend = 0, y = -length_y, yend = length_y), linewidth = 0.5) +
    geom_segment(data = tick_x_frame, aes(x = ticks, xend = ticks, y = 0, yend = -0.05)) +
    geom_segment(data = tick_y_frame, aes(x = 0, xend = -0.05, y = ticks, yend = ticks)) + 
    geom_text(data = lab_x_frame, aes(x = lab, y = -0.4, label = lab), size = 3.5) +
    geom_text(data = lab_y_frame, aes(x = -0.4, y = lab, label = lab), size = 3.5) +
    geom_rect(aes(xmin = -length_x, xmax = length_x, ymin = -length_y, ymax = length_y), fill = NA, color = "black") +
    theme_minimal() + 
    theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank())
}

equation_parsed <- paste0("italic(y) == ", round(slope, 2), "*italic(x) + ", round(intercept, 2))
stats_text <- paste0("R^2 == ", round(r_squared, 3), "~~italic(P) == ", format(p_value, digits = 3))

#########################################################################
# 4. 绘图 (Plotting)
#########################################################################

# Adjusting axis range to match your coord_cartesian better
p1 <- draw_axis_line_with_edges(7, 7) + 
  geom_point(data = data, aes(x, y, color = term), size = 2.2, alpha = 1) +
  
  geom_text_repel(
    data = label_data, 
    aes(x, y, label = Gene.Symbol),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.4,
    segment.color = "black",
    force = 12,
    min.segment.length = 0
  ) +
  
  geom_smooth(data = data, aes(x = x, y = y), 
              method = "lm", se = FALSE, linetype = "dashed", color = "black", linewidth = 0.7) +
  
  scale_color_manual(values = final_palette) + 
  labs(x = col_x, y = col_y, color = "Groups") +
  coord_cartesian(xlim = c(-2.5, 6.5), ylim = c(-2.5, 6.5)) +
  theme(legend.position = "right", 
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10)) +
  
  annotate("text", x = 6.2, y = -1.5, label = equation_parsed, parse = TRUE, hjust = 1, size = 4.5) +
  annotate("text", x = 6.2, y = -2.1, label = stats_text, parse = TRUE, hjust = 1, size = 4.5)

print(p1)

# Save files
ggsave("7_RoR_color_by_groups_V2_WTKO_m1KO_V2.png", plot = p1, height = 8, width = 12, dpi = 300)
ggsave("7_RoR_color_by_groups_V2_WTKO_m1KO_V2.pdf", plot = p1, height = 8, width = 11, device = "pdf", useDingbats = FALSE)

