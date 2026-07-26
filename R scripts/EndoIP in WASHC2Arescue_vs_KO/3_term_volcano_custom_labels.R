rm(list = ls())
options(stringsAsFactors = FALSE)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(ggrepel)

# Update only this section ==========================================
setwd("/Users/taofu/HMS Dropbox/Fu Tao/tao/3_MassSpec/20250510_ALFA endoIP in FAM21A KO and AP1B1 KO cells with rescue repeat for 20250418/R_plot/")

file_path <- "./all_stats_results_sa02564_02565_SDS.csv"

# --- NEW: CUSTOM LABEL LIST ---
# Add the Gene Symbols you want to see on the plot here
genes_to_label <- c(
  "WASHC2A", "CAPZB", "CCDC93", "CCDC22", "CAPZA2", "VPS35L", "CAPZA1", 
  "FKBP15", "SEZ6L2", "ABCA2", "LRP1", "COMMD3", "SNX7", "ATP6V0D1", 
  "ATP6V0C", "SNX4",  "WASHC1",  "LEPROTL1", "PIKFYVE",  "LEPROT", 
  "WASHC3", "ECE1",  "ATP6V0A1"
)

# Define the x-axis label as a variable
x_axis_label <- TeX("$Log_2 \\ FC\\ (232res\\ vs.\\ 232KO)$")
log2fc_column <- "log2FC_232res_D13.232KO_D13"
pvalue_column <- "q.val_232res_D13.232KO_D13"
# ===================================================================

# 读取数据：
data <- read.csv(file_path)
term_data <- read.csv("3_color_by_groups_refined.csv")

if (nrow(data) == 0 | nrow(term_data) == 0) {
  stop("One of the input dataframes is empty.")
}

data <- na.omit(data)

# Flip the plot
#data[[log2fc_column]] <- data[[log2fc_column]] * -1 

# 添加GO一列：
data$GO_term <- "others"
term_data <- term_data[term_data$Gene.names %in% data$Gene.Symbol,]
data$GO_term[data$Gene.Symbol %in% term_data$Gene.names] <- term_data$term[match(data$Gene.Symbol[data$Gene.Symbol %in% term_data$Gene.names], term_data$Gene.names)]

target_order <- unique(term_data$term)
data$GO_term <- factor(data$GO_term, levels = c(target_order, "others"))

unique_terms <- target_order 
specified_colors <- c("#dd4124","#b565a7","#49c2c6", "#fbcbcc", "#eef0ac", "#b1daa7", "#d0d0a0", "#6b5b95", "#88b04b", "#d65076", "#92a8d1","#955251",  "#009b77")

colors <- rep(specified_colors, length.out = length(unique_terms))
names(colors) <- unique_terms

# 计算上调下调数目：
Down_num <- length(which(data[[pvalue_column]] < 0.01 & data[[log2fc_column]] < -log2(1.5)))
Up_num <- length(which(data[[pvalue_column]] < 0.01 & data[[log2fc_column]] > log2(1.5)))

# --- UPDATED LABEL LOGIC ---
# Only create a label if the gene is in your custom 'genes_to_label' list
data$label <- ifelse(data$Gene.Symbol %in% genes_to_label, data$Gene.Symbol, NA)

ggplot(data[which(data$GO_term != "others"), ],
       aes_string(x = log2fc_column, y = paste0("-log10(", pvalue_column, ")"), fill = "GO_term")) +
  geom_point(data = data[which(data$GO_term == "others"), ],
             aes_string(x = log2fc_column, y = paste0("-log10(", pvalue_column, ")")),
             size = 0.5, color = "#999999") +
  geom_point(size = 3, shape = 21, color = "black") +
  scale_fill_manual(values = colors) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "#999999") +
  labs(x = x_axis_label, y = TeX("$-Log_{10} \\ p-value$")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "#fefde2", colour = "black", size = 0.2),
        legend.key = element_rect(fill = "#fefde2"),
        legend.key.size = unit(15, "pt"),
        legend.title = element_blank(),
        title = element_text(size = 15), 
        text = element_text(size = 15)) +
  annotate("text", label = "bolditalic(Down)", parse = TRUE, x = -0.9, y = 21, size = 4) +
  annotate("text", label = "bolditalic(Up)", parse = TRUE, x = 0.9, y = 21, size = 4) +
  annotate("text", label = Down_num, x = -0.9, y = 20, size = 3) +
  annotate("text", label = Up_num, x = 0.9, y = 20, size = 3) +
  # Adding labels only for selected genes
  geom_text_repel(aes(label = label), 
                  size = 3, 
                  max.overlaps = Inf, # High value ensures your custom labels always show
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "black")

ggsave("3_term_volcano_custom_labels.pdf", height = 6, width = 6)

