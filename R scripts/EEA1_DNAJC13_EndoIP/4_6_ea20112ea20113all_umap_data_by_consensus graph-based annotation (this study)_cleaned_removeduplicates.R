rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(ggpmisc)

# Load the data
data1 <- read.csv("./ea20112ea20113all_umap_data_by_consensus graph-based annotation (this study)_cleaned_removeduplicates.csv", header = TRUE)

# 🔑 KEY MODIFICATION 1: Define the order for organelles and colors
#selected_organelles <- c("Cytosol", "Nuclear","nucleolus", "mitochondria","peroxisome","plasma_membrane","ER", "ERGIC","early_endosome","recycling_endosome","lysosome", "Golgi", "trans-Golgi","unclassified","peroxisome","actin-binding_protein")
selected_organelles <- c("Cytosol", "Nuclear","mitochondria","plasma membrane","ER", "ERGIC","early_endosome","recycling_endosome","lysosome", "Golgi", "trans-Golgi")
# 🔑 KEY MODIFICATION 1: Define the order for organelles and colors
#selected_organelles <- c("cytoplasm", "nucleus", "mitochondrion","peroxisome","endosome","early.endosome", "late.endosome","Recycling.endosome","lysosome",  "endoplasmic.reticulum", "golgi", "single.pass.membrane.protein", "multi.pass.membrane.protein", "trans.Golgi.network","vesicle")

# Define the full color palette (Must have at least as many colors as selected_organelles)
my_color_full <- c('#171C1D','#092B42','#3f6f76','#BEC4DB','#7632AC','#1B9553','#61c1d7','#085990','#ea8d47','#a2d94d','#c94c4c')

#my_color_full <- c('#092B42', '#105186','#FDB338','#171C1D','#E86349','#9E8774','#425D5F','#7632AC','#BEC4DB','#1B9553','#E35D24','#E21E2E','#5EC4BE','#3f6f76','#a2d94d')



# Filter the data for selected organelles
data1_filtered <- data1[data1$Organelle %in% selected_organelles, ]

# --- NEW CODE: Create and save the list of filtered data frames ---

# Select only the required columns from the filtered data
selected_cols_data <- data1_filtered[c("Gene.Symbol", 
                                       "log2FC_EEA1_rescue.EEA1KO", 
                                       "q.val_EEA1_rescue.EEA1KO", 
                                       "log2FC_DNAJC13rescue.DNAJC13KO", 
                                       "q.val_DNAJC13rescue.DNAJC13KO", 
                                       "Organelle")]

# Use the split function to create a list where each element is a data frame 
# for a single organelle. The names of the list elements will be the organelle names.
organelle_data_list <- split(selected_cols_data, f = selected_cols_data$Organelle)

# Save the list to an RDS file
saveRDS(organelle_data_list, "4_6_ea20112ea20113all_umap_data_by_consensus graph-based annotation (this study)_cleaned_removeduplicates.rds")

# Create a new directory to store the CSVs to keep files organized
dir.create("4_6_ea20112ea20113all_umap_data_by_consensus graph-based annotation (this study)_cleaned_removeduplicates", showWarnings = FALSE)

# Loop through the list and save each data frame as a CSV
for (organelle_name in names(organelle_data_list)) {
  # Construct the file name
  filename <- paste0("./4_6_ea20112ea20113all_umap_data_by_consensus graph-based annotation (this study)_cleaned_removeduplicates/", organelle_name, "_stats.csv")
  
  # Write the data frame to a CSV file (row.names = FALSE prevents an extra column)
  write.csv(organelle_data_list[[organelle_name]], file = filename, row.names = FALSE)
  
  # Optional: Print confirmation
  cat("Saved:", filename, "\n")
}
# -------------------------------------------------------------------------

# Calculate mean and standard deviation for the groups (This part remains for the plot)
log2FC_EG <- aggregate(log2FC_EEA1_rescue.EEA1KO ~ Organelle, data1_filtered, function(x) c(mean = mean(x), sd = sd(x)))
log2FC_C1 <- aggregate(log2FC_DNAJC13rescue.DNAJC13KO ~ Organelle, data1_filtered, function(x) c(mean = mean(x), sd = sd(x)))

# Combine means and SDs
data2 <- data.frame(
  Organelle = log2FC_EG$Organelle,
  log2FC_EG_mean = log2FC_EG$log2FC_EEA1_rescue.EEA1KO[, 1],
  log2FC_EG_sd = log2FC_EG$log2FC_EEA1_rescue.EEA1KO[, 2],
  log2FC_C1_mean = log2FC_C1$log2FC_DNAJC13rescue.DNAJC13KO[, 1],
  log2FC_C1_sd = log2FC_C1$log2FC_DNAJC13rescue.DNAJC13KO[, 2]
)

# --- Gene Count Calculation & Labeling ---
gene_counts <- aggregate(data1_filtered$Organelle, by = list(Organelle = data1_filtered$Organelle), FUN = length)
names(gene_counts)[2] <- "Gene_Count"
data2 <- merge(data2, gene_counts, by = "Organelle")

data2$Gene_Count_Char <- as.character(data2$Gene_Count)
# Create the labeled factor for the legend (e.g., cytosol (465))
data2$Organelle_Label <- factor(paste0(data2$Organelle, " (", data2$Gene_Count_Char, ")"))

# 🔑 KEY MODIFICATION 2: Set the levels of the new label factor based on selected_organelles
# 1. Create the target order of the *Labels* based on the order of *Organelles*
target_label_order <- data2$Organelle_Label[match(selected_organelles, data2$Organelle)]
# 2. Set the factor levels to this explicit order
data2$Organelle_Label <- factor(data2$Organelle_Label, levels = target_label_order)

# -------------------------------------------------------------------

# 🔑 KEY MODIFICATION 3: Map the color palette to the new label order
# 1. Get the list of the final *labeled* organelles present in data2, in the desired order
final_labels <- levels(data2$Organelle_Label)
# 2. Trim the color palette to match the number of *present* organelles
num_clades <- length(final_labels)
assigned_colors <- my_color_full[1:num_clades] # Use only the first N colors
# 3. Name the color vector using the final_labels (Ensures colors map correctly)
names(assigned_colors) <- final_labels

# Create the plot using the FULL data frame: data2
p_all <- ggplot(data2, aes(log2FC_EG_mean, log2FC_C1_mean)) +
  # Use Organelle_Label for color, Gene_Count for size
  geom_errorbar(aes(xmin = log2FC_EG_mean - log2FC_EG_sd, xmax = log2FC_EG_mean + log2FC_EG_sd, color = Organelle_Label), linewidth = 0.8, width = 0) +  # X error bars
  geom_errorbar(aes(ymin = log2FC_C1_mean - log2FC_C1_sd, ymax = log2FC_C1_mean + log2FC_C1_sd, color = Organelle_Label), linewidth = 0.8, width = 0) +  # Y error bars
  geom_point(aes(color = Organelle_Label, size = Gene_Count)) +  # Points with size and color based on label
  
  # Scales: Pass the named color vector to scale_color_manual
  scale_color_manual(values = assigned_colors) +
  scale_size_continuous(name = "Number of Genes", range = c(2, 6)) +  
  
  # Regression and Stats (Calculated on ALL points in data2)
  geom_smooth(method = "lm", color = "black", fill = "#91D1C2FF", se = TRUE,  
              formula = y ~ x,
              linetype = 1, alpha = 0.2) +
  stat_poly_eq(formula = y ~ x,  
               aes(label = paste(eq.label, after_stat(adj.rr.label), ..p.value.label.., sep = "~~~~")), parse = TRUE,
               size = 4.5, label.x = "left") +
  
  # Labels and Theme
  labs(x = "log2FC_EG.293T", y = "log2FC_C1_aD13mSG.Clone1_KO", color = NULL) +
  theme_bw() +
  
  # Add horizontal and vertical lines at 0, 0 for reference
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  
  # Legend adjustments
  theme(legend.box = "horizontal",
        legend.position = "right")

print(p_all)

# Save the plot
ggsave("4_6_ea20112ea20113all_umap_data_by_consensus graph-based annotation (this study)_cleaned_removeduplicates.pdf", plot = p_all, height = 6, width = 9)
ggsave("4_6_ea20112ea20113all_umap_data_by_consensus graph-based annotation (this study)_cleaned_removeduplicates.eps", plot = p_all, height = 6, width = 9)
ggsave("4_6_ea20112ea20113all_umap_data_by_consensus graph-based annotation (this study)_cleaned_removeduplicates.png", plot = p_all, height = 6, width = 9)

