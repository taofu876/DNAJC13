rm(list = ls())
library(tidyverse)
library(svglite)
library(ggsci)
library(tools)
setwd("/Volumes/SeagateHub/NIC_imaging_data/Tao/20260422/R_plot/V6_figure3")
# --- 1. 设置路径 ---
folder_path <- "/Volumes/SeagateHub/NIC_imaging_data/Tao/20260422/R_plot/V6_figure3"

# --- 2. 自动获取并读取所有 _filter.csv 文件 ---
file_list <- list.files(path = folder_path, pattern = "_filter\\.csv$", full.names = TRUE)

df_all <- file_list %>%
  map_df(function(x) {
    group_name <- sub(".*Group_(.*)_filter.*", "\\1", basename(x))
    read.csv(x) %>%
      mutate(Group = group_name)
  })

# --- 3. 计算各组样本量 (n=) ---
vesicle_counts <- df_all %>%
  group_by(Group) %>%
  summarise(n_vesicles = n_distinct(VesicleId), .groups = 'drop')

# --- 4. 统计处理：改用 Norm0 计算 Mean 和 SEM ---
plot_data <- df_all %>%
  group_by(FrameId, Group) %>%
  summarise(
    mean_val = mean(Norm0, na.rm = TRUE),  # <--- 已修改为 Norm0
    sem_val = sd(Norm0, na.rm = TRUE) / sqrt(n()), # <--- 已修改为 Norm0
    .groups = 'drop'
  ) %>%
  left_join(vesicle_counts, by = "Group") %>%
  mutate(
    LegendLabel = paste0(Group, " (n=", n_vesicles, ")")
  )

# --- 5. 定义配色 ---
unique_groups <- unique(plot_data$Group)
custom_colors <- setNames(pal_npg()(length(unique_groups)), unique_groups)

# --- 6. 绘图 ---
p <- ggplot(plot_data, aes(x = FrameId, y = mean_val, 
                           color = Group, 
                           fill = Group)) +
  geom_ribbon(aes(ymin = mean_val - sem_val, ymax = mean_val + sem_val), 
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.5) + 
  theme_classic(base_size = 12) + 
  labs(
    title = "Normalized Intensity Kinetics (Norm0)", # <--- 标题更新
    subtitle = "Mean ± s.e.m. (Normalized to Frame 0)",
    x = "Frame Number",
    y = "Normalized Intensity (F/F0)", # <--- Y轴标签更新
    color = NULL, 
    fill = NULL
  ) +
  scale_color_manual(values = custom_colors, 
                     labels = setNames(plot_data$LegendLabel, plot_data$Group)) +
  scale_fill_manual(values = custom_colors, 
                    labels = setNames(plot_data$LegendLabel, plot_data$Group)) +
  theme(
    legend.position = c(0.05, 0.95), 
    legend.justification = c("left", "top"),
    legend.background = element_blank(),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 18),
    axis.title = element_text(size = 16, face = "bold")
  )

print(p)

# --- 7. 导出 ---
# 建议修改导出文件名，以免覆盖之前的 Net0 结果
output_name <- "Combined_Norm0_Kinetics_Plot"
ggsave(paste0(output_name, ".pdf"), plot = p, width = 8, height = 6)
ggsave(paste0(output_name, ".png"), plot = p, width = 8, height = 6, dpi = 300)

