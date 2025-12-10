library(readxl)
library(tidyverse)
library(ggplot2)
library(rstatix)
library(viridis) 
library(ggpubr)
library(magrittr)


setwd('/Users/liumoting/R_Works/Valeric_acid/codes/')
countsrdyInf_organoids <- read.csv("../results/Inflame_score/countsrdyInf_organoids.csv", header = TRUE, row.names = 1)
df <- countsrdyInf_organoids

# 因子化变量
df$Cytokines <- factor(df$Cytokines, levels = c("No", "Yes"))
df$Pretreat <- factor(df$Pretreat, levels = c("No", "Yes"))
# 模型中用连续变量建模
df$Concentration <- as.numeric(df$Concentration)
# 绘图时作为 factor 显示剂量点
df$Conc_factor <- factor(df$Concentration, levels = c(0, 0.3, 1, 3, 10))
df$Colon_organoids_cell_line <- factor(df$Colon_organoids_cell_line)
df$Timepoint <- factor(df$Timepoint)

# Blk VS CKs
df_sub_BlkCK <- df[df$Concentration == 0.0, ]

p_1 <- ggplot(df_sub_BlkCK, aes(x = Cytokines, y = oxKeNf, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Colon_organoids_cell_line), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("13" = "#E64B35", "18" = "#4DBBD5", "1045" = "#FDB462")) + 
  # facet_wrap(~ Colon_organoids_cell_line, scales = "free_y") +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("No", "Yes")),
    label = "p.signif",
    label.y = max(df_sub_BlkCK$oxKeNf, na.rm = TRUE) + 0.05
  ) +
  labs(
    x = "Cytokines",
    y = "Oxidative Score (Keap1/Nrf2)",
    title = "Oxidative Score by Cytokines and Cell Line",
    # fill = "Cytokine Treatment",
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  guides(fill = "none")
print(p_1)

ggsave("../results/Inflame_score/p_oxKeNf_boxplot_2.pdf", p_1, bg = "transparent", width = 10, height = 5, dpi = 300)
ggsave("../results/Inflame_score/p_oxKeNf_boxplot_3.pdf", p_1, bg = "transparent", width = 6, height = 5, dpi = 300)
