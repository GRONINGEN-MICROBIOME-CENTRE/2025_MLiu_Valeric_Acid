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


df$Cytokines <- factor(df$Cytokines, levels = c("No", "Yes"))
df$Pretreat <- factor(df$Pretreat, levels = c("No", "Yes"))

df$Concentration <- as.numeric(df$Concentration)

df$Conc_factor <- factor(df$Concentration, levels = c(0, 0.3, 1, 3, 10))
df$Colon_organoids_cell_line <- factor(df$Colon_organoids_cell_line)
df$Timepoint <- factor(df$Timepoint)


df_sub_BlkCK <- df[df$Concentration == 0.0, ]


p_bMISIBD <- ggplot(df_sub_BlkCK, aes(x = Conc_factor, y = bMIS.IBD, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  
  geom_jitter(aes(color = Cytokines), width = 0.2, alpha = 0.6, size = 1.5) +
  
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("No" = "#3C548880", "Yes" = "#00A08780")) + 

  facet_wrap(~ Colon_organoids_cell_line, scales = "free_y") +
  labs(
    x = "Valeric Acid (mM)",
    y = "IBD Score",
    title = "IBD Score by Valeric Acid Concentration, Cytokines, and Cell Line",
    
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
print(p_bMISIBD)

ggsave("../results/Inflame_score/p_bMIS_IBD_boxplot.pdf", p_bMISIBD, bg = "transparent", width = 10, height = 5, dpi = 300)

