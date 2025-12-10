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

# ----------------------------------------------------------------------------------------
##########################################################################################
# visulization compare
# Keap1/Nrf2 oxidative stress
#########################

p_oxKeNf <- ggplot(df, aes(x = Conc_factor, y = oxKeNf, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Cytokines), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("No" = "#3C548880", "Yes" = "#00A08780")) + 
  facet_wrap(~ Colon_organoids_cell_line, scales = "free_y") +
  labs(
    x = "Valeric Acid (mM)",
    y = "oxidative stress (Keap1/Nrf2)",
    title = "OxKeap1/Nrf2 by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
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
print(p_oxKeNf)

ggsave("../results/Inflame_score/p_oxKeNf_boxplot.pdf", p_oxKeNf, bg = "transparent", width = 10, height = 5, dpi = 300)

# cell lines together
p_oxKeNf_1 <- ggplot(df, aes(x = Conc_factor, y = oxKeNf, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Colon_organoids_cell_line), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("13" = "#E64B35", "18" = "#4DBBD5", "1045" = "#FDB462")) + 
  labs(
    x = "Valeric Acid (mM)",
    y = "oxidative stress (Keap1/Nrf2)",
    title = "OxKeap1/Nrf2 by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
    color = "Organoid Cell Line"
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )
print(p_oxKeNf_1)

ggsave("../results/Inflame_score/p_oxKeNf_boxplot_1.pdf", p_oxKeNf_1, bg = "transparent", width = 10, height = 5, dpi = 300)


# visulization compare
# oxidative stress
#########################

p_oxScore2007 <- ggplot(df, aes(x = Conc_factor, y = oxScore2007, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Cytokines), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("No" = "#3C548880", "Yes" = "#00A08780")) + 
  facet_wrap(~ Colon_organoids_cell_line, scales = "free_y") +
  labs(
    x = "Valeric Acid (mM)",
    y = "oxScore2007",
    title = "OxScore2007 by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
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
print(p_oxScore2007)

ggsave("../results/Inflame_score/p_oxScore2007_boxplot.pdf", p_oxScore2007, bg = "transparent", width = 10, height = 5, dpi = 300)

# cell lines together
p_oxScore2007_1 <- ggplot(df, aes(x = Conc_factor, y = oxScore2007, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Colon_organoids_cell_line), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("13" = "#E64B35", "18" = "#4DBBD5", "1045" = "#FDB462")) + 
  labs(
    x = "Valeric Acid (mM)",
    y = "oxScore2007",
    title = "OxScore2007 by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
    color = "Organoid Cell Line"
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )
print(p_oxScore2007_1)

ggsave("../results/Inflame_score/p_oxScore2007_boxplot_1.pdf", p_oxScore2007_1, bg = "transparent", width = 10, height = 5, dpi = 300)

# Inflmae IBD score
#########################

p_bMISIBD <- ggplot(df, aes(x = Conc_factor, y = bMIS.IBD, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Cytokines), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("No" = "#3C548880", "Yes" = "#00A08780")) + 
  facet_wrap(~ Colon_organoids_cell_line, scales = "free_y") +
  labs(
    x = "Valeric Acid (mM)",
    y = "IBD Score",
    title = "IBD Score by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
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

# cell lines together
p_bMISIBD_1 <- ggplot(df, aes(x = Conc_factor, y = bMIS.IBD, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Colon_organoids_cell_line), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("13" = "#E64B35", "18" = "#4DBBD5", "1045" = "#FDB462")) + 
  labs(
    x = "Valeric Acid (mM)",
    y = "IBD Score",
    title = "IBD Score by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
    color = "Organoid Cell Line"
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )
print(p_bMISIBD_1)

ggsave("../results/Inflame_score/p_bMISIBD_boxplot_1.pdf", p_bMISIBD_1, bg = "transparent", width = 10, height = 5, dpi = 300)

### Inflmae CD score
#########################

p_bMISCD <- ggplot(df, aes(x = Conc_factor, y = bMIS.CD, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Cytokines), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("No" = "#3C548880", "Yes" = "#00A08780")) + 
  facet_wrap(~ Colon_organoids_cell_line, scales = "free_y") +
  labs(
    x = "Valeric Acid (mM)",
    y = "CD Score",
    title = "CD Score by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
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
print(p_bMISCD)

ggsave("../results/Inflame_score/p_bMIS_CD_boxplot.pdf", p_bMISCD, bg = "transparent", width = 10, height = 5, dpi = 300)

# cell lines together
p_bMISCD_1 <- ggplot(df, aes(x = Conc_factor, y = bMIS.CD, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Colon_organoids_cell_line), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("13" = "#E64B35", "18" = "#4DBBD5", "1045" = "#FDB462")) + 
  labs(
    x = "Valeric Acid (mM)",
    y = "CD Score",
    title = "CD Score by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
    color = "Organoid Cell Line"
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )
print(p_bMISCD_1)

ggsave("../results/Inflame_score/p_bMISCD_boxplot_1.pdf", p_bMISCD_1, bg = "transparent", width = 10, height = 5, dpi = 300)

### Inflmae UC score
#########################

p_bMISUC <- ggplot(df, aes(x = Conc_factor, y = bMIS.UC, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Cytokines), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("No" = "#3C548880", "Yes" = "#00A08780")) + 
  facet_wrap(~ Colon_organoids_cell_line, scales = "free_y") +
  labs(
    x = "Valeric Acid (mM)",
    y = "UC Score",
    title = "UC Score by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
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
print(p_bMISUC)

ggsave("../results/Inflame_score/p_bMIS_UC_boxplot.pdf", p_bMISUC, bg = "transparent", width = 10, height = 5, dpi = 300)

# cell lines together
p_bMISUC_1 <- ggplot(df, aes(x = Conc_factor, y = bMIS.UC, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Colon_organoids_cell_line), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("13" = "#E64B35", "18" = "#4DBBD5", "1045" = "#FDB462")) + 
  labs(
    x = "Valeric Acid (mM)",
    y = "UC Score",
    title = "UC Score by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
    color = "Organoid Cell Line"
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )
print(p_bMISUC_1)

ggsave("../results/Inflame_score/p_bMISUC_boxplot_1.pdf", p_bMISUC_1, bg = "transparent", width = 10, height = 5, dpi = 300)


# Barrier Score 
#########################
# with Claudins
p_barrier_Cld <- ggplot(df, aes(x = Conc_factor, y = Barrier.score.Cld, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Cytokines), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("No" = "#3C548880", "Yes" = "#00A08780")) + 
  facet_wrap(~ Colon_organoids_cell_line, scales = "free_y") +
  labs(
    x = "Valeric Acid (mM)",
    y = "Barrier Score",
    title = "Barrier Score by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
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
print(p_barrier_Cld)

ggsave("../results/Inflame_score/p_barrier_Cld_boxplot.pdf", p_barrier_Cld, bg = "transparent", width = 10, height = 5, dpi = 300)

# cell lines together
p_barrier_Cld_1 <- ggplot(df, aes(x = Conc_factor, y = Barrier.score.Cld, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Colon_organoids_cell_line), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("13" = "#E64B35", "18" = "#4DBBD5", "1045" = "#FDB462")) + 
  labs(
    x = "Valeric Acid (mM)",
    y = "Barrier Score",
    title = "Barrier Score by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
    color = "Organoid Cell Line"
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )
print(p_barrier_Cld_1)

ggsave("../results/Inflame_score/p_barrier_Cld_boxplot_1.pdf", p_barrier_Cld_1, bg = "transparent", width = 10, height = 5, dpi = 300)



# without Claudins
p_barrier_noCld <- ggplot(df, aes(x = Conc_factor, y = Barrier.score.noCld, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Cytokines), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("No" = "#3C548880", "Yes" = "#00A08780")) + 
  facet_wrap(~ Colon_organoids_cell_line, scales = "free_y") +
  labs(
    x = "Valeric Acid (mM)",
    y = "Barrier Score without Claudins",
    title = "Barrier Score (noClaudins) by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
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
print(p_barrier_noCld)

ggsave("../results/Inflame_score/p_barrier_noCld_boxplot.pdf", p_barrier_noCld, bg = "transparent", width = 10, height = 5, dpi = 300)

# cell lines together
p_barrier_noCld_1 <- ggplot(df, aes(x = Conc_factor, y = Barrier.score.noCld, fill = Cytokines)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  # Color jitter points by Colon_organoids_cell_line
  geom_jitter(aes(color = Colon_organoids_cell_line), width = 0.2, alpha = 0.6, size = 1.5) +
  # Custom fill colors for the box plots (Cytokines)
  scale_fill_manual(values = c("No" = "#3C5488", "Yes" = "#00A087")) + 
  scale_color_manual(values = c("13" = "#E64B35", "18" = "#4DBBD5", "1045" = "#FDB462")) + 
  labs(
    x = "Valeric Acid (mM)",
    y = "Barrier Score without Claudins",
    title = "Barrier Score (noClaudins) by Valeric Acid Concentration, Cytokines, and Cell Line",
    # fill = "Cytokine Treatment",
    color = "Organoid Cell Line"
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )
print(p_barrier_noCld_1)

ggsave("../results/Inflame_score/p_barrier_noCld_boxplot_1.pdf", p_barrier_noCld_1, bg = "transparent", width = 10, height = 5, dpi = 300)

#


# ----------------------------------------------------------------------------------------
#lm #########################################################################################
# lm 
model_ox <- lm(oxScore2007 ~ Cytokines * Concentration + Colon_organoids_cell_line + Pretreat + Timepoint, data = df)
summary(model_ox)

# 在不同细胞因子处理 + 不同戊酸浓度下的 adjusted marginal means
library(emmeans)
em_ox <- emmeans(model_ox, ~ Cytokines * Concentration, at = list(Concentration = c(0, 0.3, 1, 3, 10)))
summary(em_ox)
pairs(em_ox)
pairs(em_ox, adjust = "fdr")
contrast(em_ox, method = "revpairwise", by = "Cytokines")


library(ggplot2)

# 将 emmeans 输出转为 data.frame
df_em_ox <- as.data.frame(em_ox)

# 可视化交互边际均值
ggplot(df_em_ox, aes(x = as.factor(Concentration), y = emmean, group = Cytokines, color = Cytokines)) +
  geom_line(aes(linetype = Cytokines), position = position_dodge(0.2)) +
  geom_point(position = position_dodge(0.2), size = 2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = position_dodge(0.2)) +
#  facet_wrap(~ Colon_organoids_cell_line) +
  labs(
    title = "Adjusted oxScore2007 under Cytokine and VA Treatment",
    x = "Valeric Acid (mM)",
    y = "Adjusted oxScore2007",
    color = "Cytokines"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )


