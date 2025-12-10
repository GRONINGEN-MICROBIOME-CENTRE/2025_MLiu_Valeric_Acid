setwd("../R_works/Valeric_acid/codes")

library(ggplot2)
library(rstatix)
library(tidyverse)
library(viridis)
library(ggpubr)  
library(readxl)



counts <- readRDS('../data/RNAseq_Merged_TMM_logCPM_ComBat_batchlocinflam.RDS')

meta <- read.table('../data/RNAseqMetadata_2022_09.csv', sep = ",", header = T)
sum(meta$ID %in% rownames(counts)) 

genes <- read_excel('../data/Orga_bulk/VA_organoids_barrier.xlsx')

sum(genes$Ensembl_ID %in% colnames(counts))  
sum(!genes$Ensembl_ID %in% colnames(counts)) 

nf <- genes$Ensembl_ID[!(genes$Ensembl_ID %in% colnames(counts))]
fnd <- genes$Ensembl_ID[(genes$Ensembl_ID %in% colnames(counts))]
genes[genes$Ensembl_ID %in% nf, ]
genes[genes$Ensembl_ID %in% fnd, ]


counts2 <- counts[, colnames(counts) %in% genes$Ensembl_ID]
counts2 <- counts2[rownames(counts2) %in% meta$ID, ] 
genes <- as.data.frame(genes)

for (cc in c(1:ncol(counts2))) {
  cc2 = colnames(counts2)[cc]
  cc3 = genes$Genes[genes$Ensembl_ID == cc2]
  colnames(counts2)[cc] = cc3
}

genes <- colnames(counts2)

counts3 <- as.data.frame(counts2)
counts3$ID <- rownames(counts3)
countsrdy <- merge(counts3, meta, by = "ID") 

countsrdy <- countsrdy[countsrdy$Inflammation != "Light", ] 


genesToTest <- genes
countsrdy2 <- countsrdy
countsrdy2$Inflammation[countsrdy2$Inflammation=="Yes"] <- "Inf."
countsrdy2$Inflammation[countsrdy2$Inflammation=="No"] <- "NInf."
countsrdy2$Inflammation <- as.factor(countsrdy2$Inflammation)
countsrdy2$Disease.Inflammation <- paste0(as.character(countsrdy2$Diagnosis), '.', as.character(countsrdy2$Inflammation))
countsrdy2$Disease.Inflammation[countsrdy2$Disease.Inflammation == 'Control.NInf.'] <- 'Ctrl'
countsrdy2$Disease.Inflammation <- as.factor(countsrdy2$Disease.Inflammation)
countsrdy2$Disease.Inflammation <- factor(countsrdy2$Disease.Inflammation,
                                          levels = c('Ctrl', 'CD.NInf.', 'CD.Inf.', 'UC.NInf.', 'UC.Inf.'))


table(countsrdy2$Disease.Inflammation[countsrdy2$Location_rough == "colon"])
table(countsrdy2$Disease.Inflammation[countsrdy2$Location_rough == "ileum"])



getAllComparisons <- function(.tbl){
  .tbl %>%
    dplyr::select(.data$group1, .data$group2) %>%
    purrr::transpose() %>%
    purrr::modify_depth(1, unlist)
}


resGenes <- NULL
for (gene in genesToTest) {
  mycolors <-c('blue', '#f0b618', '#ff470a', '#7d52ff', '#ff0586') 
  globalText <- 16 
  dotSize <- 2 
  barSize <- 0.95 
  minE <- min(countsrdy2[[gene]]) 
  maxE <- max(countsrdy2[[gene]])*1.5 
  
  cntsIleum <- countsrdy2[countsrdy2$Location_rough == "ileum", ]
  g1 <- ggplot(cntsIleum,
               aes_string(x = "Disease.Inflammation",
                          y = gene,
                          col = "Disease.Inflammation")) +
    geom_boxplot(outlier.alpha = 0, size = barSize, na.rm = T) + 
    geom_jitter(alpha = 0.33, width = 0.15, height = 0.01, size = dotSize, na.rm = T) + 
    xlab(paste0(gene)) + ylab("Ileum") +
    scale_color_manual(values = mycolors) + 
    ylim(minE, maxE) +
    theme_minimal() + 
    theme(legend.position = "none") + 
    theme(text = element_text(size = globalText))
  
  cntsIleum$toTest <- cntsIleum[[gene]]
  tsts <- compare_means(toTest ~ Disease.Inflammation, 
                        data = cntsIleum,
                        p.adjust.method = "fdr")
  tstsSig <- tsts[tsts$p < 0.05, ]
  testIleum <- tsts
  testIleum$Location <- "Ileum"
  testIleum$Gene <- gene
  my_comparisons <- getAllComparisons(tstsSig)
  g1 <- g1 + stat_compare_means(comparisons = my_comparisons,
                                label = "p.signif",
                                hide.ns = T)
  
  cntsColon <- countsrdy2[countsrdy2$Location_rough == "colon", ]
  g2 <- ggplot(cntsColon,
               aes_string(x = "Disease.Inflammation",
                          y = gene,
                          col = "Disease.Inflammation")) +
    geom_boxplot(outlier.alpha = 0, size = barSize, na.rm = T) + 
    geom_jitter(alpha = 0.33, width = 0.15, height = 0.01, size = dotSize, na.rm = T) +
    xlab(paste0(gene)) + ylab("Colon") +
    scale_color_manual(values = mycolors) + 
    ylim(minE, maxE) +
    theme_minimal() +
    theme(legend.position = "none") +
    theme(text = element_text(size = globalText))
  
  cntsColon$toTest <- cntsColon[[gene]]
  tsts <- compare_means(toTest ~ Disease.Inflammation, 
                        data = cntsColon,
                        p.adjust.method = "fdr")
  tstsSig <- tsts[tsts$p < 0.05, ]
  testColon <- tsts
  testColon$Location <- "Colon"
  testColon$Gene <- gene
  
  my_comparisons <- getAllComparisons(tstsSig)
  g2 <- g2 + stat_compare_means(comparisons = my_comparisons,
                                label = "p.signif",
                                hide.ns = T)
  
  g12 <- ggarrange(g1, g2)
  print(g12)
  ggsave(plot = g12, paste0('../results/human_1000IBD/VA_0810/plots_barrier/plot_', gene, '.pdf'),
         width = 9, height = 6, bg = 'white')
  resGenes <- rbind.data.frame(resGenes, testIleum)
  resGenes <- rbind.data.frame(resGenes, testColon)
}

resGenes$FDR <- p.adjust(resGenes$p)
write.table(resGenes, '../results/human_1000IBD/VA_0810/plots_barrier/genes_alltests.csv', sep = ',', row.names = F)


############################### Plot only Colon data ###########################

library(dplyr)
library(ggplot2)
library(ggpubr)
library(purrr)



my_comparisons <- list(
  c("CD.NInf.", "CD.Inf."),
  c("UC.NInf.", "UC.Inf."),
  c("Ctrl", "CD.NInf."),
  c("Ctrl", "CD.Inf."),
  c("Ctrl", "UC.NInf."),
  c("Ctrl", "UC.Inf.")
)

resGenes <- NULL
level_order <- c("Ctrl", "CD.NInf.", "CD.Inf.", "UC.NInf.", "UC.Inf.")  


filter_significant_comparisons <- function(comparisons, significant_tests) {
  significant_comparisons <- list()
  for (comp in comparisons) {
    
    if (any(significant_tests$group1 == comp[1] & significant_tests$group2 == comp[2])) {
      significant_comparisons[[length(significant_comparisons) + 1]] <- comp
    }
  }
  return(significant_comparisons)
}


for (gene in genesToTest) {
  mycolors <- c('#8B9E68', '#5499C7', '#1F618D', '#C2736C', '#C0382B')
  globalText <- 16
  dotSize <- 2
  barSize <- 0.95
  
  
  cntsColon <- countsrdy2[countsrdy2$Location_rough == "colon", ]
  if (!gene %in% colnames(cntsColon)) next
  if (all(is.na(cntsColon[[gene]]))) next
  
  
  cntsColon$Disease.Inflammation <- factor(
    cntsColon$Disease.Inflammation,
    levels = level_order
  )
  cntsColon <- droplevels(cntsColon[cntsColon$Disease.Inflammation %in% unique(unlist(my_comparisons)), ])
  
  
  minE <- min(cntsColon[[gene]], na.rm = TRUE)
  maxE <- max(cntsColon[[gene]], na.rm = TRUE)
  rng  <- ifelse(is.finite(maxE - minE) && (maxE - minE) > 0, maxE - minE, 1)
  
  base <- maxE + 0.05 * rng       
  step <- 0.08 * rng              
  ypos <- base + step * seq_along(my_comparisons)
  upper <- tail(ypos, 1) + 0.06 * rng  
  
  
  cntsColon$toTest <- cntsColon[[gene]]
  tsts <- ggpubr::compare_means(
    toTest ~ Disease.Inflammation,
    data = cntsColon,
    method = "wilcox.test",
    p.adjust.method = "fdr"
  )
  
  tsts$Location <- "Colon"
  tsts$Gene <- gene
  
  tstsSig <- tsts[tsts$p < 0.05, ]
  significant_comparisons <- filter_significant_comparisons(my_comparisons, tstsSig)
  
  
  g2 <- ggplot(
    cntsColon,
    aes_string(x = "Disease.Inflammation", y = gene, col = "Disease.Inflammation")
  ) +
    geom_boxplot(outlier.alpha = 0, size = barSize, na.rm = TRUE) +
    geom_jitter(alpha = 0.33, width = 0.15, height = 0.01, size = dotSize, na.rm = TRUE) +
    xlab(gene) + ylab("Gene Expression Level") +
    scale_color_manual(values = mycolors) +
    theme_minimal() +
    theme(
      legend.position = "none",
      text = element_text(size = globalText),
      plot.margin = margin(t = 14, r = 14, b = 12, l = 12),
      axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank()
    ) +
    
    coord_cartesian(ylim = c(minE, upper), clip = "off")
  
  
  g2 <- g2 + stat_compare_means(
    comparisons = significant_comparisons,
    label = "p.signif",
    tip.length = 0,
    method = "wilcox.test",
    p.adjust.method = "fdr",
    y.position = ypos,
    size = 5,
  )
  
  print(g2)
  ggsave(plot = g2,
         filename = paste0('../results/human_1000IBD/VA_0810/plots_barrier/plot_colon_', gene, '.pdf'),
         width = 6, height = 7, bg = 'white')
  
  resGenes <- rbind(resGenes, tsts)
}


resGenes$FDR <- p.adjust(resGenes$p)
write.table(resGenes, '../results/human_1000IBD/VA_0810/plots_barrier/genes_colon_barrier_comparisons.csv', sep = ',', row.names = F)

