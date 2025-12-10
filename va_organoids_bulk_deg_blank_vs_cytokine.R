
library(dplyr)
library(readxl)
library(readr)
library(reshape2)
library(ggrepel)
library(ggplot2)
library(plotly)
library(matrixStats)
library(pheatmap)
library(edgeR)
library(vegan)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

setwd("/Users/liumoting/R_Works/Valeric_acid/codes/")
#-------------------------------------------------------------------------------------------------------
count_data <- read.csv("../results/Orga_bulkrna/count_data_clean.csv", header = TRUE, row.names = 1)
rownames(count_data) <- gsub("\\..*", "", rownames(count_data))

metadata_BlkCK <- read.csv("../results/Orga_bulkrna/metadata_BlkCK.csv", header = TRUE, row.names = 1)

count_data_BlkCK <- count_data[, rownames(metadata_BlkCK), drop = FALSE]

# ----------------------------------------------------------------------------------------------------------------------
# DEG
metadata_BlkCK_6h <- subset(metadata_BlkCK, Timepoint == 6)
count_data_BlkCK_6h <- count_data_BlkCK[, rownames(metadata_BlkCK_6h)]
# DGE
metadata <- metadata_BlkCK_6h
count_data <- count_data_BlkCK_6h

rownames(metadata) <- colnames(count_data)

# Create DGEList object
y <- DGEList(counts = count_data)

keep <- filterByExpr(y, design = model.matrix(~ Colon_organoids_cell_line + Cytokines, data = metadata))
y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y, method = "TMM")

metadata$Colon_organoids_cell_line <- as.factor(metadata$Colon_organoids_cell_line)
metadata$Cytokines <- as.factor(metadata$Cytokines)

design <- model.matrix(~ Colon_organoids_cell_line + Cytokines, data=metadata)

print(design)

y <- estimateDisp(y, design, robust = TRUE)
fit <- glmFit(y, design, robust = TRUE)
contrast <- makeContrasts(CytokinesYes, levels=design)

fit2 <- glmLRT(fit, contrast=contrast)

deg_results <- topTags(fit2, n=Inf)

deg_results_df <- as.data.frame(deg_results)

write.csv(deg_results_df, file = "../results/Orga_bulkrna/VA_0810/DEG_results_Cytokines_6h.csv", row.names=TRUE)

# ----------------------------------------------------------------------------------------------------------------------
logCPM <- cpm(y, log=TRUE, prior.count=1)
pca <- prcomp(t(logCPM))

pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], 
                     CellLine = metadata$Colon_organoids_cell_line, 
                     Cytokines = metadata$Cytokines)

pp <- ggplot(pca_df, aes(x = PC1, y = PC2, color = CellLine, shape = Cytokines)) +
  geom_point(size = 5, alpha = 0.7, stroke = 1.5) +  # Increased size, transparency, and stroke width
  scale_color_manual(values = c("13" = "#E69F00", "18" = "#56B4E9", "1045" = "#009E73")) +  # Custom colors for CellLine
  scale_shape_manual(values = c(16, 17)) +  # Different shapes for Cytokines
  theme_minimal(base_size = 15) +  # Larger base font size for better readability
  theme(legend.position = "right",  # Place legend on the right
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Centered and bold title
        axis.title = element_text(size = 14),  # Larger axis titles
        axis.text = element_text(size = 12),  # Larger axis labels
        panel.grid.minor = element_blank()) +  # Remove minor grid lines
  labs(title = "PCA of Colon Organoids RNA-seq Data (6h)", 
       x = "PCA 1", 
       y = "PCA 2")  # Add axis labels
print(pp)

ggsave(filename = "../results/Orga_bulkrna/VA_0810/01_BlkCK/PCA_BlkCK_6h.pdf", plot = pp, width = 8, height = 6, dpi = 300)

### Volcano
deg_results_df$Gene <- rownames(deg_results_df)
gene_mapping <- bitr(deg_results_df$Gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)

deg_results_df <- left_join(deg_results_df, gene_mapping, by=c("Gene"="ENSEMBL"))

deg_results_df <- deg_results_df %>% mutate(SYMBOL = ifelse(is.na(SYMBOL), Gene, SYMBOL))

sig_genes <- deg_results_df %>% 
  filter(FDR < 0.05 & abs(logFC) > 1) %>% 
  dplyr::select(SYMBOL, logFC, FDR)

pv <- ggplot(deg_results_df, aes(x = logFC, y = -log10(FDR), color = case_when(
  FDR < 0.05 & logFC > 0 ~ "Upregulated",   # Upregulated (logFC > 0, FDR < 0.05)
  FDR < 0.05 & logFC < 0 ~ "Downregulated", # Downregulated (logFC < 0, FDR < 0.05)
  TRUE ~ "Not significant"                   # Non-significant (FDR >= 0.05)
))) +
  geom_point(alpha = 0.8, size = 2.5, shape = 16) +  # Adjusted point size and shape
  theme_minimal(base_size = 14) +  # Larger base font for better readability
  scale_color_manual(values = c("Upregulated" = "#F09148", "Downregulated" = "#456990", "Not significant" = "gray")) +  # Adjusted colors for better contrast
  labs(title = "Volcano Plot", 
       x = "log2 Fold Change", 
       y = "-log10(FDR)") +
  theme(legend.position = "none",  # Remove legend
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centered and bold title
        plot.subtitle = element_text(hjust = 0.5, size = 12),  # Subtitle size
        axis.title = element_text(size = 14),  # Axis titles size
        axis.text = element_text(size = 12),  # Axis text size
        #panel.grid.major = element_blank(),  # Removed major grid lines
        panel.grid.minor = element_blank()) +  # Removed minor grid lines
  geom_text_repel(data = sig_genes, aes(label = SYMBOL), 
                  size = 4, 
                  box.padding = 0.35, 
                  min.segment.length = 1.5, 
                  color = "black",  # Color for labels
                  fontface = "bold",  # Bold labels
                  max.overlaps = 10)  # Limit number of overlapping labels

# Print the plot
print(pv)

ggsave(filename = "../results/Orga_bulkrna/VA_0810/01_BlkCK/Volcano_BlkCK_6h.pdf", plot = pv, width = 6, height = 6, dpi = 300)
### Heatmap
sig_genes <- deg_results_df[deg_results_df$FDR < 0.05, ]$Gene

heatmap_data <- logCPM[sig_genes, ]

ph <- pheatmap(heatmap_data, annotation_col=metadata %>% dplyr::select(-Timepoint), show_rownames=FALSE)

ggsave(filename = "../results/Orga_bulkrna/VA_0810/01_BlkCK/Heatmap_BlkCK_6h.pdf", plot = ph, width = 8, height = 6, dpi = 300)

# ----------------------------------------------------------------------------------------------------------------------
### Enrichment
deg_threshold <- 0.05  # FDR threshold
logFC_threshold <- 1   # log2FC threshold

up_genes <- deg_results_df %>%
  filter(FDR < deg_threshold & logFC > logFC_threshold) %>%
  pull(SYMBOL)

down_genes <- deg_results_df %>%
  filter(FDR < deg_threshold & logFC < -logFC_threshold) %>%
  pull(SYMBOL)

all_degs <- c(up_genes, down_genes)

### KEGG
all_genes_entrez <- bitr(all_degs, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

kegg_res <- enrichKEGG(gene = all_genes_entrez$ENTREZID,
                       organism = "hsa",  
                       pvalueCutoff = 0.05)

# KEGG dot plot for all DEGs
pkegg <- dotplot(kegg_res, showCategory=20, title="KEGG Pathway Enrichment")

print(pkegg)

# Upregulated genes
up_genes_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

kegg_res_up <- enrichKEGG(gene = up_genes_entrez$ENTREZID,
                          organism = "hsa",
                          pvalueCutoff = 0.05)

# KEGG dot plot for upregulated genes
pkegg_up <- dotplot(kegg_res_up, showCategory=20, title="KEGG Pathway Enrichment (Up)")

print(pkegg_up)

ggsave(filename = "../results/Orga_bulkrna/VA_0810/01_BlkCK/KEGGup_BlkCK_6h.pdf", plot = pkegg_up, width = 8, height = 8, dpi = 300)
# Downregulated genes
down_genes_entrez <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

kegg_res_down <- enrichKEGG(gene = down_genes_entrez$ENTREZID,
                            organism = "hsa",
                            pvalueCutoff = 0.2)

# KEGG dot plot for downregulated genes
pkegg_down <- dotplot(kegg_res_down, showCategory=20, title="KEGG Pathway Enrichment (Down)")
print(pkegg_down)

ggsave(filename = "../results/Orga_bulkrna/VA_0810/01_BlkCK/KEGGdown_BlkCK_6h.pdf", plot = pkegg_down, width = 8, height = 8, dpi = 300)
