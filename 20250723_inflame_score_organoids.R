
library(readxl)
library(org.Hs.eg.db) # gene name converter
library(tidyverse)
library(ggplot2)
library(rstatix)
library(viridis) 
library(ggpubr)

setwd('/Users/liumoting/R_Works/Valeric_acid/codes/')

# load normalized expression data
expNormalized <- read.csv("../results/Orga_bulkrna/normalized_counts_clean.csv", header = TRUE, row.names = 1)
rownames(expNormalized) <- gsub("\\..*", "", rownames(expNormalized))
meta <- read.csv("../results/Orga_bulkrna/metadata_clean.csv", header = TRUE, row.names = 1)
# > convert ENSEMBL gene names to gene symbols
genesAnnot <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys=rownames(expNormalized),
                                    keytype = "ENSEMBL",
                                    multiVals="first",
                                    columns=c("ENSEMBL","SYMBOL","GENENAME"))

genesAnnot <- genesAnnot[!duplicated(genesAnnot$ENSEMBL),]
colnames(genesAnnot) <- c("Gene","Gene.symbol","Annotation")

### Save clean counts data with gene SYMBOLS and annotation
# ēxpNormalized <- expNormalized %>%
#   tibble::rownames_to_column(var = "Gene")
# expNormalized_with_symbols <- left_join(ēxpNormalized, genesAnnot, by = "Gene")
# 
# write.csv(expNormalized_with_symbols, file = "../results/Orga_bulkrna/normalized_counts_clean_with_symbols.csv")

# > load bMIS genes
inBMIS <- read_xlsx('../data/gut_2021_biopsy_blood_transcriptomics_disease_scores_2021.xlsx',1)

# > which genes are in our data? [IBD score]
genesToGetIBD <- inBMIS$bMIS.IBD.IvsNI[!is.na(inBMIS$bMIS.IBD.IvsNI)]
sum(genesToGetIBD %in% genesAnnot$Gene.symbol) # 279
genesToGetIBD <- genesToGetIBD[genesToGetIBD %in% genesAnnot$Gene.symbol]
genesEnsToGetIBD <- genesAnnot$Gene[genesAnnot$Gene.symbol %in% genesToGetIBD]
# > which genes are in our data? [UC score]
genesToGetUC <- inBMIS$bMIS.UC.IvsNI[!is.na(inBMIS$bMIS.UC.IvsNI)]
sum(genesToGetUC %in% genesAnnot$Gene.symbol) # 429
genesToGetUC <- genesToGetUC[genesToGetUC %in% genesAnnot$Gene.symbol]
genesEnsToGetUC <- genesAnnot$Gene[genesAnnot$Gene.symbol %in% genesToGetUC]
# > which genes are in our data? [CD score]
genesToGetCD <- inBMIS$bMIS.CD.IvsNI[!is.na(inBMIS$bMIS.CD.IvsNI)]
sum(genesToGetCD %in% genesAnnot$Gene.symbol) # 201
genesToGetCD <- genesToGetCD[genesToGetCD %in% genesAnnot$Gene.symbol]
genesEnsToGetCD <- genesAnnot$Gene[genesAnnot$Gene.symbol %in% genesToGetCD]

# > extra scores
ox2007genes <- c('CDKN1A','GDF15','PLK3','ATF3','TRP53INP1','DDIT4','GADD45A','BTG2','NDRG1')
genesEnsToGetExtra1 <- genesAnnot$Gene[genesAnnot$Gene.symbol %in% ox2007genes]

# barrier score
Barrier_genes <- c("LBP", "CD14", "FABP2", "FABP6", "CLDN1", "CLDN3", "TJP1")
genesBarrier <- Barrier_genes[Barrier_genes %in% genesAnnot$Gene.symbol]
genesEnsBarrier <- genesAnnot$Gene[genesAnnot$Gene.symbol %in% genesBarrier]
# "TJP1" %in% genesAnnot$Gene.symbol

# Keap1/Nrf2 Oxidative Score
oxKeapNrf <- read_xlsx('../data/keap1_nrf2_oxidative_markers.xlsx',1)
# > which genes are in our data? [keap1/nerf2 score]
genesToGetKeNf <- oxKeapNrf$Gene[!is.na(oxKeapNrf$Gene)]
sum(genesToGetKeNf %in% genesAnnot$Gene.symbol) # 63
genesToGetKeNf <- genesToGetKeNf[genesToGetKeNf %in% genesAnnot$Gene.symbol]
genesEnsToGetKeNf <- genesAnnot$Gene[genesAnnot$Gene.symbol %in% genesToGetKeNf]

genesEnsToGetAll <- unique(c(genesEnsToGetIBD,genesEnsToGetCD,genesEnsToGetUC, genesEnsBarrier, genesEnsToGetExtra1, genesEnsToGetKeNf))


# > subset genes to inBMIS
expNormalizedScore <- expNormalized[genesEnsToGetAll,]

expNormalizedScore <- merge(x=expNormalizedScore,
                            y=genesAnnot[,c('Gene','Gene.symbol')],
                            by.x='row.names',by.y='Gene')
expNormalizedScore <- expNormalizedScore[!duplicated(expNormalizedScore$Gene.symbol),]
row.names(expNormalizedScore) <- expNormalizedScore$Gene.symbol
expNormalizedScore$Gene.symbol <- NULL; expNormalizedScore$Row.names <- NULL

# =======================
## ADD METADATA
# =======================

# > bit of cleaning
meta$ID <- rownames(meta)

# merge
expNormalizedScoret <- as.data.frame(t(expNormalizedScore))
expNormalizedScoret$ID <- row.names(expNormalizedScoret)
countsrdy <- merge(expNormalizedScoret,meta,by="ID")
rownames(countsrdy) <- countsrdy$ID
# =======================
## CALCULATE SCORE(S)
# =======================
countsrdy$bMIS.IBD <- apply(MARGIN = 1,X = countsrdy[,colnames(countsrdy) %in% genesToGetIBD],
                            FUN = function(x) {log(sum(x))} )
countsrdy$bMIS.CD <- apply(MARGIN = 1,X = countsrdy[,colnames(countsrdy) %in% genesToGetCD],
                           FUN = function(x) {log(sum(x))} )
countsrdy$bMIS.UC <- apply(MARGIN = 1,X = countsrdy[,colnames(countsrdy) %in% genesToGetUC],
                           FUN = function(x) {log(sum(x))} )
countsrdy$oxScore2007 <- apply(MARGIN = 1,X = countsrdy[,colnames(countsrdy) %in% ox2007genes],
                               FUN = function(x) {log(sum(x))} )

countsrdy$oxKeNf <- apply(MARGIN = 1,X = countsrdy[,colnames(countsrdy) %in% genesToGetKeNf],
                               FUN = function(x) {log(sum(x))} )

countsrdy$Barrier.score.Cld <- countsrdy$TJP1 + countsrdy$FABP2 + countsrdy$FABP6 - countsrdy$CD14 - countsrdy$LBP + countsrdy$CLDN1 + countsrdy$CLDN3
countsrdy$Barrier.score.noCld <- countsrdy$TJP1 + countsrdy$FABP2 + countsrdy$FABP6 - countsrdy$CD14 - countsrdy$LBP


countsrdyInf_organoids <- countsrdy

write.csv(countsrdyInf_organoids, file = "../results/Inflame_score/countsrdyInf_organoids.csv")


