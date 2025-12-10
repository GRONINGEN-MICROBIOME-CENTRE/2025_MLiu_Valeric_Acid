
# GLM summary to p-value table
# ==============================================
glmToTbl <- function(fittedGLM) {
  gSmr <- summary(fittedGLM)
  gCoef <- as.data.frame(gSmr$coefficients)
  colnames(gCoef) <- c("Estimate","Std.Error","T.value","P.value")
  gCoef <- gCoef[row.names(gCoef) != "(Intercept)",]
  gCoef2 <- t(gCoef)
  colnames(gCoef2) <- paste0('pV.',colnames(gCoef2) )
  gCoef3 <- as_tibble(gCoef2)[4,]
}

# GLM summary to p-value table
# ==============================================
glmToTbl2 <- function(fittedGLM) {
  gSmr <- summary(fittedGLM)
  gCoef <- as.data.frame(gSmr$coefficients)
  colnames(gCoef) <- c("Estimate","Std.Error","T.value","P.value")
  gCoef <- gCoef[row.names(gCoef) != "(Intercept)",]
  gCoef2 <- t(gCoef)
  gCoef2 <- gCoef2[c("Estimate","P.value"),]
  gCoefPV <- as.tibble(gCoef2)[2,]
  colnames(gCoefPV) <- paste0(colnames(gCoefPV),'.pV' )
  gCoefEst <- as.tibble(gCoef2)[1,]
  colnames(gCoefEst) <- paste0(colnames(gCoefEst),'.Est' )
  gCoef3 <- cbind(gCoefPV,gCoefEst)
  gCoef3 <- gCoef3[,order(colnames(gCoef3))]
  gCoef3
}


# ==================================================
# ==================================================
#      biopsy Molecular Inflammation Score (bMIS)
# ==================================================
# ==================================================
library(readxl)
library(org.Hs.eg.db) # gene name converter
library(tidyverse)
library(ggplot2)
library(rstatix)
library(viridis) 
library(ggpubr)
library(readxl)
library(GSVA)

setwd('/Users/liumoting/R_Works/Valeric_acid/codes/')
# load normalized expression data
expNormalized <- readRDS('../data/Merged.normalized.RDS')

# > convert ENSEMBL gene names to gene symbols
genesAnnot <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys=colnames(expNormalized),
                                    keytype = "ENSEMBL",
                                    multiVals="first",
                                    columns=c("ENSEMBL","SYMBOL","GENENAME"))
genesAnnot <- genesAnnot[!duplicated(genesAnnot$ENSEMBL),]
colnames(genesAnnot) <- c("Gene","Gene.symbol","Annotation")

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


genesEnsToGetAll <- unique(c(genesEnsToGetIBD,genesEnsToGetCD,genesEnsToGetUC, genesEnsBarrier, genesEnsToGetExtra1))

# > subset genes to inBMIS
expNormalizedBMIS <- expNormalized[,genesEnsToGetAll]
expNormalizedBMISt <- as.data.frame(t(expNormalizedBMIS))
expNormalizedBMISt <- merge(x=expNormalizedBMISt,
                            y=genesAnnot[,c('Gene','Gene.symbol')],
                            by.x='row.names',by.y='Gene')
expNormalizedBMISt <- expNormalizedBMISt[!duplicated(expNormalizedBMISt$Gene.symbol),]
row.names(expNormalizedBMISt) <- expNormalizedBMISt$Gene.symbol
expNormalizedBMISt$Gene.symbol <- NULL; expNormalizedBMISt$Row.names <- NULL



# =======================
## ADD METADATA
# =======================
meta <- read.table('../data/RNAseqMetadata_cleaned_2023_02.csv',sep=',',header=T)
# > bit of cleaning
meta$Diagnosis <- as.factor(meta$Diagnosis)

meta$Diagnosis <- factor(meta$Diagnosis,levels = c('UC','CD'))
# merge
expNormalizedBMIStt <- as.data.frame(t(expNormalizedBMISt))
expNormalizedBMIStt$ID <- row.names(expNormalizedBMIStt)
countsrdy <- merge(expNormalizedBMIStt,meta,by="ID")
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

countsrdy$Barrier.score.Cld <- countsrdy$TJP1 + countsrdy$FABP2 + countsrdy$FABP6 - countsrdy$CD14 - countsrdy$LBP + countsrdy$CLDN1 + countsrdy$CLDN3
countsrdy$Barrier.score.noCld <- countsrdy$TJP1 + countsrdy$FABP2 + countsrdy$FABP6 - countsrdy$CD14 - countsrdy$LBP


countsrdyInf <- countsrdy

write.csv(countsrdyInf, file = "../results/Inflame_score/countsrdyInf.csv")











