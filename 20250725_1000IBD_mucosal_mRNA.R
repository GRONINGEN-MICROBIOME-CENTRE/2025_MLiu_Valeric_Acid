
setwd('/Users/liumoting/R_works/Valeric_acid/codes/')
library(ggplot2)
library(rstatix)
library(tidyverse)
library(viridis) 
library(ggpubr)
library(readxl)

counts <- readRDS('/Users/liumoting/Library/CloudStorage/OneDrive-UMCG/2024_07_09_backup/R_works/Valeric_acid/data/Merged.normalized.RDS')
# load metadata
meta <- read.table('/Users/liumoting/Library/CloudStorage/OneDrive-UMCG/2024_07_09_backup/R_works/Valeric_acid/data/RNAseqMetadata_2022_09.csv', sep = ",", header = T)
# reality check
sum(meta$ID %in% rownames(counts))

#function####################################################

# 设定参数
mycolors <- c("#695F62", '#E6C857', '#996633', '#7E4598', '#330066') # 颜色
globalText <- 16  # 全局文本大小
dotSize <- 2  # 散点图点的大小
barSize <- 0.95  # 箱线图线条宽度

# 过滤显著性比较
filter_significant_comparisons <- function(my_comparisons, tstsSig) {
  significant_comparisons <- list()
  for (comp in my_comparisons) {
    if (any(tstsSig$group1 == comp[1] & tstsSig$group2 == comp[2])) {
      significant_comparisons <- append(significant_comparisons, list(comp))
    }
  }
  return(significant_comparisons)
}

# 单个基因的绘图和统计分析函数
plot_gene_expression <- function(gene, data, output_path) {
  # 设置 y 轴范围
  minE <- min(data[[gene]], na.rm = TRUE)
  maxE <- max(data[[gene]], na.rm = TRUE) * 1.3
  
  # Ileum 分析
  cntsIleum <- data[data$Location_rough == "ileum", ]
  cntsIleum$toTest <- cntsIleum[[gene]]
  
  g1 <- ggplot(cntsIleum, aes(x = Disease.Inflammation, y = !!sym(gene), col = Disease.Inflammation)) +
    geom_boxplot(outlier.alpha = 0, size = barSize, position = position_dodge(width = 0.2), na.rm = TRUE) +
    geom_jitter(alpha = 0.33, width = 0.15, height = 0.01, size = dotSize, na.rm = TRUE) +
    xlab(paste0(gene)) + ylab("Ileum") +
    scale_color_manual(values = mycolors, guide = "none") + 
    ylim(minE, maxE) +
    theme_minimal() +
    theme(legend.position = "none",
          text = element_text(size = globalText),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent", colour = NA),
          axis.line = element_line(color = "black", linewidth = 0.5))
  
  tsts <- compare_means(toTest ~ Disease.Inflammation, data = cntsIleum, p.adjust.method = "fdr")
  tstsSig <- tsts[tsts$p < 0.05, ]
  testIleum <- tsts
  testIleum$Location <- "Ileum"
  testIleum$Gene <- gene
  
  my_comparisons <- list(c("CD.NInf.", "CD.Inf."), c("UC.NInf.", "UC.Inf."), c("Ctrl", "CD.NInf."), 
                         c("Ctrl", "CD.Inf."), c("Ctrl", "UC.NInf."), c("Ctrl", "UC.Inf."))
  significant_comparisons <- filter_significant_comparisons(my_comparisons, tstsSig)
  
  g1 <- g1 + stat_compare_means(comparisons = significant_comparisons, label = "p.signif", 
                                hide.ns = TRUE, tip.length = 0, vjust = 0.5, step.increase = 0.08)
  
  # Colon 分析
  cntsColon <- data[data$Location_rough == "colon", ]
  cntsColon$toTest <- cntsColon[[gene]]
  
  g2 <- ggplot(cntsColon, aes(x = Disease.Inflammation, y = !!sym(gene), col = Disease.Inflammation)) +
    geom_boxplot(outlier.alpha = 0, size = barSize, position = position_dodge(width = 0.2), na.rm = TRUE) +
    geom_jitter(alpha = 0.33, width = 0.15, height = 0.01, size = dotSize, na.rm = TRUE) +
    xlab(paste0(gene)) + ylab("Colon") +
    scale_color_manual(values = mycolors, guide = "none") + 
    ylim(minE, maxE) +
    theme_minimal() +
    theme(legend.position = "none",
          text = element_text(size = globalText),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent", colour = NA),
          axis.line = element_line(color = "black", linewidth = 0.5))
  
  tsts <- compare_means(toTest ~ Disease.Inflammation, data = cntsColon, p.adjust.method = "fdr")
  tstsSig <- tsts[tsts$p < 0.05, ]
  testColon <- tsts
  testColon$Location <- "Colon"
  testColon$Gene <- gene
  
  significant_comparisons <- filter_significant_comparisons(my_comparisons, tstsSig)
  
  g2 <- g2 + stat_compare_means(comparisons = significant_comparisons, label = "p.signif", 
                                hide.ns = TRUE, tip.length = 0, vjust = 0.5, step.increase = 0.08)
  
  # 合并两个图
  g12 <- ggarrange(g1, g2)
  print(g12)
  
  # 保存图像
  ggsave(plot = g12, filename = paste0(output_path, "/plot_", gene, ".png"),
         width = 9, height = 6, bg = 'white')
  
  # 结果数据
  resGenes <- rbind.data.frame(testIleum, testColon)
  return(resGenes)
}


# KEAP1 ENSG00000079999
# NFE2L2 ENSG00000116044
# 1. 提取感兴趣基因的表达数据
counts_CASP <- counts[, colnames(counts) == "ENSG00000116044", drop = FALSE]  # 只保留该基因

# 2. 确保行名（样本 ID）匹配 meta 数据
counts_CASP <- counts_CASP[rownames(counts_CASP) %in% meta$ID, , drop = FALSE]

# 3. 转换为数据框，并添加 ID 列
counts_CASP <- as.data.frame(counts_CASP)
counts_CASP$ID <- rownames(counts_CASP)

# 4. 合并 `meta` 数据
countsrdy_CASP <- merge(counts_CASP, meta, by = "ID")

# 5. 进一步清理（去掉 Inflammation == "Light" 的样本）
countsrdy_CASP <- countsrdy_CASP[countsrdy_CASP$Inflammation != "Light", ]

# 6. 修改列名，使基因列名为 `WDR1`
colnames(countsrdy_CASP)[colnames(countsrdy_CASP) == "ENSG00000116044"] <- "NFE2L2"

# 复制数据框，避免修改原数据
countsrdy_CASP_2 <- countsrdy_CASP

# 1. 重新编码 Inflammation 列
countsrdy_CASP_2$Inflammation[countsrdy_CASP_2$Inflammation == "Yes"] <- "Inf."
countsrdy_CASP_2$Inflammation[countsrdy_CASP_2$Inflammation == "No"] <- "NInf."

# 2. 转换 Inflammation 为因子类型
countsrdy_CASP_2$Inflammation <- as.factor(countsrdy_CASP_2$Inflammation)

# 3. 创建 `Disease.Inflammation` 变量（基于 `Diagnosis` 和 `Inflammation`）
countsrdy_CASP_2$Disease.Inflammation <- paste0(as.character(countsrdy_CASP_2$Diagnosis), '.', as.character(countsrdy_CASP_2$Inflammation))

# 4. 将 'Control.NInf.' 重命名为 'Ctrl'
countsrdy_CASP_2$Disease.Inflammation[countsrdy_CASP_2$Disease.Inflammation == 'Control.NInf.'] <- 'Ctrl'

# 5. 转换 `Disease.Inflammation` 为因子
countsrdy_CASP_2$Disease.Inflammation <- as.factor(countsrdy_CASP_2$Disease.Inflammation)

# 6. 设置 `Disease.Inflammation` 因子的水平顺序
countsrdy_CASP_2$Disease.Inflammation <- factor(countsrdy_CASP_2$Disease.Inflammation,
                                                levels = c('Ctrl', 'CD.NInf.', 'CD.Inf.', 'UC.NInf.', 'UC.Inf.'))

# 查看最终结果
head(countsrdy_CASP_2)


gene_name <- "NFE2L2"  # 这里替换成你想要分析的基因
output_dir <- "../results/human_1000IBD/plots"

resGenes <- plot_gene_expression(gene_name, countsrdy_CASP_2, output_dir)

# **保存结果到 CSV**
# write.table(resGenes, file.path(output_dir, "genes_alltests3_Jessica.csv"), sep = ',', row.names = FALSE)


# Get all tests
get_gene_stats <- function(gene, data) {
  # 提取 ileum 和 colon 数据
  ileum_data <- data[data$Location_rough == "ileum", ]
  colon_data <- data[data$Location_rough == "colon", ]
  
  ileum_data$toTest <- ileum_data[[gene]]
  colon_data$toTest <- colon_data[[gene]]
  
  # 设定比较组
  my_comparisons <- list(
    c("CD.NInf.", "CD.Inf."), 
    c("UC.NInf.", "UC.Inf."), 
    c("Ctrl", "CD.NInf."), 
    c("Ctrl", "CD.Inf."), 
    c("Ctrl", "UC.NInf."), 
    c("Ctrl", "UC.Inf.")
  )
  
  # 所有两两比较
  testIleum <- compare_means(toTest ~ Disease.Inflammation, 
                             data = ileum_data, 
                             method = "wilcox.test", 
                             p.adjust.method = "fdr",
                             comparisons = my_comparisons)
  testIleum$Location <- "Ileum"
  testIleum$Gene <- gene
  
  testColon <- compare_means(toTest ~ Disease.Inflammation, 
                             data = colon_data, 
                             method = "wilcox.test", 
                             p.adjust.method = "fdr",
                             comparisons = my_comparisons)
  testColon$Location <- "Colon"
  testColon$Gene <- gene
  
  # 合并结果
  all_results <- rbind(testIleum, testColon)
  return(all_results)
}


# 获取统计结果
stats_df <- get_gene_stats(gene_name, countsrdy_CASP_2)

# 保存成 CSV
write.table(stats_df, file.path(output_dir, paste0("stats_", gene_name, ".csv")), sep = ",", row.names = FALSE)
