library(dplyr)

# Ingore the part
input_data <- read.table("../data/Orga_bulk/deseq2_input.txt", header = TRUE, row.names = 1)
# Round the value
input_data <-round(input_data,digits = 0)
# Change the colnames to sample id
colnames(input_data) <- gsub("^X", "Orga", colnames(input_data))
colnames(input_data) <- gsub("_rsem\\.genes\\.results$", "", colnames(input_data)) 

# Transform to the expression matrix
input_data <- as.matrix(input_data)
# check the sample names are the same
all(rownames(metadata) %in% colnames(input_data))

# Make sure the order is the same
input_data <- input_data[, rownames(metadata)]

# Due to large difference in total counts, should choose TMM 
# save the count data
write.csv(input_data, file = "../results/Orga_bulkrna/count_data.csv", row.names = TRUE)


# Begin here----------------------------------------------------------------------------------------------------------------------
library(edgeR)

count <- read.csv("../results/Orga_bulkrna/count_data.csv", header = TRUE, row.names = 1)

count_data_clean <- count %>%
  dplyr::select(-Orga1045_7, -Orga13_26)

write.csv(count_data_clean, file = "../results/Orga_bulkrna/count_data_clean.csv", row.names = TRUE)
# Normalization by TMM method
dge <- DGEList(counts = count_data_clean)
dge <- calcNormFactors(dge, method = "TMM")
logCPM_counts_clean <- cpm(dge, log = TRUE, normalized.lib.sizes = TRUE)
# 保存edger标准化后的文件
write.csv(logCPM_counts_clean, file = "../results/Orga_bulkrna/normalized_counts_clean_logcpm.csv", row.names = TRUE)
