
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
BiocManager::install("apeglm")
#install.packages("pheatmap")
BiocManager::install("org.Hs.eg.db") #change this if using a species other than human
#install.packages("ComplexHeatmap")
#install.packages("tibble")
#BiocManager::install("fgsea")
#install.packages("tidyverse")

library(readxl)
library(openxlsx)
library(tidyverse)
library(reshape2)
library(DESeq2)
library(apeglm)
library(pheatmap)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(tibble)
library(fgsea)
library(ggrepel)
library(ggfortify)


#Loading in data and metadata
setwd("/Users/alexis/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/CEMALB_DataAnalysisPM/Projects/P1016. THC Adducts/P1016.3. Analyses/P1016.1. DESeq2")
cur_date = '110526'
Output = "/Users/alexis/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/CEMALB_DataAnalysisPM/Projects/P1016. THC Adducts/P1016.3. Analyses/P1016.1. DESeq2/Output"
cts <- data.frame(read_excel("Input/RawCounts_D8_D8Q_16HBEs.xlsx", sheet = 2))
coldata <- data.frame(read_excel("Input/D8_D8Q_12_coldata.xlsx"))

#run the following lines and check that order of samples is the same in cts and coldata
head(cts)
coldata

#constructing a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition, tidy = TRUE)
dds

#pre-filtering genes
smallestGroupSize <- 3 #change this to the minimum number of samples you had in any of your treatment groups
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize #this requires genes to have 3 or more samples with counts of 10 or more. 
dds <- dds[keep,]

#specifying the reference level (control)
dds$condition <- relevel(dds$condition, ref = "VC")

#differential expression analysis for all groups combined 
dds <- DESeq(dds)

#transformed values
vsd <- vst(dds, blind=FALSE) #variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #regularized log transformation
ntd <- normTransform(dds) #normalized counts transformation
head(assay(vsd), 3)

#PCA of all groups
tiff(filename = "D8THC_D8THCQ_PCA_12.tiff", width = 1000, height = 2000, res = 300)
plotPCA(vsd, intgroup=c("condition")) #run this line only to visualize in R without saving
dev.off()


pcaData <- plotPCA(vsd, intgroup = "condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$condition <- factor(pcaData$condition, levels=c("VC", "D8THC", "D8THCQ"), labels=c("VC", "D8THC", "D8THCQ"))

#better looking PCA
png(filename = "D8THC_D8THCq_12hr_PCA.png", width = 5500, height = 4000, res = 900)
pca <- ggplot(pcaData, aes(PC1, PC2))
pca + geom_point(size=3, aes(color=condition, fill=condition, group=condition)) + 
  ggtitle("Principle Component Analysis 12 hours") + 
  xlab(paste0("PC1: ", percentVar[1],"% varience")) + 
  ylab(paste0("PC2: ", percentVar[2],"% varience")) + 
  theme(legend.title=element_blank()) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("VC" = "lightblue", "D8THC" = "#B59FDC", "D8THCQ" = "#6D2E46")) +
  coord_fixed() 
dev.off()


#Getting human gene symbols
columns(org.Hs.eg.db)

#Extracting results for each separate comparison
#comparison 1 (D8THC vs VC)
control_thc_df <- data.frame(results(dds, contrast = c("condition","D8THC","VC"))) #reference or control should be second
control_thc_df$Gene <- mapIds(org.Hs.eg.db, keys = row.names(control_thc_df), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
control_thc_df = control_thc_df[,c(7,2,5,6)] %>%
  rownames_to_column(var = 'ENSEMBL_ID') %>%
  relocate(ENSEMBL_ID, .before = Gene)
colnames(control_thc_df)[3:5] = c('log2FC', 'P Value', 'P Adj')

#comparison 2 (D8THCQ vs VC)
control_thcq_df <- data.frame(results(dds, contrast = c("condition","D8THCQ","VC")))
control_thcq_df$Gene <- mapIds(org.Hs.eg.db, keys = row.names(control_thcq_df), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
control_thcq_df = control_thcq_df[,c(7,2,5,6)] %>%
  rownames_to_column(var = 'ENSEMBL_ID') %>%
  relocate(ENSEMBL_ID, .before = Gene)
colnames(control_thcq_df)[3:5] = c('log2FC', 'P Value', 'P Adj')

#comparison 3 (D8THCQ vs D8THC)
thc_thcq_df <- data.frame(results(dds, contrast = c("condition","D8THCQ","D8THC")))
thc_thcq_df$Gene <- mapIds(org.Hs.eg.db, keys = row.names(thc_thcq_df), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
thc_thcq_df = thc_thcq_df[,c(7,2,5,6)] %>%
  rownames_to_column(var = 'ENSEMBL_ID') %>%
  relocate(ENSEMBL_ID, .before = Gene)
colnames(thc_thcq_df)[3:5] = c('log2FC', 'P Value', 'P Adj')

# exporting
write.xlsx(control_thc_df, paste0(Output,"/", "THC_Gene_Results_", cur_date, ".xlsx"), 
           rowNames = FALSE)
write.xlsx(control_thcq_df, paste0(Output,"/", "THCQ_Gene_Results_", cur_date, ".xlsx"), 
           rowNames = FALSE)
write.xlsx(thc_thcq_df, paste0(Output,"/", "THC_THCQ_Gene_Results_", cur_date, ".xlsx"), 
           rowNames = FALSE)