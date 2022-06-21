#load libraries
library(DESeq2)
library(tidyverse)
library(vsn)

#read in common variables
source("common_variables.R")

#create a variable for what the treatment is----
treatment <- "EFT226"

#read in the most abundant transcripts per gene csv file----
most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

#read in data----
#the following for loop reads in each final CDS counts file and renames the counts column by the sample name and saves each data frame to a list
data_list <- list()
for (sample in RPF_sample_names) {
  df <- read_csv(file = file.path(parent_dir, "Analysis/CDS_counts", paste0(sample, "_pc_final_counts_all_frames.csv")), col_names = T)
  colnames(df) <- c("transcript", sample)
  data_list[[sample]] <- df
}

#merge all data within the above list using reduce
#Using full join retains all transcripts, but the NAs need to be replaced with 0 as these are transcipts that had 0 counts in that sample
#DESeq2 needs the transcripts to be as rownames, not as a column
data_list %>%
  reduce(full_join, by = "transcript") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("transcript") -> RPF_counts

#create a data frame with the condition/batch information----
#you need to make sure this data frame is correct for your samples, the below creates one for a n=3 with EFT226 treatment.
sample_info <- data.frame(row_names = c('Ctrl_RPFs_1', 'Ctrl_RPFs_2', 'Ctrl_RPFs_3', 'EFT226_RPFs_1', 'EFT226_RPFs_2', 'EFT226_RPFs_3'),
                          condition = factor(c(rep("Ctrl", 3), rep(treatment, 3))),
                          batch = factor(c(rep(1:3, 2))))

#print the data frame to visually check it has been made as expected
sample_info

#make a DESeq data set from imported data----
DESeq2data <- DESeqDataSetFromMatrix(countData = RPF_counts,
                                     colData = sample_info,
                                     design = ~ batch + condition)

#pre-filter to remove genes with less than an average of 10 counts across all samples----
keep <- rowMeans(counts(DESeq2data)) >= 10
table(keep)
DESeq2data <- DESeq2data[keep,]

#make sure levels are set appropriately so that Ctrl is "untreated"
DESeq2data$condition <- relevel(DESeq2data$condition, ref = "Ctrl")

#run DESeq on DESeq data set----
dds <- DESeq(DESeq2data)

#extract results for each comparison----
res <- results(dds, contrast=c("condition", treatment, "Ctrl"))

#summarise results----
summary(res)

#apply LFC shrinkage for each comparison----
lfc_shrink <- lfcShrink(dds, coef=paste0("condition_", treatment, "_vs_Ctrl"), type="apeglm")

#write reslts to csv----
as.data.frame(lfc_shrink[order(lfc_shrink$padj),]) %>%
  rownames_to_column("transcript") %>%
  inner_join(most_abundant_transcripts, by = "transcript") -> DEseq2_output
write_csv(DEseq2_output, file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", treatment, "_DEseq2_apeglm_LFC_shrinkage.csv")))

#extract normalised counts and plot SD vs mean----
ntd <- normTransform(dds) #this gives log2(n + 1)
vsd <- vst(dds, blind=FALSE) #Variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #Regularized log transformation

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#Regularized log transformation looks preferable for this data. Check for your own data and select the appropriate one
#write out normalised counts data----
as.data.frame(assay(rld)) %>%
  rownames_to_column("transcript") %>%
  inner_join(most_abundant_transcripts, by = "transcript") -> normalised_counts
write_csv(normalised_counts, file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", treatment, "_normalised_counts.csv")))

#plot PCA----
pcaData <- plotPCA(rld, intgroup=c("condition", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png(filename = file.path(parent_dir, "plots/PCAs", paste0(treatment, "_RPFs_PCA.png")), width = 400, height = 350)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  ggtitle(paste(treatment, "RPFs"))
dev.off()

#apply batch correct and re-plot heatmap and PCA----
mat <- assay(rld)
mat <- limma::removeBatchEffect(mat, rld$batch)
assay(rld) <- mat

#PCA
pcaData <- plotPCA(rld, intgroup=c("condition", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png(filename = file.path(parent_dir, "plots/PCAs", paste0(treatment, "_RPFs_batch_corrected_PCA.png")), width = 400, height = 350)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  ggtitle(paste(treatment, "RPFs batch corrected"))
dev.off()


