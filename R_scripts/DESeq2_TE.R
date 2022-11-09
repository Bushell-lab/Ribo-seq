#load libraries
library(tidyverse)

#read in common variables
source("common_variables.R")

#create a variable for what the treatment is----
treatment <- "KO"
control <- "WT"

#read in the most abundant transcripts per gene csv file----
most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

#read in data----
## Input files
# Calculating differential translation genes (DTGs) requires the count matrices from Ribo-seq and RNA-seq.
# These should be the raw counts, they should not be normalized or batch corrected.
# It also requires a sample information file which should be in the same order as samples in the count matrices.
# It should include information on sequencing type, treatment, batch or any other covariate you need to model.

#RPFs
#the following for loop reads in each final CDS counts file and renames the counts column by the sample name and saves each data frame to a list
data_list <- list()
for (sample in RPF_sample_names) {
  df <- read_csv(file = file.path(parent_dir, "Analysis/CDS_counts", paste0(sample, "_pc_final_counts_all_frames.csv")), col_names = T, col_types = c("c","i"))
  colnames(df) <- c("transcript", sample)
  data_list[[sample]] <- df
}

#merge all data within the above list using reduce
#Using full join retains all transcript but NAs need to be be replaced with 0 as these are transcipts that had 0 counts in that sample
#DESeq2 needs the transcripts to be as rownames, not as a column
data_list %>%
  reduce(full_join, by = "transcript") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  left_join(most_abundant_transcripts, by = "transcript") %>%
  relocate("gene") %>%
  select(-c("gene_sym", "transcript")) -> RPF_counts

#Totals
#load libraries
library(tximport)
library(DESeq2)

#read transcript to gene ID file and rename and reorder
read_tsv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_gene_IDs.txt", col_names = F) %>%
  dplyr::rename(GENEID = X1,
                TXNAME = X2) %>%
  select(TXNAME, GENEID) -> tx2gene

#import rsem data
#set directory where rsem output is located
rsem_dir <- file.path(parent_dir, 'rsem')

#create a named vector of files (with path)
files <- file.path(rsem_dir, paste0(Total_sample_names, ".isoforms.results"))
names(files) <- Total_sample_names

#import data with txi
txi <- tximport(files, type="rsem", tx2gene=tx2gene)

#make sure counts are as integers
as.data.frame(txi$counts) %>%
  mutate_all(~as.integer(.)) %>%
  rownames_to_column("gene") -> Total_counts

#create a data frame with the condition/replicate information----
#you need to make sure this data frame is correct for your samples
sample_info <- data.frame(row.names = c(RPF_sample_names, Total_sample_names),
                          Condition = factor(c(rep(control, 3), rep(treatment, 3))),
                          SeqType = factor(c(rep("RPFs", 6), rep("RNA", 6))),
                          replicate = factor(c(1:3,1:3)))

print(sample_info)

#merged RPF counts and Totals
RPF_counts %>%
  inner_join(Total_counts, by = "gene") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("gene") -> merged_data

summary(merged_data)

#check rownames of sample info are in same order as the column names of the merged data
table(colnames(merged_data) == row.names(sample_info))

## Detecting differential translation regulation
### DESeq2 object with batch and interaction term in the design

batch <- 1 #set to 1 if batch effect or 0 if not batch effect

if(batch == 1){
  ddsMat <- DESeqDataSetFromMatrix(countData = merged_data,
                                   colData = sample_info, design =~ replicate + Condition + SeqType + Condition:SeqType)
  
}else if(batch == 0){
  ddsMat <- DESeqDataSetFromMatrix(countData = merged_data,
                                   colData = sample_info, design =~ Condition + SeqType + Condition:SeqType)
}else{
  stop("Batch presence should be indicated by 0 or 1 only", call.=FALSE)
}

head(counts(ddsMat))
keep <- rowMeans(counts(ddsMat)) >= 10
table(keep)
ddsMat <- ddsMat[keep,]

#make sure seqtype and condition are set correctly
ddsMat$SeqType = relevel(ddsMat$SeqType,"RNA")
ddsMat$Condition = relevel(ddsMat$Condition,control)

#run DESeq
ddsMat <- DESeq(ddsMat)

# Choose the term you want to look at from resultsNames(ddsMat) 
resultsNames(ddsMat)

# paste0("Condition", treatment, ".SeqTypeRPFs") will test for a change in Ribo-seq levels in the treatment vs control
# accounting for changes in RNA-seq levels in treatment vs control
# alpha sets the adjusted p-value
res <- results(ddsMat, name=paste0("Condition", treatment, ".SeqTypeRPFs"), alpha = 0.1)
summary(res)

#output summary
sink(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("TE_", treatment, "_summary.txt")))
summary(res)
sink()

as.data.frame(res[order(res$padj),]) %>%
  rownames_to_column("gene") %>%
  inner_join(most_abundant_transcripts, by = "gene") -> TE_output
write_csv(TE_output, file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("TE_", treatment, "_DEseq2.csv")))

#output gene groups
TE_output %>%
  filter(padj < 0.1 & log2FoldChange < 0) %>%
  pull(gene_sym) -> A1_dep_geneIDs

TE_output %>%
  filter(padj < 0.1 & log2FoldChange > 0) %>%
  pull(gene_sym) -> A1_anti_geneIDs

write.table(file = file.path(parent_dir, "Analysis/transcript_IDs/A1_dep_gene_IDs.txt"), col.names = F, row.names = F, quote = F, A1_dep_geneIDs)
write.table(file = file.path(parent_dir, "Analysis/transcript_IDs/A1_antidep_gene_IDs.txt"), col.names = F, row.names = F, quote = F, A1_anti_geneIDs)
