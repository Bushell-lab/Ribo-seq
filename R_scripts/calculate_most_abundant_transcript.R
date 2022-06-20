#load libraries
library(DESeq2)
library(tidyverse)
library(tximport)

#read in common variables
source("common_variables.R")

#read in transcript and gene IDs
transcript_to_gene_ID <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_gene_IDs.csv", col_names = c("transcript", "gene", "gene_sym"))

#create a vector of rsem isoform file names
rsem_dir <- file.path(parent_dir, 'rsem')
files <- file.path(rsem_dir, paste0(Total_sample_names, ".isoforms.results"))
names(files) <- Total_sample_names

#read in rsem isoform files
isoform_expression <- tximport(files, type = "rsem", txOut = TRUE)

#calculate mean tpm
isoform_tpm <- as.data.frame(isoform_expression$abundance)
isoform_tpm$mean_tpm <- rowMeans(isoform_tpm)

#calculate most abundant transcript across all samples
isoform_tpm %>%
  rownames_to_column("transcript") %>%
  inner_join(transcript_to_gene_ID, by = "transcript") %>%
  group_by(gene) %>%
  top_n(n = 1, wt = mean_tpm) %>%
  sample_n(size = 1) %>% #some genes (particularly those with 0 counts) will have more than one transcript with the joint highest tpm, so we therefore select the transcript by random for these genes
  select(transcript, gene, gene_sym) -> most_abundant_transcripts

#check for no duplicates
nrow(most_abundant_transcripts) == n_distinct(most_abundant_transcripts$gene)

#write out as csv
write_csv(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"), most_abundant_transcripts)
write.table(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts.txt"), most_abundant_transcripts$transcript, row.names = F, col.names = F, quote = F)
