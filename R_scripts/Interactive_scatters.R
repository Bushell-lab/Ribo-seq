#load libraries
library(tidyverse)
library(Glimma)

#read in common variables
source("common_variables.R")

#create a variable for what the treatment is----
control <- "WT"
treatment <- "KO"

#read in DESeq2 output----
read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2.csv")) %>%
  filter(!is.na(TE_group)) %>%
  mutate(TE_group = factor(TE_group),
         groupings = case_when(TE_group == "TE down" ~ -1,
                               TE_group == "TE up" ~ 1,
                               TE_group == "no change" | TE_group == "NS" ~ 0 )) -> DESeq2_data
  
summary(DESeq2_data)

RPFs_norm_counts <- read.csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", treatment, "_normalised_counts.csv")))
totals_norm_counts <- read.csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", treatment, "_normalised_counts.csv")))

#make gene ID/sym annotation table
DESeq2_data %>%
  select(gene, gene_sym) %>%
  dplyr::rename(GeneID = gene) %>%
  as.data.frame() -> gene_anno
rownames(gene_anno) <- gene_anno$GeneID

#merge norm counts with gene annotation so that row names are in the same order and then select all counts columns in preferred order and make geneIDs row names
gene_anno %>%
  inner_join(RPFs_norm_counts, by = c("GeneID" = "gene", "gene_sym")) %>%
  inner_join(totals_norm_counts, by = c("GeneID" = "gene", "gene_sym", "transcript")) %>%
  column_to_rownames("GeneID") %>%
  select(-c("transcript", "gene_sym")) -> merged_norm_counts

#check row names match up
all(rownames(gene_anno) == DESeq2_data$gene)
all(rownames(gene_anno) == row.names(merged_norm_counts))

#make html files
glXYPlot(x = DESeq2_data$totals_log2FC,
         y = DESeq2_data$RPFs_log2FC,
         xlab = "RNA log2FC",
         ylab = "RPFs log2FC",
         main = treatment,
         status = DESeq2_data$groupings,
         counts = merged_norm_counts,
         side.xlab = "Sample",
         side.ylab = "norm counts",
         sample.cols = rep(c("#F8766D", "#00BA38", "#619CFF"),4), #these are the colours for each replicate within the norm counts, vector needs to be the same length as groups vector
         groups = factor(c(rep(paste(control, "RPF"), 3), rep(paste(treatment, "RPF"), 3),
                           rep(paste(control, "RNA"), 3), rep(paste(treatment, "RNA"), 3)), 
                         levels = c(paste(control, "RPF"), paste(treatment, "RPF"),
                                    paste(control, "RNA"), paste(treatment, "RNA")), ordered = T),
         anno = gene_anno,
         path = file.path(parent_dir, "plots"),
         folder = "Interactive_scatters",
         html = paste(treatment, "TE", sep = "_"))

