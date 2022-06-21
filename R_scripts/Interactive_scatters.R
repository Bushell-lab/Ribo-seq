#load libraries
library(tidyverse)
library(Glimma)

#read in common variables
source("common_variables.R")

#make variables----
treatment <- "EFT226"

#read in DESeq2 output----
RPFs <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", treatment, "_DEseq2_apeglm_LFC_shrinkage.csv")))
totals <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", treatment, "_DEseq2_apeglm_LFC_shrinkage.csv")))
RPFs_norm_counts <- read.csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", treatment, "_normalised_counts.csv")))
totals_norm_counts <- read.csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", treatment, "_normalised_counts.csv")))

#select just the Ctrl and treatment normalised counts----
RPFs_norm_counts %>%
  column_to_rownames("gene") %>%
  select(matches(c("Ctrl", treatment))) %>%
  rownames_to_column("GeneID") -> filtered_RPFs_norm_counts

totals_norm_counts %>%
  column_to_rownames("gene") %>%
  select(matches(c("Ctrl", treatment))) %>%
  rownames_to_column("GeneID") -> filtered_totals_norm_counts

#merge data----
RPFs %>%
  dplyr::select(transcript, gene, gene_sym, log2FoldChange, padj) %>%
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) %>%
  inner_join(totals[,c("transcript", "log2FoldChange", "padj")], by = "transcript") %>%
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj) %>%
  mutate(group = factor(case_when(RPFs_padj < 0.1 & RPFs_log2FC < 0 ~ "RPFs down",
                                  RPFs_padj < 0.1 & RPFs_log2FC > 0 ~ "RPFs up")),
         groupings = case_when(group == "RPFs down" ~ -1,
                               group == "RPFs up" ~ 1,
                               is.na(group)  ~ 0 )) -> merged_data
  
summary(merged_data)

#make gene ID/sym annotation table
merged_data %>%
  select(gene, gene_sym, transcript) %>%
  dplyr::rename(GeneID = gene) %>%
  as.data.frame() -> gene_anno
rownames(gene_anno) <- gene_anno$GeneID

#merge norm counts with gene annotation so that row names are in the same order and then select all counts columns in preferred order and make geneIDs row names
gene_anno %>%
  inner_join(filtered_RPFs_norm_counts, by = "GeneID") %>%
  inner_join(filtered_totals_norm_counts, by = "GeneID") %>%
  column_to_rownames("GeneID") %>%
  select(-c("transcript", "gene_sym")) -> merged_norm_counts

#check row names match up
all(rownames(gene_anno) == merged_data$gene)
all(rownames(gene_anno) == row.names(merged_norm_counts))

#make html files
glXYPlot(x = merged_data$totals_log2FC,
         y = merged_data$RPFs_log2FC,
         xlab = "RNA log2FC",
         ylab = "RPFs log2FC",
         main = treatment,
         status = merged_data$groupings,
         counts = merged_norm_counts,
         side.xlab = "Sample",
         side.ylab = "norm counts",
         sample.cols = rep(c("#F8766D", "#00BA38", "#619CFF"),4), #these are the colours for each sample within the norm counts, needs to be the same length as groups
         groups = factor(c(rep("Ctrl RPF", 3), rep(paste(treatment, "RPF"), 3),
                           rep("Ctrl RNA", 3), rep(paste(treatment, "RNA"), 3)), 
                         levels = c("Ctrl RPF",paste(treatment, "RPF"), "Ctrl RNA", paste(treatment, "RNA")), ordered = T),
         anno = gene_anno,
         path = file.path(parent_dir, "plots"),
         folder = "Interactive_scatters",
         html = paste(treatment, "TE", sep = "_"))

