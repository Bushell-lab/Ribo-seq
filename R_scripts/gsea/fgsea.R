#This is written for mouse data, will need to read in human pathways and edit pathway names if to be run on human data

#load libraries----
library(tidyverse)
library(fgsea)
library(Glimma)

#read in common variables----
source("common_variables.R")

#create a variable for what the treatment is----
control <- "WT"
treatment <- "KO"

#themes----
mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.position = "none")

#functions----
run_fgsea <- function(named_vector, pathway) {
  results <- fgsea(pathways = pathway,
                   stats=named_vector,
                   minSize = 20,
                   maxSize = 1000)
  return(results)
}

make_plot <- function(fgsea_result, padj_threshold, title) {
  plot <- ggplot(data = fgsea_result[fgsea_result$padj < padj_threshold], aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = padj)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=title) + 
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5))
  return(plot)
}

#read in DESeq2 output----
DESeq2_data <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2.csv"))

#make named vectors----
### This will make a list of named vectors to allow gsea to be carried out separately on the RPFs, Totals and TE log2FC
DESeq2_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(RPFs_log2FC)) %>%
  deframe() -> RPFs_named_vector

DESeq2_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(totals_log2FC)) %>%
  deframe() -> totals_named_vector

DESeq2_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(TE_log2FC)) %>%
  deframe() -> TE_named_vector

named_vectors <- list(RPFs = RPFs_named_vector, totals = totals_named_vector, TE = TE_named_vector)

#read in pathways----
source("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/read_mouse_GSEA_pathways.R") #This may need to be changed to human

#hallmark----
#carry out fgsea
hallmark_results <- lapply(named_vectors, run_fgsea, pathway = pathways.hallmark)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/hallmark_results.Rdata"), hallmark_results)

padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_hallmark.png", sep = "_")), width = 500, height = 500)
make_plot(fgsea_result = hallmark_results$RPFs, padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Hallmark gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_hallmark.png", sep = "_")), width = 500, height = 500)
make_plot(fgsea_result = hallmark_results$totals, padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Hallmark gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_hallmark.png", sep = "_")), width = 500, height = 300)
make_plot(fgsea_result = hallmark_results$TE, padj_threshold = padj, title = paste(treatment, "TE\nGSEA Hallmark gene sets"))
dev.off()

#biological processes----
#carry out fgsea
bio_processes_results <- lapply(named_vectors, run_fgsea, pathway = pathways.bio_processes)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/bio_processes_results.Rdata"), bio_processes_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_bio_processes.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = bio_processes_results$RPFs, padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Biological Processes gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_bio_processes.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = bio_processes_results$totals, padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Biological Processes gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_bio_processes.png", sep = "_")), width = 1500, height = 800)
make_plot(fgsea_result = bio_processes_results$TE, padj_threshold = padj, title = paste(treatment, "TE\nGSEA Biological Processes gene sets"))
dev.off()

#molecular functions----
#carry out fgsea
mol_funs_results <- lapply(named_vectors, run_fgsea, pathway = pathways.mol_funs)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/mol_funs_results.Rdata"), mol_funs_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_mol_funs.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = mol_funs_results$RPFs, padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Molecular Functions gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_mol_funs.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = mol_funs_results$totals, padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Molecular Functions gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_mol_funs.png", sep = "_")), width = 2000, height = 800)
make_plot(fgsea_result = mol_funs_results$TE, padj_threshold = padj, title = paste(treatment, "TE\nGSEA Molecular Functions gene sets"))
dev.off()

#cellular component----
#carry out fgsea
cell_comp_results <- lapply(named_vectors, run_fgsea, pathway = pathways.cell_comp)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/cell_comp_results.Rdata"), cell_comp_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_cell_comp.png", sep = "_")), width = 1000, height = 1200)
make_plot(fgsea_result = cell_comp_results$RPFs, padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Cellular Component gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_cell_comp.png", sep = "_")), width = 1000, height = 1200)
make_plot(fgsea_result = cell_comp_results$totals, padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Cellular Component gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_cell_comp.png", sep = "_")), width = 1000, height = 800)
make_plot(fgsea_result = cell_comp_results$TE, padj_threshold = padj, title = paste(treatment, "TE\nGSEA Cellular Component gene sets"))
dev.off()

#Curated----
### the curated is specific to the mouse GSEA lists
#read in pathways
#carry out fgsea
curated_results <- lapply(named_vectors, run_fgsea, pathway = pathways.curated)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/curated_results.Rdata"), curated_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_curated.png", sep = "_")), width = 1000, height = 2000)
make_plot(fgsea_result = curated_results$RPFs, padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Curated gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_curated.png", sep = "_")), width = 1000, height = 2000)
make_plot(fgsea_result = curated_results$totals, padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Curated gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_curated.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = curated_results$TE, padj_threshold = padj, title = paste(treatment, "TE\nGSEA Curated gene sets"))
dev.off()
