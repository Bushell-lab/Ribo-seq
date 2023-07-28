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

collapse_pathways <- function(fgsea_results, named_vector, fgsea_pathway, padj_threshold) {
  collapsedPathways <- collapsePathways(fgsea_results[order(pval)][padj < padj_threshold], 
                                        fgsea_pathway, named_vector)
  mainPathways <- fgsea_results[pathway %in% collapsedPathways$mainPathways][
    order(-NES), pathway]
  
  return(mainPathways)
}

extract_pathways <- function(fgsea_results, named_vectors, gsea_set, padj, groups = c("RPFs", "Totals", "TE")) {
  collapsed_pathways_list <- list()
  if ("RPFs" %in% groups) {
    collapsed_pathways_list$RPFs <- collapse_pathways(fgsea_results = fgsea_results[[1]],
                                                named_vector = named_vectors[[1]],
                                                gsea_set,
                                                padj_threshold = padj)
  }
  
  if ("Totals" %in% groups) {
    collapsed_pathways_list$Totals <- collapse_pathways(fgsea_results = fgsea_results[[2]],
                                                   named_vector = named_vectors[[2]],
                                                   gsea_set,
                                                   padj_threshold = padj)
  }
  
  if ("TE" %in% groups) {
    collapsed_pathways_list$TE <- collapse_pathways(fgsea_results = fgsea_results[[3]],
                                               named_vector = named_vectors[[3]],
                                               gsea_set,
                                               padj_threshold = padj)
  }
  
  return(collapsed_pathways_list)
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

plot_scatters <- function(df, gsea_set, pathway, dir) {
  gene_names <- gsea_set[[pathway]]
  
  df %>%
    mutate(group = factor(gene_sym %in% gene_names),
           alpha_score = case_when(group == T ~ 1,
                                   group == F ~ 0.1)) %>%
    arrange(group) %>%
    ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = group, alpha = alpha_score))+
    geom_point()+
    scale_colour_manual(values=c("grey", "red"))+
    scale_alpha(guide = "none")+
    geom_abline(lty = 2)+
    geom_hline(yintercept = 0, lty = 2)+
    geom_vline(xintercept = 0, lty = 2)+
    ylim(c(-2.5,2.5))+
    xlim(c(-2.5,2.5))+
    mytheme+
    xlab("Total RNA log2FC")+
    ylab("RPFs log2FC")+
    ggtitle(paste(treatment, str_replace_all(pathway, "_", " "), sep = "\n")) -> scatter_plot
  
  png(filename = file.path(parent_dir, "plots/fgsea/scatters", dir, paste(treatment, pathway, "TE_scatter_plot.png", sep = "_")), width = 500, height = 500)
  print(scatter_plot)
  dev.off()
}

make_interactive_scatter <- function(gsea_set, pathway, df, dir) {
  gene_names <- gsea_set[[pathway]]
  
  df %>%
    filter(gene_sym %in% gene_names) %>%
    mutate(groupings = factor(case_when(TE_log2FC < 0 & TE_padj < 0.1 ~ -1,
                                        TE_log2FC > 0 & TE_padj < 0.1 ~ 1,
                                        TE_padj >= 0.1 ~ 0))) -> filtered_data 
  
  
  #make gene ID/sym annotation table
  filtered_data %>%
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
  if (all(rownames(gene_anno) == filtered_data$gene) & all(rownames(gene_anno) == row.names(merged_norm_counts))) {
    #make html files
    glXYPlot(x = filtered_data$totals_log2FC,
             y = filtered_data$RPFs_log2FC,
             xlab = "RNA log2FC",
             ylab = "RPFs log2FC",
             main = paste(treatment, pathway),
             status = filtered_data$groupings,
             #cols = col_pal[1:3],
             counts = merged_norm_counts,
             side.xlab = "Sample",
             side.ylab = "norm counts",
             #sample.cols = rep(c("#F8766D", "#00BA38", "#619CFF"),4), #these are the colours for each sample within the norm counts, needs to be the same length as groups
             groups = factor(c(rep(paste(control, "RPF"), 3), rep(paste(treatment, "RPF"), 3),
                               rep(paste(control, "RNA"), 3), rep(paste(treatment, "RNA"), 3)), 
                             levels = c(paste(control, "RPF"), paste(treatment, "RPF"),
                                        paste(control, "RNA"), paste(treatment, "RNA")), ordered = T),
             anno = gene_anno,
             path = file.path(parent_dir, "plots/fgsea/Interactive_scatters"),
             folder = dir,
             html = paste(treatment, pathway, sep = "_"))
  }
}

#read in DESeq2 output----
DESeq2_data <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2.csv"))

#read in normalised counts
RPFs_norm_counts <- read.csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", treatment, "_normalised_counts.csv")))
totals_norm_counts <- read.csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", treatment, "_normalised_counts.csv")))

#make named vectors----
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

named_vectors <- list(RPFs_named_vector, totals_named_vector, TE_named_vector)

#read in pathways----
source("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/read_mouse_GSEA_pathways.R") #This may need to be changed to human

#hallmark----
#carry out fgsea
hallmark_results <- lapply(named_vectors, run_fgsea, pathway = pathways.hallmark)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/hallmark_results.Rdata"), hallmark_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_hallmark.png", sep = "_")), width = 500, height = 500)
make_plot(fgsea_result = hallmark_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Hallmark gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_hallmark.png", sep = "_")), width = 500, height = 500)
make_plot(fgsea_result = hallmark_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Hallmark gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_hallmark.png", sep = "_")), width = 500, height = 300)
make_plot(fgsea_result = hallmark_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Hallmark gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = hallmark_results, named_vectors = named_vectors, gsea_set = pathways.hallmark, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/scatters/hallmark")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/scatters/hallmark"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.hallmark, dir = "hallmark")

#plot interactive scatters
lapply(all_pathways$TE, make_interactive_scatter, gsea_set = pathways.hallmark, df = DESeq2_data, dir = "hallmark")

#biological processes----
#carry out fgsea
bio_processes_results <- lapply(named_vectors, run_fgsea, pathway = pathways.bio_processes)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/bio_processes_results.Rdata"), bio_processes_results)

#set adjusted p-value
padj <- 0.001

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_bio_processes.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = bio_processes_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Biological Processes gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_bio_processes.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = bio_processes_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Biological Processes gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_bio_processes.png", sep = "_")), width = 1500, height = 800)
make_plot(fgsea_result = bio_processes_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Biological Processes gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = bio_processes_results, named_vectors = named_vectors, gsea_set = pathways.bio_processes, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/scatters/bio_processes")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/scatters/bio_processes"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.bio_processes, dir = "bio_processes")

#plot interactive scatters
lapply(all_pathways$TE, make_interactive_scatter, gsea_set = pathways.bio_processes, df = DESeq2_data, dir = "bio_processes")

#molecular functions----
#carry out fgsea
mol_funs_results <- lapply(named_vectors, run_fgsea, pathway = pathways.mol_funs)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/mol_funs_results.Rdata"), mol_funs_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_mol_funs.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = mol_funs_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Molecular Functions gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_mol_funs.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = mol_funs_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Molecular Functions gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_mol_funs.png", sep = "_")), width = 2000, height = 800)
make_plot(fgsea_result = mol_funs_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Molecular Functions gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = mol_funs_results, named_vectors = named_vectors, gsea_set = pathways.mol_funs, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/scatters/mol_funs")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/scatters/mol_funs"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.mol_funs, dir = "mol_funs")

#plot interactive scatters
lapply(all_pathways$TE, make_interactive_scatter, gsea_set = pathways.mol_funs, df = DESeq2_data, dir = "mol_funs")

#cellular component----
#carry out fgsea
cell_comp_results <- lapply(named_vectors, run_fgsea, pathway = pathways.cell_comp)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/cell_comp_results.Rdata"), cell_comp_results)

#set adjusted p-value
padj <- 0.01

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_cell_comp.png", sep = "_")), width = 1000, height = 1200)
make_plot(fgsea_result = cell_comp_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Cellular Component gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_cell_comp.png", sep = "_")), width = 1000, height = 1200)
make_plot(fgsea_result = cell_comp_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Cellular Component gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_cell_comp.png", sep = "_")), width = 1000, height = 800)
make_plot(fgsea_result = cell_comp_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Cellular Component gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = cell_comp_results, named_vectors = named_vectors, gsea_set = pathways.cell_comp, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/scatters/cell_comp")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/scatters/cell_comp"))
}

lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.cell_comp, dir = "cell_comp")

#plot interactive scatters
lapply(all_pathways$TE, make_interactive_scatter, gsea_set = pathways.cell_comp, df = DESeq2_data, dir = "cell_comp")

#tumour----
#carry out fgsea
tumour_phen_onto_results <- lapply(named_vectors, run_fgsea, pathway = pathways.tumour_phen_onto)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/tumour_phen_onto_results.Rdata"), tumour_phen_onto_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_tumour_phen_onto.png", sep = "_")), width = 1000, height = 700)
make_plot(fgsea_result = tumour_phen_onto_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Tumour phenotype ontology gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_tumour_phen_onto.png", sep = "_")), width = 1000, height = 700)
make_plot(fgsea_result = tumour_phen_onto_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Tumour phenotype ontology gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_tumour_phen_onto.png", sep = "_")), width = 500, height = 200)
make_plot(fgsea_result = tumour_phen_onto_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Tumour phenotype ontology gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = tumour_phen_onto_results, named_vectors = named_vectors, gsea_set = pathways.tumour_phen_onto, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/scatters/tumour_phen_onto")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/scatters/tumour_phen_onto"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.tumour_phen_onto, dir = "tumour_phen_onto")

#plot interactive scatters
lapply(all_pathways$TE, make_interactive_scatter, gsea_set = pathways.tumour_phen_onto, df = DESeq2_data, dir = "tumour_phen_onto")

#Curated----
#read in pathways
#carry out fgsea
curated_results <- lapply(named_vectors, run_fgsea, pathway = pathways.curated)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/curated_results.Rdata"), curated_results)

#set adjusted p-value
padj <- 0.001

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_curated.png", sep = "_")), width = 1000, height = 2000)
make_plot(fgsea_result = curated_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Curated gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_curated.png", sep = "_")), width = 1000, height = 2000)
make_plot(fgsea_result = curated_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Curated gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_curated.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = curated_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Curated gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = curated_results, named_vectors = named_vectors, gsea_set = pathways.curated, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/scatters/curated")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/scatters/curated"))
}

lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.curated, dir = "curated")

#plot interactive scatters
lapply(all_pathways$TE, make_interactive_scatter, gsea_set = pathways.curated, df = DESeq2_data, dir = "curated")

#cell_type_sig----
#carry out fgsea
cell_type_sig_results <- lapply(named_vectors, run_fgsea, pathway = pathways.cell_type_sig)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/cell_type_sig_results.Rdata"), cell_type_sig_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_cell_type_sig.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = cell_type_sig_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Cell-type signature gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_cell_type_sig.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = cell_type_sig_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Cell-type signature gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_cell_type_sig.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = cell_type_sig_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Cell-type signature gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = cell_type_sig_results, named_vectors = named_vectors, gsea_set = pathways.cell_type_sig, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/scatters/cell_type_sig")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/scatters/cell_type_sig"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.cell_type_sig, dir = "cell_type_sig")

#plot interactive scatters
lapply(all_pathways$TE, make_interactive_scatter, gsea_set = pathways.cell_type_sig, df = DESeq2_data, dir = "cell_type_sig")

#transcription_factors----
#carry out fgsea
transcription_factors_results <- lapply(named_vectors, run_fgsea, pathway = pathways.transcription_factors)

#save results
save(file = file.path(parent_dir, "Analysis/fgsea/transcription_factors_results.Rdata"), transcription_factors_results)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_transcription_factors.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = transcription_factors_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Transcription factors gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_transcription_factors.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = transcription_factors_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Transcription factors gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_transcription_factors.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = transcription_factors_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Transcription factors gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = transcription_factors_results, named_vectors = named_vectors, gsea_set = pathways.transcription_factors, padj = padj, groups = "TE")

#plot overlaid scatters
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/scatters/transcription_factors")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/scatters/transcription_factors"))
}
lapply(all_pathways$TE, plot_scatters, df = DESeq2_data, gsea_set = pathways.transcription_factors, dir = "transcription_factors")

#plot interactive scatters
lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.transcription_factors, df = DESeq2_data, dir = "transcription_factors")

