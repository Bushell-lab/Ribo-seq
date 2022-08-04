#load libraries----
library(tidyverse)
library(fgsea)
library(Glimma)

#read in common variables----
source("common_variables.R")

#variables----
treatment <- "EFT226"

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

extract_pathways <- function(fgsea_results, named_vectors, gsea_set, padj) {
  RPF_collapsed_pathways <- collapse_pathways(fgsea_results = fgsea_results[[1]],
                                              named_vector = named_vectors[[1]],
                                              gsea_set,
                                              padj_threshold = padj)
  
  totals_collapsed_pathways <- collapse_pathways(fgsea_results = fgsea_results[[2]],
                                                 named_vector = named_vectors[[2]],
                                                 gsea_set,
                                                 padj_threshold = padj)
  
  TE_collapsed_pathways <- collapse_pathways(fgsea_results = fgsea_results[[3]],
                                             named_vector = named_vectors[[3]],
                                             gsea_set,
                                             padj_threshold = padj)
  
  
  all_pathways <- unique(c(RPF_collapsed_pathways, totals_collapsed_pathways, TE_collapsed_pathways))
  
  return(all_pathways)
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
    mutate(group = factor(Human_gene_name %in% gene_names),
           alpha_score = case_when(group == T ~ 1,
                                   group == F ~ 0.1)) %>%
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
    filter(Human_gene_name %in% gene_names) %>%
    mutate(groupings = factor(case_when(RPFs_log2FC < 0 & RPFs_padj < 0.1 ~ -1,
                                        RPFs_log2FC > 0 & RPFs_padj < 0.1 ~ 1,
                                        RPFs_padj >= 0.1 ~ 0))) -> filtered_data 
  
  
  #make gene ID/sym annotation table
  filtered_data %>%
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
             sample.cols = rep(c("#F8766D", "#00BA38", "#619CFF"),4), #these are the colours for each sample within the norm counts, needs to be the same length as groups
             groups = factor(c(rep("Ctrl RPF", 3), rep(paste(treatment, "RPF"), 3),
                               rep("Ctrl RNA", 3), rep(paste(treatment, "RNA"), 3)), 
                             levels = c("Ctrl RPF",paste(treatment, "RPF"), "Ctrl RNA", paste(treatment, "RNA")), ordered = T),
             anno = gene_anno,
             path = file.path(parent_dir, "plots/fgsea/Interactive_scatters"),
             folder = dir,
             html = paste(treatment, pathway, sep = "_"))
  }
}

#make a mouse to human gene conversion table----
#this is only applicable for mouse data as the gsea gene names are all human
read_tsv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/useful_tables/mouse_to_human_gene_IDs.tsv") %>%
  dplyr::select(Gene_stable_ID_version, Human_gene_name) %>%
  filter(!(is.na(Human_gene_name))) %>%
  group_by(Gene_stable_ID_version) %>%
  sample_n(size = 1) %>%
  dplyr::rename(gene = Gene_stable_ID_version) -> Mouse2HumanTable

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
  mutate(TE = RPFs_log2FC - totals_log2FC) %>%
  inner_join(Mouse2HumanTable, by = "gene") %>%
  filter(!(is.na(RPFs_padj)) & !(is.na(totals_padj))) -> merged_data

#make named vectors----
merged_data %>%
  group_by(Human_gene_name) %>%
  summarise(stat = mean(RPFs_log2FC)) %>%
  deframe() -> RPFs_named_vector

merged_data %>%
  group_by(Human_gene_name) %>%
  summarise(stat = mean(totals_log2FC)) %>%
  deframe() -> totals_named_vector

merged_data %>%
  group_by(Human_gene_name) %>%
  summarise(stat = mean(TE)) %>%
  deframe() -> TE_named_vector

named_vectors <- list(RPFs_named_vector, totals_named_vector, TE_named_vector)

#hallmark----
#read in pathway
pathways.hallmark <- gmtPathways("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/h.all.v7.1.symbols.gmt")

#carry out fgsea
hallmark_results <- lapply(named_vectors, run_fgsea, pathway = pathways.hallmark)

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
all_pathways <- extract_pathways(fgsea_results = hallmark_results, named_vectors = named_vectors, gsea_set = pathways.hallmark, padj = padj)

#plot overlaid scatters
dir.create(file.path(parent_dir, "plots/fgsea/scatters/hallmark"))
lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.hallmark, dir = "hallmark")

#plot interactive scatters
lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.hallmark, df = merged_data, dir = "hallmark")

#biological processes----
#read in pathway
pathways.bio_processes <- gmtPathways("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/c5.bp.v7.1.symbols.gmt")

#carry out fgsea
bio_processes_results <- lapply(named_vectors, run_fgsea, pathway = pathways.bio_processes)

#set adjusted p-value
padj <- 0.001

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_bio_processes.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = bio_processes_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Biological Processes gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_bio_processes.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = bio_processes_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Biological Processes gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_bio_processes.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = bio_processes_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Biological Processes gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = bio_processes_results, named_vectors = named_vectors, gsea_set = pathways.bio_processes, padj = padj)

#plot overlaid scatters
dir.create(file.path(parent_dir, "plots/fgsea/scatters/bio_processes"))
lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.bio_processes, dir = "bio_processes")

#plot interactive scatters
lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.bio_processes, df = merged_data, dir = "bio_processes")

#molecular functions----
#read in pathway
pathways.mol_funs <- gmtPathways("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/c5.mf.v7.1.symbols.gmt")

#carry out fgsea
mol_funs_results <- lapply(named_vectors, run_fgsea, pathway = pathways.mol_funs)

#set adjusted p-value
padj <- 0.001

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_mol_funs.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = mol_funs_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Molecular Functions gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_mol_funs.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = mol_funs_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Molecular Functions gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_mol_funs.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = mol_funs_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Molecular Functions gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = mol_funs_results, named_vectors = named_vectors, gsea_set = pathways.mol_funs, padj = padj)

#plot overlaid scatters
dir.create(file.path(parent_dir, "plots/fgsea/scatters/mol_funs"))
lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.mol_funs, dir = "mol_funs")

#plot interactive scatters
lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.mol_funs, df = merged_data, dir = "mol_funs")

#cellular component----
#read in pathway
pathways.cell_comp <- gmtPathways("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/c5.cc.v7.1.symbols.gmt")

#carry out fgsea
cell_comp_results <- lapply(named_vectors, run_fgsea, pathway = pathways.cell_comp)

#set adjusted p-value
padj <- 0.001

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_cell_comp.png", sep = "_")), width = 1000, height = 1200)
make_plot(fgsea_result = cell_comp_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Cellular Component gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_cell_comp.png", sep = "_")), width = 1000, height = 1200)
make_plot(fgsea_result = cell_comp_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Cellular Component gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_cell_comp.png", sep = "_")), width = 1000, height = 1200)
make_plot(fgsea_result = cell_comp_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Cellular Component gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = cell_comp_results, named_vectors = named_vectors, gsea_set = pathways.cell_comp, padj = padj)

#plot overlaid scatters
dir.create(file.path(parent_dir, "plots/fgsea/scatters/cell_comp"))
lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.cell_comp, dir = "cell_comp")

#plot interactive scatters
lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.cell_comp, df = merged_data, dir = "cell_comp")

#onco----
#read in pathways
pathways.onco <- gmtPathways("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/c6.all.v7.1.symbols.gmt")

#carry out fgsea
onco_results <- lapply(named_vectors, run_fgsea, pathway = pathways.onco)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_onco.png", sep = "_")), width = 1000, height = 700)
make_plot(fgsea_result = onco_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Oncogenic Signature gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_onco.png", sep = "_")), width = 1000, height = 700)
make_plot(fgsea_result = onco_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Oncogenic Signature gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_onco.png", sep = "_")), width = 1000, height = 500)
make_plot(fgsea_result = onco_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Oncogenic Signature gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = onco_results, named_vectors = named_vectors, gsea_set = pathways.onco, padj = padj)

#plot overlaid scatters
dir.create(file.path(parent_dir, "plots/fgsea/scatters/onco"))
lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.onco, dir = "onco")

#plot interactive scatters
lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.onco, df = merged_data, dir = "onco")

#immunologic signature----
#read in pathways
pathways.immune <- gmtPathways("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/c7.all.v7.1.symbols.gmt")

#carry out fgsea
immune_results <- lapply(named_vectors, run_fgsea, pathway = pathways.immune)


#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_immune.png", sep = "_")), width = 1000, height = 2000)
make_plot(fgsea_result = immune_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Immunologic Signature gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_immune.png", sep = "_")), width = 1000, height = 2000)
make_plot(fgsea_result = immune_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Immunologic Signature gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_immune.png", sep = "_")), width = 1000, height = 1000)
make_plot(fgsea_result = immune_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Immunologic Signature gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = immune_results, named_vectors = named_vectors, gsea_set = pathways.immune, padj = padj)

#plot overlaid scatters
dir.create(file.path(parent_dir, "plots/fgsea/scatters/immune"))
lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.immune, dir = "immune")

#plot interactive scatters
lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.immune, df = merged_data, dir = "immune")

#kegg----
#read in pathways
pathways.kegg <- gmtPathways("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/c2.cp.kegg.v7.1.symbols.gmt")

#carry out fgsea
kegg_results <- lapply(named_vectors, run_fgsea, pathway = pathways.kegg)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "RPFs_kegg.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = kegg_results[[1]], padj_threshold = padj, title = paste(treatment, "RPFs\nGSEA Kegg gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "totals_kegg.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = kegg_results[[2]], padj_threshold = padj, title = paste(treatment, "Total RNA\nGSEA Kegg gene sets"))
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/", paste(treatment, "TE_kegg.png", sep = "_")), width = 700, height = 700)
make_plot(fgsea_result = kegg_results[[3]], padj_threshold = padj, title = paste(treatment, "TE\nGSEA Kegg gene sets"))
dev.off()

#extract pathways
all_pathways <- extract_pathways(fgsea_results = kegg_results, named_vectors = named_vectors, gsea_set = pathways.kegg, padj = padj)

#plot overlaid scatters
dir.create(file.path(parent_dir, "plots/fgsea/scatters/kegg"))
lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.kegg, dir = "kegg")

#plot interactive scatters
lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.kegg, df = merged_data, dir = "kegg")

