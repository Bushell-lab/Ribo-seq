#load packages----
library(tidyverse)
library(grid)
library(gridExtra)
library(viridis)

#read in common variables----
source("\\\\data.beatson.gla.ac.uk/data/JWALDRON/Ribosome_profiling/Organoids/4AIn/Scripts/R_scripts/common_variables.R")
control <- "Ctrl"
treatment <- "EFT226"

#read in functions----
source("\\\\data.beatson.gla.ac.uk/data/JWALDRON/Scripts/R/Functions/binning_RiboSeq_functions.R")

#create themes----
my_theme <- theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_blank())

UTR5_theme <- my_theme+
  theme(legend.position="none",
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank())

CDS_theme <- my_theme+
  theme(legend.position="none",
        axis.ticks = element_blank(),
        axis.text = element_blank())

UTR3_theme <- my_theme+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_blank())

#read in data----
load(file = file.path(parent_dir, "Counts_files/R_objects/binned_list.Rdata"))
summary(binned_list[[1]])
head(binned_list[[1]])

load(file = file.path(parent_dir, "Counts_files/R_objects/single_nt_list.Rdata"))
summary(single_nt_list[[1]])
head(single_nt_list[[1]])

region_lengths <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_region_lengths.csv", col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))

#set whether you want to plot binned_cpm (binned RPF counts per million (normalised by read depth only)) or binned_normalised_cpm (binned RPF CPMs normalised to total RNA TPMs)
#hash out which ever you do not want to plot
binned_value <- "binned_cpm"
#binned_value <- "binned_normalised_cpm"

#all transcripts----
#summarise
#summarise within each sample
summarised_binned_list <- lapply(binned_list, summarise_data, value = binned_value, grouping = "bin")
summary(summarised_binned_list[[1]])
print(summarised_binned_list[[1]])

#summarise within each condition (across replicates)
do.call("rbind", summarised_binned_list) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_binned
summary(summarised_binned)
print(summarised_binned)

#plot lines
binned_line_plots <- plot_binned_lines(df = summarised_binned, SD = T, control = control, treatment = treatment)

png(filename = file.path(parent_dir, "plots/binned_plots", paste(treatment, "all transcripts", binned_value, "lines.png")), width = 1000, height = 300)
grid.arrange(binned_line_plots[[1]], binned_line_plots[[2]], binned_line_plots[[3]], nrow = 1, widths = c(1,2,1.5))
dev.off()

binned_line_plots_all_replicates <- plot_binned_all_replicates(summarised_binned_list, control = control, treatment = treatment)

png(filename = file.path(parent_dir, "plots/binned_plots", paste(treatment, "all transcripts", binned_value, "lines all replicates.png")), width = 1000, height = 300)
grid.arrange(binned_line_plots_all_replicates[[1]], binned_line_plots_all_replicates[[2]], binned_line_plots_all_replicates[[3]], nrow = 1, widths = c(1,2,1.5))
dev.off()

#calculate and plot delta
binned_delta_data <- calculate_binned_delta(binned_list, value = binned_value, control = control, treatment = treatment, paired_data = T)
binned_delta_plots <- plot_binned_delta(binned_delta_data)

png(filename = file.path(parent_dir, "plots/binned_plots/", paste(treatment, "all transcripts", binned_value, "delta.png")), width = 1000, height = 200)
grid.arrange(binned_delta_plots[[1]], binned_delta_plots[[2]], binned_delta_plots[[3]],
             nrow = 1, widths = c(1,2,1))
dev.off()

#positional----
#normalise within each transcript
positional_list <- lapply(binned_list, calculate_positional_counts)

#summarise within each sample
summarised_positional_list <- lapply(positional_list, summarise_data, value = "positional_counts", grouping = "bin")

#summarise within each condition (across replicates)
do.call("rbind", summarised_positional_list) %>%
  group_by(grouping, condition) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_positional

#plot
positional_line_plots <- plot_positional_lines(df = summarised_positional, SD = T, control = control, treatment = treatment)

png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "_all_transcripts_binned_positional_lines.png")), width = 500, height = 200)
print(positional_line_plots)
dev.off()

#calculate and plot delta
binned_positional_delta <- calculate_positional_delta(positional_list, control = control, treatment = treatment, paired_data = F)
positional_binned_delta_plots <- plot_positional_delta(binned_positional_delta)

png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "_all_transcripts_binned_positional_delta.png")), width = 500, height = 200)
print(positional_binned_delta_plots)
dev.off()

#single nt----
#set whether you want to plot binned_cpm (binned RPF counts per million (normalised by read depth only)) or binned_normalised_cpm (binned RPF CPMs normalised to total RNA TPMs)
#hash out which ever you do not want to plot
single_nt_value <- "single_nt_cpm"
#single_nt_value <- "single_nt_normalised_cpm"

#summarise
summarised_single_nt_list <- lapply(single_nt_list, summarise_data, value = single_nt_value, grouping = "window")
summary(summarised_single_nt_list[[1]])
print(summarised_single_nt_list[[1]])

do.call("rbind", summarised_single_nt_list) %>%
  group_by(grouping, condition, region) %>%
  summarise(average_counts = mean(mean_counts),
            sd_counts = sd(mean_counts)) %>%
  ungroup() -> summarised_single_nt
summary(summarised_single_nt)
print(summarised_single_nt)

#plot
single_nt_line_plots <- plot_single_nt_lines(summarised_single_nt, SD=T, plot_ends=F, control = control, treatment = treatment)

png(filename = file.path(parent_dir, "plots/binned_plots", paste0(treatment, "_all_transcripts_", single_nt_value, "s_lines.png")), width = 1300, height = 300)
grid.arrange(single_nt_line_plots[[1]], single_nt_line_plots[[2]], single_nt_line_plots[[3]], single_nt_line_plots[[4]], nrow = 1, widths = c(1,2,2,1.5))
dev.off()

#calculate and plot delta
single_nt_delta_data <- calculate_single_nt_delta(single_nt_list, value = single_nt_value, control = control, treatment = treatment)
single_nt_delta_plots <- plot_single_nt_delta(single_nt_delta_data, SD = T)

png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "_all_transcripts_", single_nt_value, "s_delta.png")), width = 1300, height = 200)
grid.arrange(single_nt_delta_plots[[1]], single_nt_delta_plots[[2]], single_nt_delta_plots[[3]], single_nt_delta_plots[[4]], nrow = 1, widths = c(1,2,2,1))
dev.off()

#Dep vs Antidep vs Indep----
#set whether you want to plot binned_cpm (binned RPF counts per million (normalised by read depth only)) or binned_normalised_cpm (binned RPF CPMs normalised to total RNA TPMs)
#hash out which ever you do not want to plot

#binned_value <- "binned_cpm"
binned_value <- "binned_normalised_cpm"

#single_nt_value <- "single_nt_cpm"
single_nt_value <- "single_nt_normalised_cpm"

#read in DESeq2 output
RPFs <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("AK_RPFs_EFT226_DEseq2_apeglm_LFC_shrinkage.csv")))
totals <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("AK_Totals_EFT226_DEseq2_apeglm_LFC_shrinkage.csv")))

#merge RPFs with totals and assign groups
RPFs %>%
  dplyr::select(transcript, gene, gene_sym, log2FoldChange, padj) %>%
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) %>%
  inner_join(totals[,c("transcript", "log2FoldChange", "padj")], by = "transcript") %>%
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj) %>%
  mutate(TE = RPFs_log2FC - totals_log2FC,
         RPF_group = factor(case_when(RPFs_padj < 0.1 & RPFs_log2FC < 0 & totals_padj >= 0.1 ~ "RPFs down",
                                  RPFs_padj < 0.1 & RPFs_log2FC > 0 & totals_padj >= 0.1 ~ "RPFs up",
                                  RPFs_padj >= 0.1 & totals_padj < 0.1 & totals_log2FC < 0 ~"Totals down",
                                  RPFs_padj >= 0.1 & totals_padj < 0.1 & totals_log2FC > 0 ~"Totals up",
                                  RPFs_padj < 0.1 & totals_padj < 0.1 & RPFs_log2FC < 0 & totals_log2FC < 0 ~ "both down",
                                  RPFs_padj < 0.1 & totals_padj < 0.1 & RPFs_log2FC > 0 & totals_log2FC > 0 ~ "both up",
                                  RPFs_padj > 0.1 & totals_padj > 0.1 ~ "no change")),
         TE_group = factor(case_when(TE < -0.2 ~ "TE down",
                                     TE > 0.2 ~ "TE up",
                                     TE > -0.1 & TE < 0.1 ~ "No TE change"))) -> merged_RPF_data
summary(merged_RPF_data)

merged_RPF_data %>%
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = RPF_group))+
  geom_point()

#extract transcript IDs
down_IDs <- merged_RPF_data$transcript[merged_RPF_data$RPF_group == "RPFs down" & !(is.na(merged_RPF_data$RPF_group))]
up_IDs <- merged_RPF_data$transcript[merged_RPF_data$RPF_group == "RPFs up" & !(is.na(merged_RPF_data$RPF_group))]
no_change_IDs <- merged_RPF_data$transcript[merged_RPF_data$RPF_group == "no change" & !(is.na(merged_RPF_data$RPF_group))]

#dep
dep_plots <- plot_subset(IDs = down_IDs, binned_value = binned_value, single_nt_value = single_nt_value)

#binned lines
png(filename = file.path(parent_dir, "plots/binned_plots", paste0(treatment, "_dep_", binned_value, "s_lines.png")), width = 1000, height = 200)
grid.arrange(dep_plots[[1]][[1]], dep_plots[[1]][[2]], dep_plots[[1]][[3]], nrow = 1, widths = c(1,2,1.5))
dev.off()

#binned delta
png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "_dep_", binned_value, "s_delta.png")), width = 1000, height = 200)
grid.arrange(dep_plots[[2]][[1]], dep_plots[[2]][[2]], dep_plots[[2]][[3]], nrow = 1, widths = c(1,2,1))
dev.off()

#single_nt lines
png(filename = file.path(parent_dir, "plots/binned_plots", paste0(treatment, "_dep_", single_nt_value, "s_lines.png")), width = 1300, height = 300)
grid.arrange(dep_plots[[3]][[1]], dep_plots[[3]][[2]], dep_plots[[3]][[3]], dep_plots[[3]][[4]], nrow = 1, widths = c(1,2,2,1.5))
dev.off()

#single_nt delta
png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "_dep_", single_nt_value, "s_delta.png")), width = 1300, height = 200)
grid.arrange(dep_plots[[4]][[1]], dep_plots[[4]][[2]], dep_plots[[4]][[3]], dep_plots[[4]][[4]], nrow = 1, widths = c(1,2,2,1))
dev.off()

#positional lines
png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "dep_binned_positional_lines.png")), width = 500, height = 200)
print(dep_plots[[5]])
dev.off()

#positional delta
png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "dep_binned_positional_delta.png")), width = 500, height = 200)
print(dep_plots[[6]])
dev.off()

#anti
anti_plots <- plot_subset(IDs = up_IDs, binned_value = binned_value, single_nt_value = single_nt_value)

#binned lines
png(filename = file.path(parent_dir, "plots/binned_plots", paste0(treatment, "_anti_", binned_value, "s_lines.png")), width = 1000, height = 200)
grid.arrange(anti_plots[[1]][[1]], anti_plots[[1]][[2]], anti_plots[[1]][[3]], nrow = 1, widths = c(1,2,1.5))
dev.off()

#binned delta
png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "_anti_", binned_value, "s_delta.png")), width = 1000, height = 200)
grid.arrange(anti_plots[[2]][[1]], anti_plots[[2]][[2]], anti_plots[[2]][[3]], nrow = 1, widths = c(1,2,1))
dev.off()

#single_nt lines
png(filename = file.path(parent_dir, "plots/binned_plots", paste0(treatment, "_anti_", single_nt_value, "s_lines.png")), width = 1000, height = 200)
grid.arrange(anti_plots[[3]][[1]], anti_plots[[3]][[2]], anti_plots[[3]][[3]], anti_plots[[3]][[4]], nrow = 1, widths = c(rep(1,3),1.2))
dev.off()

#single_nt delta
png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "_anti_", single_nt_value, "s_delta.png")), width = 1000, height = 200)
grid.arrange(anti_plots[[4]][[1]], anti_plots[[4]][[2]], anti_plots[[4]][[3]], anti_plots[[4]][[4]], nrow = 1)
dev.off()

#positional lines
png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "anti_binned_positional_lines.png")), width = 500, height = 200)
print(anti_plots[[5]])
dev.off()

#positional delta
png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "anti_binned_positional_delta.png")), width = 500, height = 200)
print(anti_plots[[6]])
dev.off()

#indep
indep_plots <- plot_subset(IDs = no_change_IDs, binned_value = binned_value, single_nt_value = single_nt_value)

#binned lines
png(filename = file.path(parent_dir, "plots/binned_plots", paste0(treatment, "_indep_", binned_value, "s_lines.png")), width = 1000, height = 200)
grid.arrange(indep_plots[[1]][[1]], indep_plots[[1]][[2]], indep_plots[[1]][[3]], nrow = 1, widths = c(1,2,1.5))
dev.off()

#binned delta
png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "_indep_", binned_value, "s_delta.png")), width = 1000, height = 200)
grid.arrange(indep_plots[[2]][[1]], indep_plots[[2]][[2]], indep_plots[[2]][[3]], nrow = 1, widths = c(1,2,1))
dev.off()

#single_nt lines
png(filename = file.path(parent_dir, "plots/binned_plots", paste0(treatment, "_indep_", single_nt_value, "s_lines.png")), width = 1000, height = 200)
grid.arrange(indep_plots[[3]][[1]], indep_plots[[3]][[2]], indep_plots[[3]][[3]], indep_plots[[3]][[4]], nrow = 1, widths = c(rep(1,3),1.2))
dev.off()

#single_nt delta
png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "_indep_", single_nt_value, "s_delta.png")), width = 1000, height = 200)
grid.arrange(indep_plots[[4]][[1]], indep_plots[[4]][[2]], indep_plots[[4]][[3]], indep_plots[[4]][[4]], nrow = 1)
dev.off()

#positional lines
png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "indep_binned_positional_lines.png")), width = 500, height = 200)
print(indep_plots[[5]])
dev.off()

#positional delta
png(filename = file.path(parent_dir, "plots/binned_plots/", paste0(treatment, "indep_binned_positional_delta.png")), width = 500, height = 200)
print(indep_plots[[6]])
dev.off()

#GSEA pathways----
library(fgsea)
#make a mouse to human gene conversion table (this is only applicable for mouse data as the gsea gene names are all human)
read_tsv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/useful_tables/mouse_to_human_gene_IDs.tsv") %>%
  dplyr::select(Gene_stable_ID_version, Human_gene_name) %>%
  filter(!(is.na(Human_gene_name))) %>%
  group_by(Gene_stable_ID_version) %>%
  sample_n(size = 1) %>%
  dplyr::rename(gene = Gene_stable_ID_version) -> Mouse2HumanTable

#read in most abundant transcripts
most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/AK_most_abundant_transcripts_IDs.csv"))

#kegg
#read in pathways
#pathways.cell_comp <- gmtPathways("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/c5.cc.v7.1.symbols.gmt")
pathways.kegg <- gmtPathways("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/c2.cp.kegg.v7.1.symbols.gmt")

#read in fgsea output
load(file = file.path(parent_dir, "Analysis/fgsea/kegg_results.Rdata"))

#extract transcript IDs
kegg_pathways <- kegg_results[[3]]$pathway[kegg_results[[3]]$padj < 0.05]

lapply(kegg_pathways, plot_GSEA_binned, GSEA_set = pathways.kegg, dir = "kegg", human = F, conversion_table = Mouse2HumanTable,
       binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
       plot_binned = T, plot_single_nt = F, plot_positional = F)



