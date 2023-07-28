#load packages----
library(tidyverse)
library(grid)
library(gridExtra)
library(viridis)

#read in common variables----
source("common_variables.R")

#set what you have called your control and treated samples. This can be a vector of strings if more than one treatment has been used.
control <- "WT"
treatment <- "KO"

#read in functions----
source("binning_RiboSeq_functions.R")

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
most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

#all transcripts----
#summarise
#summarise within each sample
summarised_binned_list <- lapply(binned_list, summarise_data, value = "binned_cpm", grouping = "bin")
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

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned lines.png"), width = 1000, height = 300)
grid.arrange(binned_line_plots[[1]], binned_line_plots[[2]], binned_line_plots[[3]], nrow = 1, widths = c(1,2,1.5))
dev.off()

binned_line_plots_all_replicates <- plot_binned_all_replicates(summarised_binned_list, control = control, treatment = treatment)

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned lines all replicates.png"), width = 1000, height = 300)
grid.arrange(binned_line_plots_all_replicates[[1]], binned_line_plots_all_replicates[[2]], binned_line_plots_all_replicates[[3]], nrow = 1, widths = c(1,2,1.5))
dev.off()

#calculate and plot delta
binned_delta_data <- calculate_binned_delta(binned_list, value = "binned_cpm", control = control, treatment = treatment, paired_data = F)
binned_delta_plots <- plot_binned_delta(binned_delta_data)

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned delta.png"), width = 1000, height = 200)
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

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned positional lines.png"), width = 500, height = 200)
print(positional_line_plots)
dev.off()

#calculate and plot delta
binned_positional_delta <- calculate_positional_delta(positional_list, control = control, treatment = treatment, paired_data = F)
positional_binned_delta_plots <- plot_positional_delta(binned_positional_delta)

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts binned positional delta.png"), width = 500, height = 200)
print(positional_binned_delta_plots)
dev.off()

#single nt----
#summarise
summarised_single_nt_list <- lapply(single_nt_list, summarise_data, value = "single_nt_cpm", grouping = "window")
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

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts single nt lines.png"), width = 1300, height = 300)
grid.arrange(single_nt_line_plots[[1]], single_nt_line_plots[[2]], single_nt_line_plots[[3]], single_nt_line_plots[[4]], nrow = 1, widths = c(1,2,2,1.5))
dev.off()

#calculate and plot delta
single_nt_delta_data <- calculate_single_nt_delta(single_nt_list, value = "single_nt_cpm", control = control, treatment = treatment, paired_data = F)
single_nt_delta_plots <- plot_single_nt_delta(single_nt_delta_data, SD = T)

png(filename = file.path(parent_dir, "plots/binned_plots/all_transcripts/all transcripts single nt delta.png"), width = 1300, height = 200)
grid.arrange(single_nt_delta_plots[[1]], single_nt_delta_plots[[2]], single_nt_delta_plots[[3]], single_nt_delta_plots[[4]], nrow = 1, widths = c(1,2,2,1))
dev.off()

#Dep vs Antidep vs Indep----
#read in DESeq2 output
read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2.csv")) %>%
  inner_join(most_abundant_transcripts, by = c("gene", "gene_sym")) %>%
  mutate(RPFs_group = factor(RPFs_group),
         TE_group = factor(TE_group)) -> DESeq2_data

summary(DESeq2_data)

#plot to check groupings
DESeq2_data %>%
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = RPFs_group))+
  geom_point()

#plot to check groupings
DESeq2_data %>%
  filter(!(is.na(TE_group))) %>%
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = TE_group))+
  geom_point()

#extract transcript IDs
RPFs_down_IDs <- DESeq2_data$transcript[DESeq2_data$RPFs_group == "RPFs down" & !(is.na(DESeq2_data$RPFs_group))]
RPFs_up_IDs <- DESeq2_data$transcript[DESeq2_data$RPFs_group == "RPFs up" & !(is.na(DESeq2_data$RPFs_group))]

TE_down_IDs <- DESeq2_data$transcript[DESeq2_data$TE_group == "TE down" & !(is.na(DESeq2_data$TE_group))]
TE_up_IDs <- DESeq2_data$transcript[DESeq2_data$TE_group == "TE up" & !(is.na(DESeq2_data$TE_group))]

no_change_IDs <- DESeq2_data$transcript[DESeq2_data$TE_group == "no change" & !(is.na(DESeq2_data$TE_group))]

#plot binned
plot_subset(IDs = RPFs_down_IDs, subset = "RPFs-down", sub_dir = "Dep",
            control = control, treatment = treatment,
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = F,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

plot_subset(IDs = RPFs_up_IDs, subset = "RPFs-up", sub_dir = "Dep",
            control = control, treatment = treatment,
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = F,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

plot_subset(IDs = TE_down_IDs, subset = "TE-down", sub_dir = "Dep",
            control = control, treatment = treatment,
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = F,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

plot_subset(IDs = TE_up_IDs, subset = "TE-up", sub_dir = "Dep",
            control = control, treatment = treatment,
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = F,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

plot_subset(IDs = no_change_IDs, subset = "no_change", sub_dir = "Dep",
            control = control, treatment = treatment,
            binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
            plot_binned = T, plot_single_nt = F, plot_positional = F,
            plot_replicates = T, plot_delta = T, SD = T, paired_data = F)

#plot heatmaps
TE_down_heatmap <- plot_binned_heatmaps(IDs = TE_down_IDs, col_lims = c(-0.02, 0.01), control = control, treatment = treatment, value = "binned_normalised_cpm")
TE_up_heatmap <- plot_binned_heatmaps(IDs = TE_up_IDs, col_lims = c(-0.01, 0.02), control = control, treatment = treatment, value = "binned_normalised_cpm")

png(filename = file.path(parent_dir, "plots/binned_plots/Dep", paste(treatment, "TE-down binned heatmap.png")), width = 1000, height = 1000)
print(TE_down_heatmap)
dev.off()

png(filename = file.path(parent_dir, "plots/binned_plots/Dep", paste(treatment, "TE-up binned heatmap.png")), width = 1000, height = 1000)
print(TE_up_heatmap)
dev.off()

#GSEA pathways----
library(fgsea)

#read in pathways
source("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/read_mouse_GSEA_pathways.R")

#hallmark
#read in fgsea output
load(file = file.path(parent_dir, "Analysis/fgsea/hallmark_results.Rdata"))

#extract transcript IDs
hallmark_pathways <- hallmark_results[[3]]$pathway[hallmark_results[[3]]$padj < 0.05]

lapply(hallmark_pathways, plot_GSEA_binned,
       GSEA_set = pathways.hallmark, sub_dir = "hallmark",
       control = control, treatment = treatment,
       binned_value = "binned_normalised_cpm", single_nt_value = "single_nt_normalised_cpm",
       plot_binned = T, plot_single_nt = F, plot_positional = F,
       plot_delta = T, SD = T, paired_data = F)



