#load packages----
library(tidyverse)
library(grid)
library(gridExtra)

#read in common variables
source("common_variables.R")

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
load(file = file.path(parent_dir, "Counts_files/R_objects/counts_list.Rdata"))

most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "Analysis/transcript_IDs/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))
region_lengths <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_region_lengths.csv", col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))

#run on an individual genes----

#dir is sub directory (within "plots/binned_plots/single_transcripts" where plots will be saved. It will create directory if it does not already exist)
#plot_binned/plot_single_nt sets whether to plot either or both of these
#SD=T will add standard deviation bars to line plots, set to false if not wanted
#plot_replicates=T will create a seperate binned plot with individual replicates rather than average
#control and treatment need to state what you have called your control and treatment samples. You can specificy more than one treatment as a vector, but you will have to set plot_delta=F
#keep paired_data=T if data is paired, otherwise set to F. This is for calculating 95% confidence intervals for the delta plots
#region_cutoffs = c(0,0,0) will mean any transcript length can be plotted. For meta plots, this is normally set to region_cutoffs = c(50,300,50)

plot_single_transcripts(gene = "Ctnnb1", dir = "candidate",
                        plot_binned = T, plot_single_nt = T,
                        SD = T, plot_replicates = T, plot_delta = T,
                        control = "Ctrl", treatment = "EFT226", paired_data = T,
                        region_cutoffs = c(0,0,0))

#read in DEseq2 output----
read_csv(file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2.csv")) %>%
  inner_join(most_abundant_transcripts, by = "gene") %>%
  mutate(RPFs_group = factor(RPFs_group),
         TE_group = factor(TE_group)) -> DESeq2_df
summary(DESeq2_df)

#extract IDs
DESeq2_df %>%
  filter(RPFs_group == "RPFs down") %>%
  top_n(n = -100, wt = RPFs_padj) %>%
  top_n(n = -50, wt = RPFs_log2FC) %>%
  pull(gene_sym) -> RPFs_down_IDs

DESeq2_df %>%
  filter(RPFs_group == "RPFs up") %>%
  top_n(n = -100, wt = RPFs_padj) %>%
  top_n(n = 50, wt = RPFs_log2FC) %>%
  pull(gene_sym) -> RPFs_up_IDs

#run on gene lists
lapply(RPFs_down_IDs, plot_single_transcripts, dir = "RPFs_down",
       plot_binned = T, plot_single_nt = F,
       SD = T, plot_replicates = F, plot_delta = F,
       control = "Ctrl", treatment = "EFT226", paired_data = T,
       region_cutoffs = c(0,0,0))

lapply(RPFs_up_IDs, plot_single_transcripts, dir = "RPFs_up",
       plot_binned = T, plot_single_nt = F,
       SD = T, plot_replicates = F, plot_delta = F,
       control = "Ctrl", treatment = "EFT226", paired_data = T,
       region_cutoffs = c(0,0,0))
