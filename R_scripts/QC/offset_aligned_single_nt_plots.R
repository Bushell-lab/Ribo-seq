#load packages----
library(tidyverse)
library(grid)
library(gridExtra)

#read in common variables----
source("common_variables.R")

#read in data----
region_lengths <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_region_lengths.csv", col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))

#create themes----
my_theme <- theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 16),
        legend.position="none")

#Read in the counts csvs, calculate Counts per Million (CPM) and splice into first and last 25/50nt of CDS/UTRs----
spliced_list <- list()
for (sample in RPF_sample_names) {
  
  #read in counts csv
  df <- read_csv(file = file.path(parent_dir, "Counts_files/csv_files", paste0(sample, "_pc_final_counts.csv")))
  
  total_counts <- sum(df$Counts)
  
  df %>%
    mutate(CPM = (Counts / total_counts) * 1000000)  %>%
    inner_join(region_lengths, by = "transcript") -> merged_data
  
  #splice the transcript
  #UTR5 end
  merged_data %>%
    filter(Position <= UTR5_len) %>%
    group_by(transcript) %>%
    top_n(n = 25, wt = Position) %>% #extracts the 3' most nts
    ungroup() %>%
    mutate(nt = (Position - UTR5_len) - 1) %>%
    group_by(nt) %>%
    summarise(mean_cpm = mean(CPM)) %>%
    ungroup() %>%
    mutate(region = rep("UTR5")) -> UTR5_end
  
  #CDS start
  merged_data %>%
    filter(Position > UTR5_len & Position <= (UTR5_len + CDS_len)) %>%
    mutate(nt = Position - UTR5_len) %>%
    filter(nt <= 50) %>%
    group_by(nt) %>%
    summarise(mean_cpm = mean(CPM)) %>%
    ungroup() %>%
    mutate(region = rep("CDS")) -> CDS_start
  
  #CDS end
  merged_data %>%
    filter(Position > UTR5_len & Position <= (UTR5_len + CDS_len)) %>%
    group_by(transcript) %>%
    top_n(n = 50, wt = Position) %>% #extracts the 3' most nts
    ungroup() %>%
    mutate(nt = (Position - (UTR5_len + CDS_len)) - 1) %>%
    group_by(nt) %>%
    summarise(mean_cpm = mean(CPM)) %>%
    ungroup() %>%
    mutate(region = rep("CDS")) -> CDS_end
  
  #UTR3 start
  merged_data %>%
    filter(Position > (UTR5_len + CDS_len)) %>%
    mutate(nt = Position - (UTR5_len + CDS_len)) %>%
    filter(nt <= 25) %>%
    group_by(nt) %>%
    summarise(mean_cpm = mean(CPM)) %>%
    ungroup() %>% 
    mutate(region = rep("UTR3")) -> UTR3_start
  
  bind_rows(UTR5_end, CDS_start, CDS_end, UTR3_start) %>%
    mutate(sample = rep(sample)) -> spliced_list[[sample]]
}

#plot samples individually----
for (sample in RPF_sample_names) {
  df <- spliced_list[[sample]]
  
  ylims <- c(0,max(df$mean_cpm))
  
  #5'UTR
  df[df$region == "UTR5" & df$nt < 0,] %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    geom_line(size = 1)+
    ylim(ylims)+
    xlab("nt\n(relative to start codon)")+
    ylab("mean CPM")+
    my_theme -> UTR5_end_plot
  
  #CDS
  df[df$region == "CDS" & df$nt > 0,] %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    geom_line(size = 1)+
    ylim(ylims)+
    xlab("nt\n(relative to start codon)")+
    my_theme+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) -> CDS_start_plot
  
  df[df$region == "CDS" & df$nt < 0,] %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    geom_line(size = 1)+
    ylim(ylims)+
    xlab("nt\n(relative to stop codon)")+
    my_theme+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) -> CDS_end_plot
  
  #3'UTR
  df[df$region == "UTR3" & df$nt > 0,] %>%
    ggplot(aes(x = nt, y = mean_cpm))+
    geom_line(size = 1)+
    ylim(ylims)+
    xlab("nt\n(relative to stop codon)")+
    my_theme+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) -> UTR3_start_plot
  
  png(filename = file.path(parent_dir, "plots/offset_aligned_single_nt_plots", paste(sample, "offset aligned single nt plot.png")), width = 1300, height = 300)
  grid.arrange(UTR5_end_plot, CDS_start_plot, CDS_end_plot, UTR3_start_plot, nrow = 1, widths = c(1,2,2,1))
  dev.off()
}

#all mean of all samples----
do.call("rbind", spliced_list) %>%
  group_by(region, nt) %>%
  summarise(mean_cpm = mean(mean_cpm)) -> all_data

ylims <- c(0,max(all_data$mean_cpm))

#5'UTR
all_data[all_data$region == "UTR5" & all_data$nt < 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to start codon)")+
  ylab("mean CPM")+
  my_theme -> UTR5_end_plot

#CDS
all_data[all_data$region == "CDS" & all_data$nt > 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to start codon)")+
  my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) -> CDS_start_plot

all_data[all_data$region == "CDS" & all_data$nt < 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to stop codon)")+
  my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) -> CDS_end_plot

#3'UTR
all_data[all_data$region == "UTR3" & all_data$nt > 0,] %>%
  ggplot(aes(x = nt, y = mean_cpm))+
  geom_line(size = 1)+
  ylim(ylims)+
  xlab("nt\n(relative to stop codon)")+
  my_theme+
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) -> UTR3_start_plot

png(filename = file.path(parent_dir, "plots/offset_aligned_single_nt_plots/all samples offset aligned single nt plot.png"), width = 1300, height = 300)
grid.arrange(UTR5_end_plot, CDS_start_plot, CDS_end_plot, UTR3_start_plot, nrow = 1, widths = c(1,2,2,1))
dev.off()
