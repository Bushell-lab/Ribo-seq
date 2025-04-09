#load libraries
library(tidyverse)

#read in common variables
source("common_variables.R")

myTheme <- theme_classic()+
  theme(axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 18),
        legend.title = element_blank())


#read in read counts summaries----
data_list <- list()
for (sample in Total_sample_names) {
  read_csv(file = file.path(parent_dir, "logs", paste0(sample, "_read_counts.csv"))) %>%
    mutate(sample = rep(sample)) %>%
    inner_join(Total_sample_info, by = "sample") -> data_list[[sample]]
}

Totals_counts <- do.call("rbind", data_list)

write_csv(file = file.path(parent_dir, "logs/Totals_reads_summary.csv"), Totals_counts)

#plot total counts----
Totals_counts %>%
  ggplot(aes(x = condition, y = cutadapt_in, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("all counts")+
  ggtitle("Totals")+
  myTheme -> Totals_counts_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/Totals_all_counts.png"), width = 500, height = 300)
print(Totals_counts_plot)
dev.off()

#plot cutadapt counts %----
Totals_counts %>%
  mutate(cutadapt_perc = (cutadapt_out / cutadapt_in) * 100) %>%
  ggplot(aes(x = condition, y = cutadapt_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("trimmed reads %")+
  ylim(c(0,100))+
  ggtitle("Totals")+
  myTheme -> Totals_cutadapt_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/Totals_cutadapt_counts.png"), width = 500, height = 300)
print(Totals_cutadapt_plot)
dev.off()

#pc %----
Totals_counts %>%
  mutate(pc_perc = (pc_out / pc_in) * 100) %>%
  ggplot(aes(x = condition, y = pc_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("pc reads %")+
  ylim(c(0,100))+
  ggtitle("Totals")+
  myTheme -> Totals_pc_alignment_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/Totals_pc_alignments.png"), width = 500, height = 300)
print(Totals_pc_alignment_plot)
dev.off()

#plot unique counts %----
Totals_counts %>%
  mutate(unique_perc = (deduplication_out / deduplication_in) * 100) %>%
  ggplot(aes(x = condition, y = unique_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("Unique pc reads %")+
  ylim(c(0,100))+
  ggtitle("Totals")+
  myTheme -> Totals_deduplication_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/Totals_deduplication_counts.png"), width = 500, height = 300)
print(Totals_deduplication_plot)
dev.off()

#final pc counts----
Totals_counts %>%
  ggplot(aes(x = condition, y = deduplication_out, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("Total unique pc reads")+
  ggtitle("Totals")+
  myTheme -> Totals_pc_count_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/Totals_pc_counts.png"), width = 500, height = 300)
print(Totals_pc_count_plot)
dev.off()

