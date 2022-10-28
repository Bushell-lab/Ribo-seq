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
for (sample in RPF_sample_names) {
  read_csv(file = file.path(parent_dir, "logs", paste0(sample, "_read_counts.csv"))) %>%
    mutate(sample = rep(sample)) %>%
    inner_join(RPF_sample_info, by = "sample") -> data_list[[sample]]
}

RPF_counts <- do.call("rbind", data_list)
write_csv(file = file.path(parent_dir, "logs/RPF_reads_summary.csv"), RPF_counts)

#input read counts----
RPF_counts %>%
  ggplot(aes(x = condition, y = cutadapt_in, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("all counts")+
  ggtitle("RPFs")+
  myTheme -> RPFs_counts_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPF_all_counts.png"), width = 500, height = 300)
print(RPFs_counts_plot)
dev.off()

#trimmed percentages----
RPF_counts %>%
  mutate(cutadapt_perc = (cutadapt_out / cutadapt_in) * 100) %>%
  ggplot(aes(x = condition, y = cutadapt_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("trimmed reads %")+
  ylim(c(0,100))+
  ggtitle("RPFs")+
  myTheme -> RPFs_cutadapt_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPFs_cutadapt_counts.png"), width = 500, height = 300)
print(RPFs_cutadapt_plot)
dev.off()

#alignment percentages----
RPF_counts %>%
  mutate(rRNA_perc = rRNA_out / UMI_clipped_out * 100,
         pc_perc = (pc_out / UMI_clipped_out) * 100,
         tRNA_perc = (tRNA_out / UMI_clipped_out) * 100,
         unaligned_perc = ((pc_in - pc_out) / UMI_clipped_out) * 100,
         sample = str_replace(sample, "_RPFs_", "\n")) %>%
  select(sample, rRNA_perc, pc_perc, tRNA_perc, unaligned_perc) %>%
  gather(key = alignment, value = percentage, rRNA_perc, pc_perc, tRNA_perc, unaligned_perc) %>%
  mutate(alignment = factor(alignment, levels = c("unaligned_perc", "tRNA_perc", "pc_perc", "rRNA_perc"), labels = c("un-aligned", "tRNA", "pc", "rRNA"), ordered = T)) %>%
  ggplot(aes(x = sample, y = percentage, fill = alignment))+
  geom_col()+
  xlab("sample")+
  ylab("% aligments")+
  ggtitle("RPFs")+
  myTheme -> RPF_aligments_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPF_aligments_plot.png"), width = 900, height = 300)
print(RPF_aligments_plot)
dev.off()

#unique pc percentages----
RPF_counts %>%
  mutate(unique_perc = (deduplication_out / deduplication_in) * 100) %>%
  ggplot(aes(x = condition, y = unique_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("Unique pc reads %")+
  ylim(c(0,100))+
  ggtitle("RPFs")+
  myTheme -> RPFs_deduplication_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPFs_deduplication_counts.png"), width = 500, height = 300)
print(RPFs_deduplication_plot)
dev.off()

#final pc counts----
RPF_counts %>%
  ggplot(aes(x = condition, y = deduplication_out, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("Total unique pc reads")+
  ggtitle("RPFs")+
  myTheme -> RPFs_pc_count_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPFs_pc_counts.png"), width = 500, height = 300)
print(RPFs_pc_count_plot)
dev.off()
