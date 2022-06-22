#load libraries
library(tidyverse)

#read in common variables
source("common_variables.R")

myTheme <- theme_classic()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


#read in fastq line counts csv----
line_counts <- read_csv(file = file.path(parent_dir, "logs/fastq_lines.csv"), col_names = c("lines", "fylename"))

#remove the non_rRNA/tRNA/mito filenames and divide file lines by 4 to get read counts----
line_counts %>%
  filter(!(str_detect(fylename, "non"))) %>% 
  mutate(counts = lines / 4) %>%
  select(-lines) -> read_counts

#extract sample info from file names----
fylename_list <- list()
for (i in 1:nrow(read_counts)) {
  df <- read_counts[i,]
  
  fylename <- df$fylename
  split_fylename <- str_remove(fylename, ".+\\/")
  split_fylename <- str_remove(split_fylename, "\\.fastq")
  split_fylename <- str_split(split_fylename, "_")
  
  fylename_list[[fylename]] <- data.frame(fylename = fylename,
                            condition = split_fylename[[1]][1],
                            library = split_fylename[[1]][2],
                            replicate = split_fylename[[1]][3],
                            stage = split_fylename[[1]][4])
  
}
fylename_df <- do.call("rbind", fylename_list)

#merge data and replace NAs with "all" (these are the starting fastq files so don't have a final value for stage)
read_counts %>%
  inner_join(fylename_df, by = "fylename") %>%
  mutate(stage = replace(stage, is.na(stage), "all")) %>%
  select(-fylename) -> merged_df

merged_df %>%
filter(library == "RPFs") %>%
  spread(value = counts, key = stage) -> RPF_counts

write_csv(file = file.path(parent_dir, "Analysis/reads_summary/RPF_reads_summary.csv"), RPF_counts)

merged_df %>%
  filter(library == "Totals") %>%
  spread(value = counts, key = stage) -> Totals_counts

write_csv(file = file.path(parent_dir, "Analysis/reads_summary/Totals_reads_summary.csv"), Totals_counts)

#plot total counts----
RPF_counts %>%
  ggplot(aes(x = condition, y = all, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("all counts")+
  ggtitle("RPFs counts")+
  myTheme -> RPFs_counts_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPFs_all_counts.png"), width = 500, height = 300)
print(RPFs_counts_plot)
dev.off()

Totals_counts %>%
  ggplot(aes(x = condition, y = all, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("all counts")+
  ggtitle("Totals counts")+
  myTheme -> Totals_counts_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/Totals_all_counts.png"), width = 500, height = 300)
print(Totals_counts_plot)
dev.off()

#plot cutadapt counts %----
RPF_counts %>%
  mutate(unique_perc = (cutadapt / all) * 100) %>%
  ggplot(aes(x = condition, y = unique_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("trimmed reads %")+
  ylim(c(0,100))+
  ggtitle("RPF cutadapt")+
  myTheme -> RPFs_cutadapt_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPFs_cutadapt_counts.png"), width = 500, height = 300)
print(RPFs_cutadapt_plot)
dev.off()

Totals_counts %>%
  mutate(unique_perc = (cutadapt / all) * 100) %>%
  ggplot(aes(x = condition, y = unique_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("trimmed reads %")+
  ylim(c(0,100))+
  ggtitle("Totals cutadapt")+
  myTheme -> Totals_cutadapt_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/Totals_cutadapt_counts.png"), width = 500, height = 300)
print(Totals_cutadapt_plot)
dev.off()

#plot unique counts %----
RPF_counts %>%
  mutate(unique_perc = (cdhitdup / cutadapt) * 100) %>%
  ggplot(aes(x = condition, y = unique_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("Unique reads %")+
  ylim(c(0,100))+
  ggtitle("RPF de-duplication")+
  myTheme -> RPFs_deduplication_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPFs_deduplication_counts.png"), width = 500, height = 300)
print(RPFs_deduplication_plot)
dev.off()

Totals_counts %>%
  mutate(unique_perc = (cdhitdup / cutadapt) * 100) %>%
  ggplot(aes(x = condition, y = unique_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("Unique reads %")+
  ylim(c(0,100))+
  ggtitle("Totals de-duplication")+
  myTheme -> Totals_deduplication_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/Totals_deduplication_counts.png"), width = 500, height = 300)
print(Totals_deduplication_plot)
dev.off()

#rRNA %----
RPF_counts %>%
  mutate(rRNA_perc = (rRNA / cdhitdup) * 100) %>%
  ggplot(aes(x = condition, y = rRNA_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("rRNA reads %")+
  ylim(c(0,100))+
  ggtitle("RPFs rRNA")+
  myTheme -> RPFs_rRNA_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPFs_rRNA_counts.png"), width = 500, height = 300)
print(RPFs_rRNA_plot)
dev.off()

Totals_counts %>%
  mutate(rRNA_perc = (rRNA / cdhitdup) * 100) %>%
  ggplot(aes(x = condition, y = rRNA_perc, fill = replicate))+
  geom_col(position = position_dodge())+
  ylab("rRNA reads %")+
  ylim(c(0,100))+
  ggtitle("Totals rRNA")+
  myTheme -> Totals_rRNA_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/Totals_rRNA_counts.png"), width = 500, height = 300)
print(Totals_rRNA_plot)
dev.off()

#aligned reads----
RPF_counts %>%
  mutate(non_rRNA_counts = cdhitdup - rRNA,
         pc_perc = (pc / non_rRNA_counts) * 100,
         tRNA_perc = (tRNA / non_rRNA_counts) * 100,
         mito_perc = (mito / non_rRNA_counts) * 100,
         unaligned_perc = (unaligned / non_rRNA_counts) * 100,
         sample = paste(condition, replicate, sep = " ")) %>%
  select(sample, pc_perc, tRNA_perc, mito_perc, unaligned_perc) %>%
  gather(key = alignment, value = percentage, pc_perc, tRNA_perc, mito_perc, unaligned_perc) %>%
  mutate(alignment = factor(alignment, levels = c("unaligned_perc", "mito_perc", "tRNA_perc", "pc_perc"), labels = c("un-aligned", "mito", "tRNA", "pc"), ordered = T)) %>%
  ggplot(aes(x = sample, y = percentage, fill = alignment))+
  geom_col()+
  xlab("sample")+
  ylab("% aligments")+
  ggtitle("RPFs non-rRNA % alignments")+
  myTheme -> RPF_aligments_plot

png(filename = file.path(parent_dir, "plots/read_counts_summary/RPF_aligments_plot.png"), width = 900, height = 300)
print(RPF_aligments_plot)
dev.off()

