#load packages----
library(tidyverse)

#read in common variables
source("common_variables.R")

#themes----
violin_theme <- theme_bw()+
  theme(legend.position='none',
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=18), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=14, face="bold"))

#read in data----
merged_DEseq_data <- read_csv(file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2.csv"))
most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

#read in feature properties----
UTR5_data <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_filtered_UTR5_composition.csv")
CDS_data <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_filtered_CDS_composition.csv")
UTR3_data <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_filtered_UTR3_composition.csv")

#UTR5----
UTR5_data %>%
  inner_join(most_abundant_transcripts, by = "transcript") %>%
  inner_join(merged_DEseq_data, by = "gene") %>%
  mutate(TE_group = factor(TE_group, levels = c("TE down", "no change", "TE up", "NS"), ordered = T),
         RPFs_group = factor(RPFs_group)) -> merged_UTR5_data

summary(merged_UTR5_data)

#length
merged_UTR5_data %>%
  filter(RPFs_group == "RPFs down" | RPFs_group == "RPFs up" | is.na(RPFs_group)) %>%
  ggplot(aes(x = RPFs_group, y = length, fill = RPFs_group))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  scale_y_log10(breaks=c(10,100,1000),limits=c(10, 2000))+
  violin_theme+
  ylab("5\'UTR length (nt)") -> UTR5_length_RPFs_group_plot

png(filename = file.path(parent_dir, "plots/feature_properties/UTR5_length_RPFs_group_plot.png"), width = 500, height = 500)
print(UTR5_length_RPFs_group_plot)
dev.off()

merged_UTR5_data %>%
  filter(!(is.na(TE_group)),
         TE_group != "NS") %>%
  ggplot(aes(x = TE_group, y = length, fill = TE_group))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA,color="white")+
  scale_y_log10(breaks=c(10,100,1000),limits=c(10, 2000))+
  violin_theme+
  ylab("5\'UTR length (nt)") -> UTR5_length_TE_group_plot

png(filename = file.path(parent_dir, "plots/feature_properties/UTR5_length_TE_group_plot.png"), width = 500, height = 500)
print(UTR5_length_TE_group_plot)
dev.off()

#GC content
merged_UTR5_data %>%
  filter(RPFs_group == "RPFs down" | RPFs_group == "RPFs up" | is.na(RPFs_group)) %>%
  mutate(GC_content = G_content + C_content) %>%
  ggplot(aes(x = RPFs_group, y = GC_content, fill = RPFs_group))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  violin_theme+
  ylab("5\'UTR GC content (%)") -> UTR5_GC_RPFs_group_plot

png(filename = file.path(parent_dir, "plots/feature_properties/UTR5_GC_RPFs_group_plot.png"), width = 500, height = 500)
print(UTR5_GC_RPFs_group_plot)
dev.off()

merged_UTR5_data %>%
  filter(!(is.na(TE_group)),
         TE_group != "NS") %>%
  mutate(GC_content = G_content + C_content) %>%
  ggplot(aes(x = TE_group, y = GC_content, fill = TE_group))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA, color="white")+
  violin_theme+
  ylab("5\'UTR GC content (%)") -> UTR5_GC_TE_group_plot

png(filename = file.path(parent_dir, "plots/feature_properties/UTR5_GC_TE_group_plot.png"), width = 500, height = 500)
print(UTR5_GC_TE_group_plot)
dev.off()

#CDS----
CDS_data %>%
  inner_join(most_abundant_transcripts, by = "transcript") %>%
  inner_join(merged_DEseq_data, by = "gene") %>%
  mutate(TE_group = factor(TE_group, levels = c("TE down", "no change", "TE up", "NS"), ordered = T),
         RPFs_group = factor(RPFs_group)) -> merged_CDS_data

summary(merged_CDS_data)

#length
merged_CDS_data %>%
  filter(RPFs_group == "RPFs down" | RPFs_group == "RPFs up" | is.na(RPFs_group)) %>%
  ggplot(aes(x = RPFs_group, y = length, fill = RPFs_group))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  scale_y_log10(breaks=c(100,1000,10000),limits=c(100, 10000))+
  violin_theme+
  ylab("CDS length (nt)") -> CDS_length_RPFs_group_plot

png(filename = file.path(parent_dir, "plots/feature_properties/CDS_length_RPFs_group_plot.png"), width = 500, height = 500)
print(CDS_length_RPFs_group_plot)
dev.off()

merged_CDS_data %>%
  filter(!(is.na(TE_group)),
         TE_group != "NS") %>%
  ggplot(aes(x = TE_group, y = length, fill = TE_group))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA, color="white")+
  scale_y_log10(breaks=c(100,1000,10000),limits=c(100, 10000))+
  violin_theme+
  ylab("CDS length (nt)") -> CDS_length_TE_group_plot

png(filename = file.path(parent_dir, "plots/feature_properties/CDS_length_TE_group_plot.png"), width = 500, height = 500)
print(CDS_length_TE_group_plot)
dev.off()

#GC content
merged_CDS_data %>%
  filter(RPFs_group == "RPFs down" | RPFs_group == "RPFs up" | is.na(RPFs_group)) %>%
  mutate(GC_content = G_content + C_content) %>%
  ggplot(aes(x = RPFs_group, y = GC_content, fill = RPFs_group))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  violin_theme+
  ylab("CDS GC content (%)") -> CDS_GC_RPFs_group_plot

png(filename = file.path(parent_dir, "plots/feature_properties/CDS_GC_RPFs_group_plot.png"), width = 500, height = 500)
print(CDS_GC_RPFs_group_plot)
dev.off()

merged_CDS_data %>%
  filter(!(is.na(TE_group)),
         TE_group != "NS") %>%
  mutate(GC_content = G_content + C_content) %>%
  ggplot(aes(x = TE_group, y = GC_content, fill = TE_group))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA, color="white")+
  violin_theme+
  ylab("CDS GC content (%)") -> CDS_GC_TE_group_plot

png(filename = file.path(parent_dir, "plots/feature_properties/CDS_GC_TE_group_plot.png"), width = 500, height = 500)
print(CDS_GC_TE_group_plot)
dev.off()

#UTR3----
UTR3_data %>%
  inner_join(most_abundant_transcripts, by = "transcript") %>%
  inner_join(merged_DEseq_data, by = "gene") %>%
  mutate(TE_group = factor(TE_group, levels = c("TE down", "no change", "TE up", "NS"), ordered = T),
         RPFs_group = factor(RPFs_group)) -> merged_UTR3_data

summary(merged_UTR3_data)

#length
merged_UTR3_data %>%
  filter(RPFs_group == "RPFs down" | RPFs_group == "RPFs up" | is.na(RPFs_group)) %>%
  ggplot(aes(x = RPFs_group, y = length, fill = RPFs_group))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  scale_y_log10(breaks=c(100,1000,10000),limits=c(10, 10000))+
  violin_theme+
  ylab("3\'UTR length (nt)") -> UTR3_length_RPFs_group_plot

png(filename = file.path(parent_dir, "plots/feature_properties/UTR3_length_RPFs_group_plot.png"), width = 500, height = 500)
print(UTR3_length_RPFs_group_plot)
dev.off()

merged_UTR3_data %>%
  filter(!(is.na(TE_group)),
         TE_group != "NS") %>%
  ggplot(aes(x = TE_group, y = length, fill = TE_group))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA, color="white")+
  scale_y_log10(breaks=c(100,1000,10000),limits=c(10, 10000))+
  violin_theme+
  ylab("3\'UTR length (nt)") -> UTR3_length_TE_group_plot

png(filename = file.path(parent_dir, "plots/feature_properties/UTR3_length_TE_group_plot.png"), width = 500, height = 500)
print(UTR3_length_TE_group_plot)
dev.off()

#GC content
merged_UTR3_data %>%
  filter(RPFs_group == "RPFs down" | RPFs_group == "RPFs up" | is.na(RPFs_group)) %>%
  mutate(GC_content = G_content + C_content) %>%
  ggplot(aes(x = RPFs_group, y = GC_content, fill = RPFs_group))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA)+
  violin_theme+
  ylab("3\'UTR GC content (%)") -> UTR3_GC_RPFs_group_plot

png(filename = file.path(parent_dir, "plots/feature_properties/UTR3_GC_RPFs_group_plot.png"), width = 500, height = 500)
print(UTR3_GC_RPFs_group_plot)
dev.off()

merged_UTR3_data %>%
  filter(!(is.na(TE_group)),
         TE_group != "NS") %>%
  mutate(GC_content = G_content + C_content) %>%
  ggplot(aes(x = TE_group, y = GC_content, fill = TE_group))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(width = 0.2, outlier.shape=NA, color="white")+
  violin_theme+
  ylab("3\'UTR GC content (%)") -> UTR3_GC_TE_group_plot

png(filename = file.path(parent_dir, "plots/feature_properties/UTR3_GC_TE_group_plot.png"), width = 500, height = 500)
print(UTR3_GC_TE_group_plot)
dev.off()

