#load libraries
library(tidyverse)

#read in common variables
source("common_variables.R")

#create a variable for what the treatment is----
treatment <- "KO"

#themes----
mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 16))

#read in DESeq2 output----
totals <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("Totals_", treatment, "_DEseq2_apeglm_LFC_shrinkage.csv")))
RPFs <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("RPFs_", treatment, "_DEseq2_apeglm_LFC_shrinkage.csv")))
TE <- read_csv((file = file.path(parent_dir, "Analysis/DESeq2_output", paste0("TE_", treatment, "_DEseq2.csv"))))

#plot volcanos----
RPFs %>%
  filter(!(is.na(padj))) %>%
  mutate(sig = factor(case_when(padj < 0.1 ~ "*",
                                padj >= 0.1 ~ "NS"))) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
  geom_point(alpha = 0.5)+
  mytheme+
  xlab("log2FC")+
  ylab("-log10(padj)")+
  ggtitle(paste(treatment, "RPFs")) -> RPFs_volcano

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("RPFs_", treatment, "_volcano.png")), width = 300, height = 400)
print(RPFs_volcano)
dev.off()

totals %>%
  filter(!(is.na(padj))) %>%
  mutate(sig = factor(case_when(padj < 0.1 ~ "*",
                                padj >= 0.1 ~ "NS"))) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
  geom_point(alpha = 0.5)+
  mytheme+
  xlab("log2FC")+
  ylab("-log10(padj)")+
  ggtitle(paste(treatment, "Cytoplasmic RNA")) -> totals_volcano

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("Totals_", treatment, "_volcano.png")), width = 300, height = 400)
print(totals_volcano)
dev.off()

#plot MAs----
RPFs %>%
  filter(!(is.na(padj))) %>%
  mutate(sig = factor(case_when(padj < 0.1 ~ "*",
                                padj >= 0.1 ~ "NS"))) %>%
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, colour = sig))+
  geom_point(alpha = 0.5)+
  mytheme+
  xlab("log10(mean expression)")+
  ylab("log2FC")+
  ggtitle(paste(treatment, "RPFs")) -> RPFs_MA

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("RPFs_", treatment, "_MA.png")), width = 400, height = 300)
print(RPFs_MA)
dev.off()

totals %>%
  filter(!(is.na(padj))) %>%
  mutate(sig = factor(case_when(padj < 0.1 ~ "*",
                                padj >= 0.1 ~ "NS"))) %>%
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, colour = sig))+
  geom_point(alpha = 0.5)+
  mytheme+
  xlab("log10(mean expression)")+
  ylab("log2FC")+
  ggtitle(paste(treatment, "Cytoplasmic RNA")) -> totals_MA

png(filename = file.path(parent_dir, "plots/DE_analysis", paste0("Totals_", treatment, "_MA.png")), width = 400, height = 300)
print(totals_MA)
dev.off()

#merge RPF with totals data----
#select apdj thresholds
TE_sig_padj <- 0.1
RPF_sig_padj <- 0.1

TE_non_sig_padj <- 0.9
RPF_non_sig_padj <- 0.5

#merged data and make groups based on RPF/Total adjusted p-values or TE adjusted p-values
RPFs %>%
  select(gene, gene_sym, log2FoldChange, padj) %>%
  rename(RPFs_log2FC = log2FoldChange,
         RPFs_padj = padj) %>%
  inner_join(totals[,c("gene", "log2FoldChange", "padj", "gene_sym")], by = c("gene", "gene_sym")) %>%
  rename(totals_log2FC = log2FoldChange,
         totals_padj = padj) %>%
  inner_join(TE[,c("gene","log2FoldChange", "padj")], by = "gene") %>%
  rename(TE_log2FC = log2FoldChange,
         TE_padj = padj) %>%
  mutate(TE_group = factor(case_when(TE_padj < TE_sig_padj & TE_log2FC < 0 ~ "TE down",
                                     TE_padj < TE_sig_padj & TE_log2FC > 0 ~ "TE up",
                                     TE_padj >= TE_non_sig_padj ~ "no change",
                                     (TE_padj >= TE_sig_padj & TE_padj <= TE_non_sig_padj) | is.na(TE_padj) ~ "NS"),
                           levels = c("TE down", "TE up", "no change", "NS"), ordered = T),
         RPFs_group = factor(case_when((RPFs_padj < RPF_sig_padj & RPFs_log2FC < -log2FC_threshold) & (totals_padj >= RPF_non_sig_padj | totals_log2FC > log2FC_threshold) ~ "RPFs down",
                                       (RPFs_padj < RPF_sig_padj & RPFs_log2FC > log2FC_threshold) & (totals_padj >= RPF_non_sig_padj | totals_log2FC < -log2FC_threshold)  ~ "RPFs up",
                                       (totals_padj < RPF_sig_padj & totals_log2FC < -log2FC_threshold) & (RPFs_padj >= RPF_non_sig_padj | RPFs_log2FC > log2FC_threshold)  ~ "Totals down",
                                       (totals_padj < RPF_sig_padj & totals_log2FC > log2FC_threshold) & (RPFs_padj >= RPF_non_sig_padj | RPFs_log2FC < -log2FC_threshold)  ~ "Totals up",
                                       RPFs_padj < RPF_sig_padj & totals_padj < RPF_sig_padj & RPFs_log2FC < -log2FC_threshold & totals_log2FC < -log2FC_threshold ~ "both down",
                                       RPFs_padj < RPF_sig_padj & totals_padj < RPF_sig_padj & RPFs_log2FC > log2FC_threshold & totals_log2FC > log2FC_threshold ~ "both up",
                                       totals_padj >= RPF_non_sig_padj & RPFs_padj >= RPF_non_sig_padj ~ "no change")),
         RPFs_group = factor(case_when(is.na(RPFs_group) ~ "NS",
                                       !(is.na(RPFs_group)) ~ RPFs_group), levels = c("RPFs down", "RPFs up", "Totals down", "Totals up", "both down", "both up", "no change", "NS"), ordered = T)) -> merged_data
summary(merged_data)


#plot TE scatters----
#based on RPFs/totals logFC
#add "n=" labels
merged_data %>%
  group_by(RPFs_group) %>%
  summarize(num = n()) %>%
  mutate(lab = paste0(RPFs_group, " (n=", num, ") ")) -> labs

merged_data %>%
  mutate(alpha_score = case_when(RPFs_group == "NS" | RPFs_group == "no change" ~ 0.1,
                                 RPFs_group != "NS" & RPFs_group != "no change" ~ 1)) %>%
  arrange(desc(RPFs_group)) %>%
  left_join(labs, by = "RPFs_group") %>%
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = lab, alpha = alpha_score))+
  geom_point()+
  scale_alpha(guide = "none")+
  scale_color_manual(values=c("#A6CEE3", "#1F78B4", "#2e4057", "grey", "#FB9A99", "#33A02C", "#B2DF8A", "#E31A1C"))+
  mytheme+
  xlab("Cytoplasmic RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(treatment, "TE scatter"))+
  xlim(c(-5,5))+
  ylim(c(-5,5))+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_vline(xintercept = 0, lty=1) -> RPF_groups_scatter_plot

png(filename = file.path(parent_dir, "plots/DE_analysis", paste(treatment, "_RPF_groups_scatter.png")), width = 500, height = 400)
print(RPF_groups_scatter_plot)
dev.off()

#based on TE
#add "n=" labels
merged_data %>%
  group_by(TE_group) %>%
  summarize(num = n()) %>%
  mutate(lab = paste0(TE_group, "\n(n=", num, ")\n")) -> labs

merged_data %>%
  filter(!(is.na(TE_group))) %>%
  arrange(desc(TE_group)) %>%
  mutate(alpha_score = case_when(TE_group =="TE down" | TE_group =="TE up" ~ 1,
                                 TE_group == "no change" | TE_group =="NS" ~ 0.1)) %>%
  left_join(labs, by = "TE_group") %>%
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = lab, alpha = alpha_score))+
  geom_point()+
  scale_colour_manual(values=c("#2e4057", "grey", "#d1495b", "purple"))+
  scale_alpha(guide = "none")+
  mytheme+
  xlab("Cytoplasmic RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(treatment, "TE scatter"))+
  xlim(c(-5,5))+
  ylim(c(-5,5))+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_hline(yintercept = 1, lty=2)+
  geom_hline(yintercept = -1, lty=2)+
  geom_vline(xintercept = 0, lty=1)+
  geom_vline(xintercept = 1, lty=2)+
  geom_vline(xintercept = -1, lty=2) -> TE_scatter_plot

png(filename = file.path(parent_dir, "plots/DE_analysis", paste(treatment, "_TE_scatter.png")), width = 500, height = 400)
print(TE_scatter_plot)
dev.off()

#write out csv
write_csv(merged_data, file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2.csv"))

