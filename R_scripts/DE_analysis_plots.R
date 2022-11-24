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
  ggtitle(paste(treatment, "Totals")) -> totals_volcano

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
  ggtitle(paste(treatment, "Totals")) -> totals_MA

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
                           levels = c("TE down", "no change", "TE up", "NS"), ordered = T),
         RPFs_group = factor(case_when(RPFs_padj < RPF_sig_padj & RPFs_log2FC < 0 & totals_padj >= RPF_non_sig_padj ~ "RPFs down",
                                  RPFs_padj < RPF_sig_padj & RPFs_log2FC > 0 & totals_padj >= RPF_non_sig_padj ~ "RPFs up",
                                  RPFs_padj >= RPF_non_sig_padj & totals_padj < RPF_sig_padj & totals_log2FC < 0 ~"Totals down",
                                  RPFs_padj >= RPF_non_sig_padj & totals_padj < RPF_sig_padj & totals_log2FC > 0 ~"Totals up",
                                  RPFs_padj < RPF_sig_padj & totals_padj < RPF_sig_padj & RPFs_log2FC < 0 & totals_log2FC < 0 ~ "both down",
                                  RPFs_padj < RPF_sig_padj & totals_padj < RPF_sig_padj & RPFs_log2FC > 0 & totals_log2FC > 0 ~ "both up"))) -> merged_data
summary(merged_data)


#plot TE scatters----
#based on RPFs/totals logFC
merged_data %>%
  mutate(alpha_score = case_when(is.na(RPFs_group) ~ 0.1,
                                 !(is.na(RPFs_group)) ~ 1)) %>%
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = RPFs_group, alpha = alpha_score))+
  geom_point()+
  scale_alpha(guide = "none")+
  mytheme+
  xlab("Total RNA log2FC")+
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
  geom_vline(xintercept = -1, lty=2) -> RPF_groups_scatter_plot

png(filename = file.path(parent_dir, "plots/DE_analysis", paste(treatment, "_RPF_groups_scatter.png")), width = 500, height = 400)
print(RPF_groups_scatter_plot)
dev.off()

#write out group sizes
write.table(file = file.path(parent_dir, "plots/DE_analysis", paste0(treatment, "_RPF_groups_scatter.txt")), summary(merged_data$RPFs_group), col.names = F, quote = F)

#based on TE
merged_data %>%
  filter(!(is.na(TE_group))) %>%
  mutate(alpha_score = case_when(TE_group =="TE down" | TE_group =="TE up" | TE_group == "no change" ~ 1,
                                 TE_group =="NS" ~ 0.1)) %>%
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = TE_group, alpha = alpha_score))+
  geom_point()+
  scale_colour_manual(values=c("red", "blue", "purple", "grey"))+
  scale_alpha(guide = "none")+
  mytheme+
  xlab("Total RNA log2FC")+
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

#write out group sizes
write.table(file = file.path(parent_dir, "plots/DE_analysis", paste0(treatment, "_TE_scatter_groups.txt")), summary(merged_data$TE_group), col.names = F, quote = F)

#write out csv
write_csv(merged_data, file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2.csv"))

