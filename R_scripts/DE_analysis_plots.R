#load libraries
library(tidyverse)

#read in common variables
source("common_variables.R")

#create a variable for what the treatment is----
treatment <- "EFT226"

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
RPFs %>%
  select(gene, log2FoldChange, padj) %>%
  rename(RPFs_log2FC = log2FoldChange,
         RPFs_padj = padj) %>%
  inner_join(totals[,c("gene", "log2FoldChange", "padj")], by = "gene") %>%
  rename(totals_log2FC = log2FoldChange,
         totals_padj = padj) %>%
  mutate(group = factor(case_when(RPFs_padj < 0.1 & RPFs_log2FC < 0 & totals_padj >= 0.1 ~ "RPFs down",
                                  RPFs_padj < 0.1 & RPFs_log2FC > 0 & totals_padj >= 0.1 ~ "RPFs up",
                                  RPFs_padj >= 0.1 & totals_padj < 0.1 & totals_log2FC < 0 ~"Totals down",
                                  RPFs_padj >= 0.1 & totals_padj < 0.1 & totals_log2FC > 0 ~"Totals up",
                                  RPFs_padj < 0.1 & totals_padj < 0.1 & RPFs_log2FC < 0 & totals_log2FC < 0 ~ "both down",
                                  RPFs_padj < 0.1 & totals_padj < 0.1 & RPFs_log2FC > 0 & totals_log2FC > 0 ~ "both up")),
         alpha_score = case_when(is.na(group) ~ 0.1,
                                 !(is.na(group)) ~ 1)) -> merged_data
summary(merged_data)


#plot TE scatter----
#caluclate axis limits
RPF_lims <- max(abs(merged_data$RPFs_log2FC))
totals_lims <- max(abs(merged_data$totals_log2FC))
lim <- max(c(RPF_lims, totals_lims))
lims <- c(-lim, lim)

merged_data %>%
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = group, alpha = alpha_score))+
  geom_point()+
  scale_alpha(guide = "none")+
  mytheme+
  xlab("Total RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(treatment, "TE scatter"))+
  xlim(lims)+
  ylim(lims)+
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
write.table(file = file.path(parent_dir, "plots/DE_analysis", paste0(treatment, "_TE_scatter_groups.txt")), summary(merged_data$group), col.names = F, quote = F)

