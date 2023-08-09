#load libraries----
library(tidyverse)
library(fgsea)
library(ggrepel)

#read in common variables
source("\\\\data.beatson.gla.ac.uk/data/JWALDRON/Ribosome_profiling/SI_epithelial_extractions/AK_A1_A2KOs/paper/GEO_upload/test_run/Scripts/R_scripts/A1/A1_common_variables.R")

#create a variable for what the treatment is----
treatment <- "A1KO"

#themes----
mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 16, hjust = 0, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.position = "none")
#functions----
plot_scatters <- function(DESeq2_df, gsea_result, gsea_set, pathway, dir) {
  
  #extract adjusted p-value and create label from it
  pval <- gsea_result$padj[gsea_result$pathway == pathway]
  plab <- myP(pval)
  
  gene_names <- gsea_set[[pathway]]
  
  DESeq2_df %>%
    mutate(group = factor(gene_sym %in% gene_names),
           alpha_score = case_when(group == T ~ 1,
                                   group == F ~ 0.1)) %>%
    arrange(group) %>%
    ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = group, alpha = alpha_score))+
    geom_point()+
    scale_colour_manual(values=c("grey", "red"))+
    scale_alpha(guide = "none")+
    geom_abline(lty = 2)+
    geom_hline(yintercept = 0, lty = 2)+
    geom_vline(xintercept = 0, lty = 2)+
    ylim(c(-2.5,2.5))+
    xlim(c(-2.5,2.5))+
    mytheme+
    xlab("Total RNA log2FC")+
    ylab("RPFs log2FC")+
    ggtitle(paste(treatment, str_replace_all(pathway, "_", " "), sep = "\n"), subtitle = plab) -> scatter_plot
  
  png(filename = file.path(parent_dir, "plots/fgsea/scatters", dir, paste(treatment, pathway, "TE_scatter_plot.png", sep = "_")), width = 500, height = 500)
  print(scatter_plot)
  dev.off()
}

myP <- function(x) {
  p <- as.numeric(x)
  if (p < 2.2e-16) {
    p_label <- "p adj < 2.2e-16"
  } else {
    if (p < 0.001) {
      rounded_p <- formatC(p, format = "e", digits = 2)
    } else {
      rounded_p <- round(p, digits = 3)
    }
    p_label <- paste("p adj =", rounded_p)
  }
  return(p_label)
}

#read in DESeq2 output----
DESeq2_data <- read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output/merged_DESeq2.csv"))

#read in pathways----
source("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/read_mouse_GSEA_pathways.R") #This may need to be changed to human

#plot all sig pathways----
### The following part will plot overlaid scatters for all significant TE pathways for the hallmarks. It can be edited to do the same for either the RPFs or totals or for different pathways

#load results (the fgsea script needs to have been run first to generate the below files)
load(file = file.path(parent_dir, "Analysis/fgsea/hallmark_results.Rdata"))

#extract pathways
padj <- 0.05 #set adjusted p-value threshold to use
sig_pathways <- hallmark_results$TE$pathway[hallmark_results$TE$padj < padj]

#plot overlaid scatters----
#create directory
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/scatters/hallmark")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/scatters/hallmark"))
}
lapply(sig_pathways, plot_scatters, DESeq2_df = DESeq2_data, gsea_result = hallmark_results$TE, gsea_set = pathways.hallmark, dir = "hallmark")

#plot individual pathways----
### The following part plots the same as above, but with significant gene names annotated, but for individual pathways, allowing some extra control for the pathways that you may want to include in a publication
pathway <- "HALLMARK_E2F_TARGETS"
gene_names <- pathways.hallmark[[pathway]]

DESeq2_data %>%
  filter(!(is.na(TE_group))) %>% #this could be changed to RPFs groups if desired
  mutate(group = factor(case_when(gene_sym %in% gene_names & TE_group == "TE up" ~ "sig up",
                                  gene_sym %in% gene_names & TE_group == "TE down" ~ "sig down",
                                  gene_sym %in% gene_names & (TE_group == "no change" | TE_group == "NS") ~ "NS",
                                  !(gene_sym %in% gene_names) ~ "non pathway"),
                        levels = c("sig down", "NS", "sig up", "non pathway"), ordered = T),
         alpha_score = case_when(group == "sig up" | group == "sig down" ~ 1,
                                 group == "NS" ~ 0.5,
                                 group == "non pathway"  ~ 0.1)) -> plot_df
summary(plot_df)

pval <- hallmark_results$TE$padj[hallmark_results[[3]]$pathway == pathway]
plab <- myP(pval)

plot_df %>%
  arrange(desc(group)) %>%
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = group, alpha = alpha_score))+
  geom_point()+
  scale_alpha(guide = "none")+
  scale_colour_manual(values=c("#d1495b", "#2e4057", "grey"))+
  geom_text_repel(data = plot_df[plot_df$group == "sig down" | plot_df$group == "sig up",], aes(label = gene_sym, colour = group), alpha = 1, size = 4, nudge_x = 1)+
  geom_abline(lty = 2)+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  ylim(c(-3.5,3.5))+
  xlim(c(-3.5,3.5))+
  mytheme+
  theme(plot.subtitle = element_text(size = 18, face = "bold"))+
  xlab("Cytoplasmic RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(treatment, str_replace_all(pathway, "_", " "), sep = "\n"), subtitle = plab) -> annotated_scatter

#create directory
if (!(dir.exists(file.path(parent_dir, "plots/fgsea/scatters/hallmark/annotated")))) {
  dir.create(file.path(parent_dir, "plots/fgsea/scatters/hallmark/annotated"))
}

#export plot
png(filename = file.path(parent_dir, "plots/fgsea/scatters/hallmark/annotated", paste(treatment, pathway, "TE_scatter_plot.png", sep = "_")), width = 400, height = 400)
print(annotated_scatter)
dev.off()

