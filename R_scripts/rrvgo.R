#load libraries----
library(rrvgo)
library(tidyverse)

#read in common variables----
source("common_variables.R")

#set threshold----
#This determines how many parent terms will be defined when reducing the GO terms
#scale of 0.1-0.9 with 0.9 have the the most reduction or lowest amount of parent terms. 0.7 is the standard analysis
threshold <- 0.7

#set adjusted p-value threshold----
#This determines which GO terms are deemed statistically significant
padj_threshold <- 0.05

#read in data----
##GO term IDs----
GO_terms <- read.table("\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/GSEA/go_term_to_id.tsv", header = T)

#subset GO_terms for the ontology
GO_terms %>%
  filter(str_detect(pathway, "GOBP")) -> BP_GO_terms
GO_terms %>%
  filter(str_detect(pathway, "GOCC")) -> CC_GO_terms
GO_terms %>%
  filter(str_detect(pathway, "GOMF")) -> MF_GO_terms

##FGSEA results----
load(file.path(parent_dir, "Analysis/fgsea/bio_processes_results.Rdata"))
load(file.path(parent_dir, "Analysis/fgsea/cell_comp_results.Rdata"))
load(file.path(parent_dir, "Analysis/fgsea/mol_funs_results.Rdata"))

#extract the TE output as a dataframe (this is the 3rd item in the fgsea results list, can change to RPFs [[1]] or Totals [[2]] if required)
BP_df <- as.data.frame(bio_processes_results[[3]])
CC_df <- as.data.frame(cell_comp_results[[3]])
MF_df <- as.data.frame(mol_funs_results[[3]])

#merge and filter data----
#merge fgsea data with GO terms, select GO_ID, adjusted p-value and NES columns and filter for significant adjusted p-values
BP_GO_terms %>%
  inner_join(BP_df, by = "pathway") %>%
  select(GO_ID, padj, NES) %>%
  filter(padj < padj_threshold) -> BP_data

CC_GO_terms %>%
  inner_join(CC_df, by = "pathway") %>%
  select(GO_ID, padj, NES) %>%
  filter(padj < padj_threshold) -> CC_data

MF_GO_terms %>%
  inner_join(MF_df, by = "pathway") %>%
  select(GO_ID, padj, NES) %>%
  filter(padj < padj_threshold) -> MF_data

#subset into up and down
BP_data %>%
  filter(NES > 0) -> BP_up

BP_data %>%
  filter(NES < 0) -> BP_down

CC_data %>%
  filter(NES > 0) -> CC_up

CC_data %>%
  filter(NES < 0) -> CC_down

MF_data %>%
  filter(NES > 0) -> MF_up

MF_data %>%
  filter(NES < 0) -> MF_down

#make sim matrices----
BP_simMatrix_up <- calculateSimMatrix(BP_up$GO_ID,
                                      orgdb="org.Mm.eg.db",
                                      ont="BP",
                                      method="Rel")

BP_simMatrix_down <- calculateSimMatrix(BP_down$GO_ID,
                                        orgdb="org.Mm.eg.db",
                                        ont="BP",
                                        method="Rel")

CC_simMatrix_up <- calculateSimMatrix(CC_up$GO_ID,
                                      orgdb="org.Mm.eg.db",
                                      ont="CC",
                                      method="Rel")

CC_simMatrix_down <- calculateSimMatrix(CC_down$GO_ID,
                                        orgdb="org.Mm.eg.db",
                                        ont="CC",
                                        method="Rel")

MF_simMatrix_up <- calculateSimMatrix(MF_up$GO_ID,
                                      orgdb="org.Mm.eg.db",
                                      ont="MF",
                                      method="Rel")

MF_simMatrix_down <- calculateSimMatrix(MF_down$GO_ID,
                                        orgdb="org.Mm.eg.db",
                                        ont="MF",
                                        method="Rel")

#create scores----
BP_scores_up <- setNames(-log10(BP_up$padj), BP_up$GO_ID)
BP_scores_down <- setNames(-log10(BP_down$padj), BP_down$GO_ID)

CC_scores_up <- setNames(-log10(CC_up$padj), CC_up$GO_ID)
CC_scores_down <- setNames(-log10(CC_down$padj), CC_down$GO_ID)

MF_scores_up <- setNames(-log10(MF_up$padj), MF_up$GO_ID)
MF_scores_down <- setNames(-log10(MF_down$padj), MF_down$GO_ID)

#determine reduced GO terms
BP_reducedTerms_up <- reduceSimMatrix(BP_simMatrix_up,
                                      BP_scores_up,
                                      threshold=threshold,
                                      orgdb="org.Mm.eg.db")

BP_reducedTerms_down <- reduceSimMatrix(BP_simMatrix_down,
                                        BP_scores_down,
                                        threshold=threshold,
                                        orgdb="org.Mm.eg.db")

CC_reducedTerms_up <- reduceSimMatrix(CC_simMatrix_up,
                                      CC_scores_up,
                                      threshold=threshold,
                                      orgdb="org.Mm.eg.db")

CC_reducedTerms_down <- reduceSimMatrix(CC_simMatrix_down,
                                        CC_scores_down,
                                        threshold=threshold,
                                        orgdb="org.Mm.eg.db")

MF_reducedTerms_up <- reduceSimMatrix(MF_simMatrix_up,
                                      MF_scores_up,
                                      threshold=threshold,
                                      orgdb="org.Mm.eg.db")

MF_reducedTerms_down <- reduceSimMatrix(MF_simMatrix_down,
                                        MF_scores_down,
                                        threshold=threshold,
                                        orgdb="org.Mm.eg.db")

#write out reduced terms----
write_csv(BP_reducedTerms_up, file = file.path(parent_dir, "Analysis/fgsea/bio_processes_reducedTerms_up.csv"))
write_csv(BP_reducedTerms_down, file = file.path(parent_dir, "Analysis/fgsea/bio_processes_reducedTerms_down.csv"))

write_csv(CC_reducedTerms_down, file = file.path(parent_dir, "Analysis/fgsea/cell_comp_reducedTerms_down.csv"))
write_csv(CC_reducedTerms_down, file = file.path(parent_dir, "Analysis/fgsea/cell_comp_reducedTerms_down.csv"))

write_csv(MF_reducedTerms_down, file = file.path(parent_dir, "Analysis/fgsea/mol_funs_reducedTerms_down.csv"))
write_csv(MF_reducedTerms_down, file = file.path(parent_dir, "Analysis/fgsea/mol_funs_reducedTerms_down.csv"))

#plot treemaps----
png(filename = file.path(parent_dir, "plots/fgsea/rrvgo/bio_processes_up_treemap.png"), width = 500, height =500 )
treemap <- treemapPlot(BP_reducedTerms_up, size = "score") #size gives size of GO term, score is according to the score you define (pval)
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/rrvgo/bio_processes_down_treemap.png"), width = 500, height =500 )
treemap <- treemapPlot(BP_reducedTerms_down, size = "score") #size gives size of GO term, score is according to the score you define (pval)
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/rrvgo/cell_comp_up_treemap.png"), width = 500, height =500 )
treemap <- treemapPlot(CC_reducedTerms_up, size = "score") #size gives size of GO term, score is according to the score you define (pval)
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/rrvgo/cell_comp_down_treemap.png"), width = 500, height =500 )
treemap <- treemapPlot(CC_reducedTerms_down, size = "score") #size gives size of GO term, score is according to the score you define (pval)
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/rrvgo/mol_funs_up_treemap.png"), width = 500, height =500 )
treemap <- treemapPlot(MF_reducedTerms_up, size = "score") #size gives size of GO term, score is according to the score you define (pval)
dev.off()

png(filename = file.path(parent_dir, "plots/fgsea/rrvgo/mol_funs_down_treemap.png"), width = 500, height =500 )
treemap <- treemapPlot(MF_reducedTerms_down, size = "score") #size gives size of GO term, score is according to the score you define (pval)
dev.off()

#make bar charts----
BP_reducedTerms_down %>%
  group_by(parentTerm) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  mutate(score = -score) -> BP_reducedTerms_down_summarised

BP_reducedTerms_up %>%
  group_by(parentTerm) %>%
  summarise(score = mean(score)) %>%
  ungroup() -> BP_reducedTerms_up_summarised

BP_reducedTerms_down_summarised %>%
  bind_rows(BP_reducedTerms_up_summarised) -> BP_reducedTerms_summarised

BP_reducedTerms_summarised %>%
  arrange(score) %>%
  pull(parentTerm) -> order_BP_terms

BP_reducedTerms_summarised%>%
  ggplot(aes(x = factor(parentTerm, levels = order_BP_terms, ordered = T), y = score))+
  geom_col()+
  coord_flip()+
  theme_classic()+
  theme(axis.title = element_blank()) -> summarised_terms_plot

png(filename = file.path(parent_dir, "plots/fgsea/rrvgo/bio_processes_summarised_terms_plot.png"), width = 450, height =350 )
print(summarised_terms_plot)
dev.off()

