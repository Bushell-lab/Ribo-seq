library(tidyverse)
library(grid)
library(gridExtra)

home <- "\\\\data.beatson.gla.ac.uk/data/JWALDRON/External_data/Bornelov"

#read in data----
data_list <- list()
for (sample in c("Self_1", "Self_2", "Self_3", "Self_4", "Diff_1", "Diff_2", "Diff_3", "Diff_4")) {
  df <- read_csv(file = file.path(home, paste0("Analysis/Codon_counts_", sample, "_RPFs_RPFcounts_27-29_len.csv")))
  
  #filter out stop codons and normalise data
  df %>%
    filter(codon != "TGA" & codon != "TAG" & codon != "TAA") %>%
    mutate(minus_4 = (log2(`-4`)) - log2(total / 7),
           minus_3 = (log2(`-3`)) - log2(total / 7),
           E_site = (log2(E)) - log2(total / 7),
           P_site = (log2(P)) - log2(total / 7),
           A_site = (log2(A)) - log2(total / 7),
           plus_1 = (log2(`1`)) - log2(total / 7),
           plus_2 = (log2(`2`)) - log2(total / 7)) %>%
    select(codon, minus_4, minus_3, E_site, P_site, A_site, plus_1, plus_2) %>%
    mutate(sample = factor(rep(sample))) -> data_list[[sample]]
}
data <- do.call("rbind", data_list)

#gather the data into tidy format----
data %>%
  pull(codon) %>%
  unique() -> codons

gathered_list <- list()
for (sample in c("Self_1", "Self_2", "Self_3", "Self_4", "Diff_1", "Diff_2", "Diff_3", "Diff_4")) {
  for (codon in codons) {
    data[data$codon == codon & data$sample == sample,] %>%
      select(-sample) %>%
      gather(key = codon, value = freq) %>%
      rename(position = codon) %>%
      mutate(codon = rep(codon),
             sample = factor(rep(sample))) -> gathered_list[[paste(codon, sample, sep = "_")]]
  }
}
gathered_data <- do.call("rbind", gathered_list)

gathered_data %>%
  mutate(position = as.numeric(case_when(position == "minus_4" ~ -4,
                                         position == "minus_3" ~ -3,
                                         position == "E_site" ~ -2,
                                         position == "P_site" ~ -1,
                                         position == "A_site" ~ 0,
                                         position == "plus_1" ~ 1,
                                         position == "plus_2" ~ 2)),
         wobble = factor(str_sub(codon, 3,3))) -> plot_data

#plot data----
plot_list <- list()
for (sample in RPF_sample_names) {
  plot_data[plot_data$sample == sample,] %>%
    ggplot(aes(x = position, y = freq, colour = codon))+
    geom_point(size = 2) +
    stat_summary(fun=mean, geom="line", size = 1)+
    theme_bw()+
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.position = "none")+
    ylab("normalised codon frequency")+
    scale_x_continuous(limits = c(-4,2), breaks = -4:2)+
    xlab("Codon position (0 = A-site)")+
    ggtitle(sample) -> plot_list[[sample]]
}
png(filename = file.path(home, "plots/all_codons_plot.png"), width = 2000, height = 1000)
grid.arrange(plot_list$Self_1, plot_list$Self_2, plot_list$Self_3, plot_list$Self_4,
             plot_list$Diff_1, plot_list$Diff_2, plot_list$Diff_3, plot_list$Diff_4, nrow = 2)
dev.off()

#calculate delta between samples----
plot_data %>%
  mutate(group = factor(case_when(sample == "Self_1" | sample == "Self_2" | sample == "Self_3" | sample == "Self_4"~ "Self",
                                  sample == "Diff_1" | sample == "Diff_2" | sample == "Diff_3" | sample == "Diff_4"~ "Diff"))) %>%
  group_by(group, codon, position) %>%
  summarise(mean_freq = mean(freq)) %>%
  ungroup() %>%
  spread(key = group, value = mean_freq) %>%
  mutate(delta_freq = Diff - Self,
         wobble = factor(str_sub(codon, 3,3))) -> delta_data

#plot delta----
delta_data %>%
  ggplot(aes(x = position, y = delta_freq, colour = wobble))+
  geom_point() +
  stat_summary(fun=mean, geom="line")+
  scale_x_continuous(limits = c(-4,2), breaks = -4:2)+
  xlab("Codon position (0 = A-site)")+
  ylab("log2FC codon enrichment (Diff/Self)")+
  theme_bw() -> delta_plot

png(filename = file.path(home, "plots/codon_enrichment_delta.png"), width = 500, height = 500)
print(delta_plot)
dev.off()

#plot scatter plots for A and P-sites----
delta_data %>%
  filter(position == 0) %>%
  ggplot(aes(x = Self, y = Diff, colour = wobble))+
  geom_point()+
  geom_abline(lty=2)+
  theme_classic()+
  xlab("A-site codon enrichment (Self)")+
  ylab("A-site codon enrichment (Diff)") -> A_site_scatter_plot

delta_data %>%
  filter(position == -1) %>%
  ggplot(aes(x = Self, y = Diff, colour = wobble))+
  geom_point()+
  geom_abline(lty=2)+
  theme_classic()+
  xlab("P-site codon enrichment (Self)")+
  ylab("P-site codon enrichment (Diff)") -> P_site_scatter_plot

png(filename = file.path(home, "plots/scatter_plots.png"), width = 800, height = 300)
grid.arrange(A_site_scatter_plot, P_site_scatter_plot, nrow = 1)
dev.off()
