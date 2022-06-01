#load libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(parallel)

#read in common variables
source("common_variables.R")

#functions
#write a function that will read in a csv file for use with parLapply
read_counts_csv <- function(k){
  df <- read.csv(file = k)
  df$fyle <- rep(k)
  return(df)
}

#write theme
myTheme <- theme_bw()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.title = element_blank())
  
#read in data----
#generate a list of file names
fyle_list <- list()
for(sample in RPF_sample_names) {
  for(i in lengths){
    fyle_list[[paste(sample, i, sep = "_")]] <- file.path(parent_dir, "Analysis/region_counts", paste0(sample, "_pc_L", i, "_Off0_region_counts.csv"))
  }
}

#read in the data with parLapply
no_cores <- detectCores() - 1 #sets the number of cores to use (all but one)
cl <- makeCluster(no_cores) #Initiates cluster
data_list <- parLapply(cl, fyle_list, read_counts_csv) #reads in the data
stopCluster(cl) #Stops cluster

#combine data_list into one data frame
all_data <- do.call("rbind", data_list)

#extract sample and read length from fylenames
all_data %>%
  mutate(read_length = str_remove(fyle, ".+pc_best_L"),
         read_length = as.numeric(str_remove(read_length, "_Off0_region_counts.csv")),
         sample = str_remove(fyle, ".+region_counts/"),
         sample = str_remove(sample, "_pc_best.+")) %>%
  select(-fyle) -> all_data

summary(all_data)

#plot data----
#plot length distribution of all counts
all_data %>%
  group_by(sample, read_length) %>%
  summarise(total_counts = sum(counts, na.rm = T)) -> total_counts

total_counts %>%
  mutate(sample = factor(sample, levels = RPF_sample_names, labels = str_replace_all(RPF_sample_names, "_", " "))) %>%
  ggplot(aes(x = read_length, y = total_counts, colour = sample))+
  geom_line(size = 1)+
  myTheme+
  ylab("total counts")+
  xlab("read length")+
  scale_x_continuous(breaks = seq(min(lengths), max(lengths),2)) -> total_counts_plot

png(filename = file.path(parent_dir, "plots/summed_counts/all_counts_lengths.png"), width = 500, height = 300)
print(total_counts_plot)
dev.off()

#for each sample, plot summed counts within each region across the length distribution
plot_list <- list()
for (sample in RPF_sample_names) {
  all_data[all_data$sample == sample,] %>%
    mutate(region = factor(region, levels = c("UTR3", "UTR5", "CDS"), labels = c("3\'UTR", "5\'UTR", "CDS"), ordered = T)) %>%
    group_by(sample, read_length, region) %>%
    summarise(total_counts = sum(counts, na.rm = T)) %>%
    ggplot(aes(x = read_length, y = total_counts, fill = region))+
    geom_col()+
    myTheme+
    xlab("read length")+
    ylab("total counts")+
    ggtitle(str_replace_all(sample, "_", " ")) -> plot_list[[sample]]
  
  png(filename = file.path(parent_dir, paste0("plots/summed_counts/", sample, "_summed_region_counts.png")), width = 500, height = 300)
  print(plot_list[[sample]])
  dev.off()
}

#for each sample, plot percentage counts within each region across the length distribution
plot_list <- list()
for (sample in RPF_sample_names) {
  all_data[all_data$sample == sample,] %>%
    inner_join(total_counts, by = c("sample", "read_length")) %>%
    mutate(perc_counts = (counts / total_counts * 100),
           region = factor(region, levels = c("UTR3", "UTR5", "CDS"), labels = c("3\'UTR", "5\'UTR", "CDS"), ordered = T)) %>%
    ggplot(aes(x = read_length, y = perc_counts, fill = region))+
    geom_col()+
    myTheme+
    xlab("read length")+
    ylab("% counts")+
    ggtitle(str_replace_all(sample, "_", " ")) -> plot_list[[sample]]
  
  png(filename = file.path(parent_dir, paste0("plots/summed_counts/", sample, "_perc_region_counts.png")), width = 500, height = 300)
  print(plot_list[[sample]])
  dev.off()
}

