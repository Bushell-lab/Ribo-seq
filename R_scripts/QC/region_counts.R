#load libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(parallel)

#read in common variables
source("common_variables.R")

#set the read lengths you wish to plot
lengths <- 25:35

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
#extract sample and read length from fyle and inner_join with sample info to get condition and replicate
do.call("rbind", data_list) %>%
  mutate(read_length = str_remove(fyle, ".+pc_L"),
         read_length = as.numeric(str_remove(read_length, "_Off0_region_counts.csv")),
         sample = str_remove(fyle, ".+region_counts/"),
         sample = factor(str_remove(sample, "_pc.+"))) %>%
  inner_join(RPF_sample_info, by = "sample") %>%
  select(-fyle) -> all_data

summary(all_data)

#plot data----
#plot length distribution of all counts
#summed counts
all_data %>%
  group_by(condition, replicate, read_length) %>%
  summarise(read_length_counts = sum(counts, na.rm = T)) %>%
  ggplot(aes(x = read_length, y = read_length_counts, colour = condition, lty = replicate, shape = replicate))+
  geom_line(size = 1)+
  geom_point(size = 2)+
  myTheme+
  ylab("total counts")+
  xlab("read length")+
  scale_x_continuous(breaks = seq(min(lengths), max(lengths),2)) -> total_counts_plot

png(filename = file.path(parent_dir, "plots/summed_counts/lengths_count_plot.png"), width = 500, height = 300)
print(total_counts_plot)
dev.off()

#percentage counts
all_data %>%
  group_by(condition, replicate) %>%
  summarise(total_counts = sum(counts, na.rm = T)) -> total_counts

all_data %>%
  group_by(condition, replicate, read_length) %>%
  summarise(read_length_counts = sum(counts, na.rm = T)) %>%
  inner_join(total_counts, by = c("condition", "replicate")) %>%
  mutate(perc_counts = (read_length_counts / total_counts) * 100) %>%
  ggplot(aes(x = read_length, y = perc_counts, colour = condition, lty = replicate, shape = replicate))+
  geom_line(size = 1)+
  geom_point(size = 2)+
  myTheme+
  ylab("% counts")+
  xlab("read length")+
  scale_x_continuous(breaks = seq(min(lengths), max(lengths),2)) -> perc_counts_plot

png(filename = file.path(parent_dir, "plots/summed_counts/lengths_perc_plot.png"), width = 500, height = 300)
print(perc_counts_plot)
dev.off()

#for each sample, plot summed counts within each region across the length distribution
for (sample in RPF_sample_names) {
  all_data[all_data$sample == sample,] %>%
    mutate(region = factor(region, levels = c("UTR3", "UTR5", "CDS"), labels = c("3\'UTR", "5\'UTR", "CDS"), ordered = T)) %>%
    group_by(read_length, region) %>%
    summarise(region_counts = sum(counts, na.rm = T)) %>%
    ggplot(aes(x = read_length, y = region_counts, fill = region))+
    geom_col()+
    myTheme+
    xlab("read length")+
    ylab("total counts")+
    scale_x_continuous(breaks = seq(min(lengths), max(lengths),2))+
    ggtitle(str_replace_all(sample, "_", " ")) -> region_counts_plot
  
  png(filename = file.path(parent_dir, paste0("plots/summed_counts/", sample, "_summed_region_counts.png")), width = 500, height = 300)
  print(region_counts_plot)
  dev.off()
}

#for each sample, plot percentage counts within each region across the length distribution
all_data %>%
  group_by(sample, read_length) %>%
  summarise(total_counts = sum(counts, na.rm = T)) -> read_length_counts

for (sample in RPF_sample_names) {
  all_data[all_data$sample == sample,] %>%
    group_by(sample, region, read_length) %>%
    summarise(region_counts = sum(counts, na.rm = T)) %>%
    ungroup() %>%
    inner_join(read_length_counts, by = c("sample", "read_length")) %>%
    mutate(perc_counts = (region_counts / total_counts * 100),
           region = factor(region, levels = c("UTR3", "UTR5", "CDS"), labels = c("3\'UTR", "5\'UTR", "CDS"), ordered = T)) %>%
    #summary()
    ggplot(aes(x = read_length, y = perc_counts, fill = region))+
    geom_col()+
    myTheme+
    xlab("read length")+
    ylab("% counts")+
    scale_x_continuous(breaks = seq(min(lengths), max(lengths),2))+
    ggtitle(str_replace_all(sample, "_", " ")) -> region_perc_plot
  
  png(filename = file.path(parent_dir, paste0("plots/summed_counts/", sample, "_perc_region_counts.png")), width = 500, height = 300)
  print(region_perc_plot)
  dev.off()
}

#plot region percentages for all samples in one plot
all_data %>%
  group_by(sample, region, read_length) %>%
  summarise(region_counts = sum(counts, na.rm = T)) %>%
  ungroup() %>%
  inner_join(read_length_counts, by = c("sample", "read_length")) %>%
  mutate(perc_counts = (region_counts / total_counts * 100),
         region = factor(region, levels = c("UTR3", "UTR5", "CDS"), labels = c("3\'UTR", "5\'UTR", "CDS"), ordered = T)) %>%
  ggplot(aes(x = read_length, y = perc_counts, colour = region))+
  geom_point()+
  #geom_boxplot(aes(fill = region))+
  myTheme+
  xlab("read length")+
  ylab("% counts")+
  scale_x_continuous(breaks = seq(min(lengths), max(lengths),2)) -> region_perc_plot

png(filename = file.path(parent_dir, "plots/summed_counts/all_samples_perc_region_counts.png"), width = 500, height = 300)
print(region_perc_plot)
dev.off()
