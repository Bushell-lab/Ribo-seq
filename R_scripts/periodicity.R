#load libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(parallel)

#read in common variables
source("common_variables.R")
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
    fyle_list[[paste(sample, i, sep = "_")]] <- file.path(parent_dir, "Analysis/periodicity", paste0(sample, "_pc_L", i, "_Off0_periodicity.csv"))
  }
}

#read in the data with parLapply
no_cores <- detectCores() - 1 #sets the number of cores to use (all but one)
cl <- makeCluster(no_cores) #Initiates cluster
data_list <- parLapply(cl, fyle_list, read_counts_csv) #reads in the data
stopCluster(cl) #Stops cluster

#combine data_list into one data frame
#extract sample and read length from fyle
do.call("rbind", data_list) %>%
  mutate(read_length = str_remove(fyle, ".+pc_L"),
         read_length = as.numeric(str_remove(read_length, "_Off0_.+")),
         sample = str_remove(fyle, ".+periodicity/"),
         sample = factor(str_remove(sample, "_pc.+"))) %>%
  select(-fyle) -> all_data

summary(all_data)

#plot data----
#individual samples
for (sample in RPF_sample_names) {
  df <- all_data[all_data$sample == sample,]
  
  df %>%
    group_by(read_length) %>%
    summarise(total_counts = sum(counts)) -> summed_counts
  
  df %>%
    group_by(read_length, frame) %>%
    summarise(frame_counts = sum(counts)) %>%
    inner_join(summed_counts, by = "read_length") %>%
    mutate(frame_perc = (frame_counts / total_counts) * 100,
           frame = factor(frame, levels = c("f2", "f1", "f0"), ordered = T)) %>%
    ggplot(aes(x = read_length, y = frame_perc, fill = frame))+
    geom_col()+
    scale_x_continuous(breaks = seq(min(lengths), max(lengths),2))+
    ylab("% counts")+
    xlab("read length")+
    ggtitle(sample)+
    myTheme -> periodicity_col_plot
  
  png(filename = file.path(parent_dir, paste0("plots/periodicity/", sample, "_periodicity.png")), width = 500, height = 300)
  print(periodicity_col_plot)
  dev.off()
  
  plot_list <- list()
  for (i in lengths) {
    df[df$read_length == i,] %>%
      group_by(read_length, transcript) %>%
      summarise(total_counts = sum(counts)) %>%
      filter(total_counts > 50) -> transcript_counts
    
    df[df$read_length == i,] %>%
      inner_join(transcript_counts, by = c("transcript", "read_length")) %>%
      mutate(frame_perc = (counts / total_counts) * 100) %>%
      ggplot(aes(x = frame, y = frame_perc, fill = frame))+
      geom_boxplot()+
      xlab("Frame")+
      ylab("% counts")+
      myTheme+
      theme(legend.position = "none")+
      ggtitle(paste("read length", i)) -> plot_list[[i]]
  }
  png(filename = file.path(parent_dir, paste0("plots/periodicity/", sample, "_periodicity_boxplots.png")), width = 1000, height = 500)
  grid.arrange(plot_list[[28]], plot_list[[29]], plot_list[[30]],
               plot_list[[31]], plot_list[[32]], plot_list[[33]], nrow = 2)
  dev.off()
}

#all samples
all_data %>%
  group_by(read_length, sample) %>%
  summarise(total_counts = sum(counts)) -> summed_counts

all_data %>%
  group_by(read_length, frame, sample) %>%
  summarise(frame_counts = sum(counts)) %>%
  inner_join(summed_counts, by = c("read_length", "sample")) %>%
  mutate(frame_perc = (frame_counts / total_counts) * 100,
         frame = factor(frame, levels = c("f2", "f1", "f0"), ordered = T)) %>%
  ggplot(aes(x = factor(read_length), y = frame_perc, fill = frame))+
  geom_boxplot(width = 0.5, outlier.shape=NA)+
  ylab("% counts")+
  xlab("read length")+
  myTheme -> periodicity_col_plot

png(filename = file.path(parent_dir, paste0("plots/periodicity/all_samples_periodicity.png")), width = 500, height = 300)
print(periodicity_col_plot)
dev.off()
