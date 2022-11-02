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

#themes
myTheme <- theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5))

#read in data----
#generate a list of file names
fyle_list <- list()
for(sample in RPF_sample_names) {
  for(i in lengths){
    fyle_list[[paste(sample, i, "start", sep = "_")]] <- file.path(parent_dir, "Analysis/spliced_counts", paste0(sample, "_pc_L", i, "_Off0_start_site.csv"))
    fyle_list[[paste(sample, i, "stop", sep = "_")]] <- file.path(parent_dir, "Analysis/spliced_counts", paste0(sample, "_pc_L", i, "_Off0_stop_site.csv"))
  }
}

#read in the data with parLapply
no_cores <- detectCores() - 1 #sets the number of cores to use (all but one)
cl <- makeCluster(no_cores) #Initiates cluster
data_list <- parLapply(cl, fyle_list, read_counts_csv) #reads in the data
stopCluster(cl) #Stops cluster

#combine data_list into one data frame
#extract sample, read length and splice position from fyle
do.call("rbind", data_list) %>%
  mutate(read_length = str_remove(fyle, ".+pc_L"),
         read_length = as.numeric(str_remove(read_length, "_Off0_.+")),
         sample = str_remove(fyle, ".+spliced_counts/"),
         sample = factor(str_remove(sample, "_pc.+")),
         splice = str_remove(fyle, ".+Off0_"),
         splice = factor(str_remove(splice, ".csv"))) %>%
  select(-fyle) -> all_data

summary(all_data)

#plot offset plots----
for (sample in RPF_sample_names) {
  
  start_site_plot_list <- list()
  stop_site_plot_list <- list()
  
  for (i in lengths) {
    start_site_data <- all_data[all_data$sample == sample & all_data$splice == "start_site" & all_data$read_length == i,]
    stop_site_data <- all_data[all_data$sample == sample & all_data$splice == "stop_site" & all_data$read_length == i,]
    
    start_site_data %>%
      ggplot(aes(x = position, y = counts)) + 
      geom_col()+
      xlab("Position relative to start codon")+
      ylab("Total counts")+
      geom_vline(xintercept = -12, colour = "red", lty=2)+
      scale_x_continuous(limits = c(-25, 25), breaks = c(-25, -12, 0, 25))+
      myTheme+
      ggtitle(paste("read length", i)) -> start_site_plot_list[[i]]
    
    stop_site_data %>%
      ggplot(aes(x = position, y = counts)) + 
      geom_col()+
      xlab("Position relative to stop codon")+
      ylab("Total counts")+
      geom_vline(xintercept = -18, colour = "red", lty=2)+
      scale_x_continuous(limits = c(-25, 25), breaks = c(-25, -18, 0, 25))+
      myTheme+
      ggtitle(paste("read length", i)) -> stop_site_plot_list[[i]]
  }
  png(filename = file.path(parent_dir, paste0("plots/offset/", sample, "_start_site_offset.png")), width = 1000, height = 500)
  grid.arrange(start_site_plot_list[[28]], start_site_plot_list[[29]], start_site_plot_list[[30]],
               start_site_plot_list[[31]], start_site_plot_list[[32]], start_site_plot_list[[33]], nrow = 2)
  dev.off()
  
  png(filename = file.path(parent_dir, paste0("plots/offset/", sample, "_stop_site_offset.png")), width = 1000, height = 500)
  grid.arrange(stop_site_plot_list[[28]], stop_site_plot_list[[29]], stop_site_plot_list[[30]],
               stop_site_plot_list[[31]], stop_site_plot_list[[32]], stop_site_plot_list[[33]], nrow = 2)
  dev.off()
}

#all samples
start_site_plot_list <- list()
stop_site_plot_list <- list()

for (i in lengths) {
  start_site_data <- all_data[all_data$splice == "start_site" & all_data$read_length == i,]
  stop_site_data <- all_data[all_data$splice == "stop_site" & all_data$read_length == i,]
  
  start_site_data %>%
    ggplot(aes(x = position, y = counts)) + 
    geom_col()+
    xlab("Position relative to start codon")+
    ylab("Total counts")+
    geom_vline(xintercept = -12, colour = "red", lty=2)+
    scale_x_continuous(limits = c(-25, 25), breaks = c(-25, -12, 0, 25))+
    myTheme+
    ggtitle(paste("read length", i)) -> start_site_plot_list[[i]]
  
  stop_site_data %>%
    ggplot(aes(x = position, y = counts)) + 
    geom_col()+
    xlab("Position relative to stop codon")+
    ylab("Total counts")+
    geom_vline(xintercept = -18, colour = "red", lty=2)+
    scale_x_continuous(limits = c(-25, 25), breaks = c(-25, -18, 0, 25))+
    myTheme+
    ggtitle(paste("read length", i)) -> stop_site_plot_list[[i]]
}
png(filename = file.path(parent_dir, paste0("plots/offset/all_samples_start_site_offset.png")), width = 1000, height = 500)
grid.arrange(start_site_plot_list[[29]], start_site_plot_list[[30]], start_site_plot_list[[31]],
             start_site_plot_list[[32]], start_site_plot_list[[33]], start_site_plot_list[[34]], nrow = 2)
dev.off()

png(filename = file.path(parent_dir, paste0("plots/offset/all_samples_stop_site_offset.png")), width = 1000, height = 500)
grid.arrange(stop_site_plot_list[[29]], stop_site_plot_list[[30]], stop_site_plot_list[[31]],
             stop_site_plot_list[[32]], stop_site_plot_list[[33]], stop_site_plot_list[[34]], nrow = 2)
dev.off()
