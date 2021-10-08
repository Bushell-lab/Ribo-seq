#load libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(parallel)
library(viridis)

#read in common variables
source("common_variables.R")

#functions
#write a function that will read in a csv file for use with parLapply
read_counts_csv <- function(k){
  df <- read.csv(file = k)
  df$fyle <- rep(k)
  return(df)
}

#exports just the legend of a plot
myLegend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#themes
myTheme <- theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.position = "none")

#read in data----
#generate a list of file names
fyle_list <- list()
for(sample in RPF_sample_names) {
  for(i in lengths){
    fyle_list[[paste(sample, i, "start", sep = "_")]] <- file.path(parent_dir, "Analysis/spliced_counts", paste0(sample, "_pc_best_L", i, "_Off0_start_site.csv"))
    fyle_list[[paste(sample, i, "stop", sep = "_")]] <- file.path(parent_dir, "Analysis/spliced_counts", paste0(sample, "_pc_best_L", i, "_Off0_stop_site.csv"))
    fyle_list[[paste(sample, i, "UTR5_start", sep = "_")]] <- file.path(parent_dir, "Analysis/spliced_counts", paste0(sample, "_pc_best_L", i, "_Off0_UTR5_start.csv"))
    fyle_list[[paste(sample, i, "UTR3_end", sep = "_")]] <- file.path(parent_dir, "Analysis/spliced_counts", paste0(sample, "_pc_best_L", i, "_Off0_UTR3_end.csv"))
  }
}

#read in the data with parLapply
no_cores <- detectCores() - 1 #sets the number of cores to use (all but one)
cl <- makeCluster(no_cores) #Initiates cluster
data_list <- parLapply(cl, fyle_list, read_counts_csv) #reads in the data
stopCluster(cl) #Stops cluster

#combine data_list into one data frame
all_data <- do.call("rbind", data_list)

#extract sample, read length and splice position from fylenames
all_data %>%
  mutate(read_length = str_remove(fyle, ".+pc_best_L"),
         read_length = as.numeric(str_remove(read_length, "_Off0_.+")),
         sample = str_remove(fyle, ".+spliced_counts/"),
         sample = str_remove(sample, "_pc_best.+"),
         splice = str_remove(fyle, ".+Off0_"),
         splice = str_remove(splice, ".csv")) %>%
  select(-fyle) -> all_data

summary(all_data)

#plot heatmaps----
for (sample in RPF_sample_names) {
  
  fill_lims <- c(0, max(all_data$counts[all_data$sample == sample]))
  
  all_data[all_data$sample == sample & all_data$splice == "start_site",] %>%
    ggplot(aes(x = position, y = read_length, fill= counts)) + 
    geom_tile()+
    scale_fill_viridis(discrete = F, limits = fill_lims)+
    xlab("Position relative to start codon")+
    ylab("Read length")+
    myTheme -> start_site_plot
  
  all_data[all_data$sample == sample & all_data$splice == "stop_site",] %>%
    ggplot(aes(x = position, y = read_length, fill= counts)) + 
    geom_tile()+
    scale_fill_viridis(discrete = F, limits = fill_lims)+
    xlab("Position relative to stop codon")+
    ylab("Read length")+
    myTheme -> stop_site_plot
  
  all_data[all_data$sample == sample & all_data$splice == "UTR5_start",] %>%
    ggplot(aes(x = position, y = read_length, fill= counts)) + 
    geom_tile()+
    scale_fill_viridis(discrete = F, limits = fill_lims)+
    xlab("Position relative to 5\' end of transcript")+
    ylab("Read length")+
    myTheme -> UTR5_start_plot
  
  all_data[all_data$sample == sample & all_data$splice == "UTR3_end",] %>%
    mutate(position = position - 50) %>%
    ggplot(aes(x = position, y = read_length, fill= counts)) + 
    geom_tile()+
    scale_fill_viridis(discrete = F, limits = fill_lims)+
    xlab("Position relative to 3\' end of transcript")+
    ylab("Read length")+
    myTheme -> UTR3_end_plot
  
  all_data %>%
    ggplot(aes(x = position, y = read_length, fill= counts)) + 
    geom_tile()+
    scale_fill_viridis(discrete = F, limits = fill_lims)+
    theme(legend.text = element_text(size = 16),
          legend.title = element_text(size = 16)) -> legend_plot
  
  legend <- myLegend(legend_plot)
  
  png(filename = file.path(parent_dir, paste0("plots/heatmaps/", sample, "_heatmap.png")), width = 2000, height = 300)
  grid.arrange(UTR5_start_plot, start_site_plot, stop_site_plot, UTR3_end_plot, legend,
               widths = c(1,2,2,1, 0.2))
  dev.off()
}
