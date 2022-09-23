#functions----
#reads in csv files  (for use with parLapply)
read_counts_csv <- function(k){
  df <- read.csv(file = k, header = T)
  df$transcript <- rep(gsub("_.+", "", k))
  return(df)
}

#Returns transcript ID (for use with parLapply)
get_transcript_ID <- function(x) {
  transcript <- str_replace(x, "_.+", "")
  return(transcript)
}

normalise_data <- function(df, tpms) {
  
  df %>%
    group_by(condition, replicate) %>%
    summarise(total_counts = sum(Counts)) -> total_counts
  
  df %>%
    inner_join(total_counts, by = c("condition", "replicate")) %>%
    inner_join(tpms, by = c("transcript", "condition", "replicate")) %>%
    filter(tpm > 0) %>%
    mutate(CPM = Counts / total_counts * 1000000, #calculates counts per million (CPM)
           normalised_CPM = CPM / tpm) %>% #normalise to total RNA-seq transcripts per million (TPM)
    select(-c("total_counts", "tpm")) -> normalised_df
  
  return(normalised_df)
}

calculate_bins <- function(x, n) {
  bins <- n
  cut_size <- 1 / bins
  breaks <- seq(0, 1, cut_size)
  start <- cut_size / 2
  stop <- 1 - start
  bin_centres <- seq(start, stop, cut_size)
  bin <- as.numeric(cut(x, breaks, include.lowest = F, labels = bin_centres))
  return(bin)
}

#the following function will bin the data (df) within the 5'UTR, CDS and 3'UTR.
#It needs a region lengths file
#region cutoffs are the minimum sizes of the regions to include. It should be supplied as a numeric vector, referring to the minimum size for the 5'UTR, CDS and 3'UTR (in that order). Default is c(100,300,100)
#bins are the number of bins within each region. This should also be a numeric vector, referring to the number of bins for the 5'UTR, CDS and 3'UTR (in that order). Default is c(25,50,25)
bin_data <- function(df, region_lengths, region_cutoffs = c(50,300,50), bins = c(25,50,25)) {
  
  #filter the transcripts with UTRs/CDSs less than the defined thresholds
  df %>%
    inner_join(region_lengths, by = "transcript") %>%
    filter(UTR5_len >= region_cutoffs[1] & CDS_len >= region_cutoffs[2] & UTR5_len >= region_cutoffs[3]) -> filtered_data
  
  #bin data
  #5'UTR
  filtered_data %>%
    filter(Position <= UTR5_len) %>%
    mutate(normalised_position = Position / UTR5_len,
           bin = calculate_bins(normalised_position, bins[1])) %>%
    group_by(transcript, bin, condition, replicate) %>%
    summarise(binned_cpm = mean(CPM),
              binned_normalised_cpm = mean(normalised_CPM)) %>%
    ungroup() %>%
    mutate(region = factor(rep("UTR5"))) -> UTR5_df
  
  #CDS
  filtered_data %>%
    filter(Position > UTR5_len & Position <= (UTR5_len + CDS_len)) %>%
    mutate(CDS_position = Position - UTR5_len,
           normalised_position = CDS_position / CDS_len,
           bin = calculate_bins(normalised_position, bins[2])) %>%
    group_by(transcript, bin, condition, replicate) %>%
    summarise(binned_cpm = mean(CPM),
              binned_normalised_cpm = mean(normalised_CPM)) %>%
    ungroup() %>%
    mutate(region = factor(rep("CDS"))) -> CDS_df
  
  #3'UTR
  filtered_data %>%
    filter(Position > (UTR5_len + CDS_len)) %>%
    mutate(UTR3_position = Position - (UTR5_len + CDS_len),
           normalised_position = UTR3_position / UTR3_len,
           bin = calculate_bins(normalised_position, bins[3])) %>%
    group_by(transcript, bin, condition, replicate) %>%
    summarise(binned_cpm = mean(CPM),
              binned_normalised_cpm = mean(normalised_CPM)) %>%
    ungroup() %>%
    mutate(region = factor(rep("UTR3"))) -> UTR3_df  
  
  return(bind_rows(UTR5_df, CDS_df, UTR3_df))
}

calculate_binned_delta <- function(alist, value, control = control, treatment = treatment, paired_data = T) {
  do.call("rbind", alist) %>%
    rename(value = value) %>%
    select(transcript, bin, replicate, region, condition, value) %>%
    filter(condition == control | condition == treatment) %>%
    spread(key = condition, value = value) %>%
    rename(control = control,
           treatment = treatment) %>%
    group_by(bin, replicate, region) %>%
    summarise(mean_ctrl = mean(control, na.rm = T),
              mean_treatment = mean(treatment, na.rm = T)) %>%
    ungroup() -> spread_data
  
  #calculate delta with 95% confidence intervals
  delta_list <- list()
  for (region in c("UTR5", "CDS", "UTR3")) {
    df <- spread_data[spread_data$region == region,]
    bins <- max(df$bin)
    for (n in 1:bins) {
      bin_n <- as.data.frame(df[df$bin == n,])
      
      if (paired_data == T) {
        t <- t.test(bin_n[,"mean_treatment"], bin_n[,"mean_ctrl"], paired = T, conf.int = T)
        
        delta_list[[paste(region, n, sep = "_")]] <- data.frame(region = region,
                                                                bin = n,
                                                                upper = t$conf.int[[1]],
                                                                lower = t$conf.int[[2]],
                                                                delta = t$estimate)
      } else {
        t <- t.test(bin_n[,"mean_treatment"], bin_n[,"mean_ctrl"], paired = F, conf.int = T)
        
        delta_list[[paste(region, n, sep = "_")]] <- data.frame(region = region,
                                                                bin = n,
                                                                upper = t$conf.int[[1]],
                                                                lower = t$conf.int[[2]],
                                                                delta = (t$estimate[1] - t$estimate[2]))
      }
      
    }
  }
  delta_data <- do.call("rbind", delta_list)
  
  #convert 95% confidence intervals to 0 if delta is exactly 0 (will be NaN)
  delta_data$upper[delta_data$delta == 0] <- 0
  delta_data$lower[delta_data$delta == 0] <- 0
  
  return(delta_data)
}

calculate_positional_counts <- function(df, region = "CDS") {
  df[df$region == region,] %>%
    group_by(transcript, region) %>%
    summarise(total_counts = sum(binned_cpm)) %>%
    ungroup() -> summed_counts
  
  df[df$region == region,] %>%
    inner_join(summed_counts, by = c("transcript", "region")) %>%
    filter(total_counts > 0) %>%
    mutate(positional_counts = binned_cpm / total_counts) -> positional_counts_df
  
  return(positional_counts_df)
}

calculate_positional_delta <- function(alist, control = control, treatment = treatment, paired_data = T) {
  do.call("rbind", alist) %>%
    spread(key = condition, value = positional_counts) %>%
    group_by(bin, replicate) %>%
    rename(control = control,
           treatment = treatment) %>%
    summarise(mean_ctrl = mean(control, na.rm = T),
              mean_treatment = mean(treatment, na.rm = T)) %>%
    ungroup() -> spread_data
  
  #calculate delta with 95% confidence intervals
  delta_list <- list()
  
  bins <- max(spread_data$bin)
  for (n in 1:bins) {
    bin_n <- as.data.frame(spread_data[spread_data$bin == n,])
    
    if (paired_data == T) {
      t <- t.test(bin_n[,"mean_treatment"], bin_n[,"mean_ctrl"], paired = T, conf.int = T)
      
      delta_list[[n]] <- data.frame(bin = n,
                                    upper = t$conf.int[[1]],
                                    lower = t$conf.int[[2]],
                                    delta = t$estimate)
    } else {
      t <- t.test(bin_n[,"mean_treatment"], bin_n[,"mean_ctrl"], paired = F, conf.int = T)
      
      delta_list[[n]] <- data.frame(bin = n,
                                    upper = t$conf.int[[1]],
                                    lower = t$conf.int[[2]],
                                    delta = (t$estimate[1] - t$estimate[2]))
    }
    
  }
  delta_data <- do.call("rbind", delta_list)
  
  #convert 95% confidence intervals to 0 if delta is exactly 0 (will be NaN)
  delta_data$upper[delta_data$delta == 0] <- 0
  delta_data$lower[delta_data$delta == 0] <- 0
  
  return(delta_data)
}

splice_single_nt <- function(df, region_lengths, codons = 50, UTR_nts = 48, region_cutoffs = c(50,300,50)) {
  nt <- codons * 3
  
  df %>%
    inner_join(region_lengths, by = "transcript") %>%
    filter(UTR5_len >= region_cutoffs[1] & CDS_len >= region_cutoffs[2] & UTR3_len >= region_cutoffs[3]) %>%
    mutate(region = factor(case_when(Position <= UTR5_len ~ "UTR5",
                                     Position > UTR5_len & Position <= (UTR5_len + CDS_len) ~ "CDS",
                                     Position > (UTR5_len + CDS_len) ~ "UTR3"))) -> filtered_data
  
  #UTR5
  filtered_data[filtered_data$region == "UTR5",] %>%
    filter(Position <= UTR_nts) %>%
    arrange(transcript, Position) %>% 
    mutate(window = rep(seq(2, (UTR_nts - 1), by = 3), each = 3, times = n_distinct(transcript))) %>%
    group_by(transcript, condition, replicate, region, window) %>%
    summarise(single_nt_normalised_cpm = mean(normalised_CPM),
              single_nt_cpm = mean(CPM)) %>%
    ungroup() -> UTR5_start
  
  filtered_data[filtered_data$region == "UTR5",] %>%
    group_by(transcript) %>%
    top_n(n = UTR_nts, wt = Position) %>% #extracts the 3' most nts
    arrange(transcript, Position) %>%
    mutate(window = rep(seq((-UTR_nts + 1), -2, by = 3), each = 3, times = n_distinct(transcript))) %>%
    ungroup() %>%
    group_by(transcript, condition, replicate, region, window) %>%
    summarise(single_nt_normalised_cpm = mean(normalised_CPM),
              single_nt_cpm = mean(CPM)) %>%
    ungroup() -> UTR5_end
  
  #CDS
  filtered_data[filtered_data$region == "CDS",] %>%
    mutate(CDS_position = Position - UTR5_len) %>%
    filter(CDS_position <= (codons * 3)) %>%
    arrange(transcript, CDS_position) %>% 
    mutate(window = rep(seq(2, ((codons * 3) - 1), by = 3), each = 3, times = n_distinct(transcript))) %>%
    group_by(transcript, condition, replicate, region, window) %>%
    summarise(single_nt_normalised_cpm = mean(normalised_CPM),
              single_nt_cpm = mean(CPM)) %>%
    ungroup() -> CDS_start
  
  filtered_data[filtered_data$region == "CDS",] %>%
    group_by(transcript) %>%
    top_n(n = (codons * 3), wt = Position) %>% #extracts the 3' most nts
    arrange(transcript, Position) %>%
    mutate(window = rep(seq((-(codons * 3) + 1), -2, by = 3), each = 3, times = n_distinct(transcript))) %>%
    ungroup() %>%
    group_by(transcript, condition, replicate, region, window) %>%
    summarise(single_nt_normalised_cpm = mean(normalised_CPM),
              single_nt_cpm = mean(CPM)) %>%
    ungroup() -> CDS_end
  
  #UTR3
  filtered_data[filtered_data$region == "UTR3",] %>%
    mutate(UTR3_position = Position - (UTR5_len + CDS_len)) %>%
    filter(UTR3_position <= UTR_nts) %>%
    arrange(transcript, UTR3_position) %>% 
    mutate(window = rep(seq(2, (UTR_nts - 1), by = 3), each = 3, times = n_distinct(transcript))) %>%
    group_by(transcript, condition, replicate, region, window) %>%
    summarise(single_nt_normalised_cpm = mean(normalised_CPM),
              single_nt_cpm = mean(CPM)) %>%
    ungroup() -> UTR3_start
  
  filtered_data[filtered_data$region == "UTR3",] %>%
    group_by(transcript) %>%
    top_n(n = UTR_nts, wt = Position) %>% #extracts the 3' most nts
    arrange(transcript, Position) %>%
    mutate(window = rep(seq((-UTR_nts + 1), -2, by = 3), each = 3, times = n_distinct(transcript))) %>%
    ungroup() %>%
    group_by(transcript, condition, replicate, region, window) %>%
    summarise(single_nt_normalised_cpm = mean(normalised_CPM),
              single_nt_cpm = mean(CPM)) %>%
    ungroup() -> UTR3_end
  
  return(bind_rows(UTR5_start, UTR5_end, CDS_start, CDS_end, UTR3_start, UTR3_end))
}

calculate_single_nt_delta <- function(alist, value, control = control, treatment = treatment, paired_data = T) {
  do.call("rbind", alist) %>%
    rename(value = value) %>%
    select(transcript, window, replicate, region, condition, value) %>%
    spread(key = condition, value = value) %>%
    group_by(window, replicate, region) %>%
    rename(control = control,
           treatment = treatment) %>%
    summarise(mean_ctrl = mean(control, na.rm = T),
              mean_treatment = mean(treatment, na.rm = T)) %>%
    ungroup() -> spread_data
  
  #calculate delta with 95% confidence intervals
  delta_list <- list()
  for (region in c("UTR5", "CDS", "UTR3")) {
    df <- spread_data[spread_data$region == region,]
    max_window <- max(df$window)
    for (n in c(seq(-max_window, -2, 3),seq(2, max_window, 3))) {
      df_n <- as.data.frame(df[df$window == n,])
      
      if (paired_data == T) {
        t <- t.test(df_n[,"mean_treatment"], df_n[,"mean_ctrl"], paired = T, conf.int = T)
        
        delta_list[[paste(region, n, sep = "_")]] <- data.frame(region = region,
                                                                window = n,
                                                                upper = t$conf.int[[1]],
                                                                lower = t$conf.int[[2]],
                                                                delta = t$estimate)
      } else {
        t <- t.test(df_n[,"mean_treatment"], df_n[,"mean_ctrl"], paired = F, conf.int = T)
        
        delta_list[[paste(region, n, sep = "_")]] <- data.frame(region = region,
                                                                window = n,
                                                                upper = t$conf.int[[1]],
                                                                lower = t$conf.int[[2]],
                                                                delta = (t$estimate[1] - t$estimate[2]))
      }
      
    }
  }
  delta_data <- do.call("rbind", delta_list)
  
  #convert 95% confidence intervals to 0 if delta is exactly 0 (will be NaN)
  delta_data$upper[delta_data$delta == 0] <- 0
  delta_data$lower[delta_data$delta == 0] <- 0
  
  return(delta_data)
}

#the following function will summarise (mean and median) across all transcripts within each sample for the given value.
#Grouping needs to be set to either "bin" (for binned_data) or "window" (for single_nt data)
summarise_data <- function(df, value, grouping) {
  df %>%
    rename(value = value,
           grouping = grouping) %>%
    group_by(region, condition, replicate, grouping) %>%
    summarise(mean_counts = mean(value),
              median_counts = median(value)) %>%
    ungroup() -> summarised_df
  
  return(summarised_df)
}

#plotting----
#the following functions will plot the binned/positional/single_nt raw data as lines
#Setting SD to True will add shaded areas to represent standard deviation
#the value to be plotted needs to be called "average_counts" within the dataframe supplied (df)

plot_binned_lines <- function(df, SD = T, control = control, treatment = treatment) {
  
  #order condition
  df$condition <- factor(df$condition, levels = c(control, treatment), ordered = T)
  
  #calculate axis limits
  if (SD == F) {
    ylims <- c(0,max(df$average_counts))
  } else {
    ylims <- c(min(df$average_counts - df$sd_counts),
               max(df$average_counts + df$sd_counts))
  }
  
  #5'UTR
  df[df$region == "UTR5",] %>%
    ggplot(aes(x = grouping, y = average_counts, colour = condition))+
    geom_line(size = 1)+
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    ylim(ylims)+
    UTR5_theme -> UTR5_plot
  
  #CDS
  df[df$region == "CDS",] %>%
    ggplot(aes(x = grouping, y = average_counts, colour = condition))+
    geom_line(size = 1)+
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    ylim(ylims)+
    CDS_theme -> CDS_plot
  
  #3'UTR
  df[df$region == "UTR3",] %>%
    ggplot(aes(x = grouping, y = average_counts, colour = condition))+
    geom_line(size = 1)+
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    ylim(ylims)+
    UTR3_theme -> UTR3_plot
  
  return(list(UTR5_plot, CDS_plot, UTR3_plot))
}

plot_binned_all_replicates <- function(alist, control = control, treatment = treatment) {
  
  df <- do.call("rbind", alist)
  
  #order condition
  df$condition <- factor(df$condition, levels = c(control, treatment), ordered = T)
  
  #calculate axis limits
  ylims <- c(0,max(df$mean_counts))
  
  #5'UTR
  df[df$region == "UTR5",] %>%
    ggplot(aes(x = grouping, y = mean_counts, colour = condition, lty = factor(replicate)))+
    geom_line(size = 1)+
    ylim(ylims)+
    UTR5_theme -> UTR5_plot
  
  #CDS
  df[df$region == "CDS",] %>%
    ggplot(aes(x = grouping, y = mean_counts, colour = condition, lty = factor(replicate)))+
    geom_line(size = 1)+
    ylim(ylims)+
    CDS_theme -> CDS_plot
  
  #3'UTR
  df[df$region == "UTR3",] %>%
    ggplot(aes(x = grouping, y = mean_counts, colour = condition, lty = factor(replicate)))+
    geom_line(size = 1)+
    ylim(ylims)+
    UTR3_theme -> UTR3_plot
  
  return(list(UTR5_plot, CDS_plot, UTR3_plot))
}

plot_positional_lines <- function(df, SD = T, control = control, treatment = treatment) {
  
  #order condition
  df$condition <- factor(df$condition, levels = c(control, treatment), ordered = T)
  
  #calculate axis limits
  if (SD == F) {
    ylims <- c(0,max(df$average_counts))
  } else {
    ylims <- c(min(df$average_counts - df$sd_counts),
               max(df$average_counts + df$sd_counts))
  }
  
  df %>%
    ggplot(aes(x = grouping, y = average_counts, colour = condition))+
    geom_line(size = 1)+
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    ylim(ylims)+
    UTR3_theme+
    theme(axis.text.y = element_text(size = 18)) -> positional_plot
  
  return(positional_plot)
}

plot_single_nt_lines <- function(df, SD = T, plot_ends = F, control = control, treatment = treatment) {
  
  #order condition
  df$condition <- factor(df$condition, levels = c(control, treatment), ordered = T)
  
  #calculate axis limits
  if (SD == F) {
    ylims <- c(0,max(df$average_counts))
  } else {
    ylims <- c(min(df$average_counts - df$sd_counts),
               max(df$average_counts + df$sd_counts))
  }
  
  #5'UTR
  df[df$region == "UTR5" & df$grouping > 0,] %>%
    ggplot(aes(x = grouping, y = average_counts, colour = condition))+
    geom_line(size = 1)+
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    ylim(ylims)+
    UTR5_theme+
    xlab("nt\n(relative to 5' end)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16)) -> UTR5_start_plot
  
  df[df$region == "UTR5" & df$grouping < 0,] %>%
    ggplot(aes(x = grouping, y = average_counts, colour = condition))+
    geom_line(size = 1)+
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    ylim(ylims)+
    UTR5_theme+
    {if(plot_ends)CDS_theme}+
    xlab("nt\n(relative to start codon)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16)) -> UTR5_end_plot
  
  #CDS
  df[df$region == "CDS" & df$grouping > 0,] %>%
    ggplot(aes(x = grouping, y = average_counts, colour = condition))+
    geom_line(size = 1)+
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    ylim(ylims)+
    CDS_theme+
    xlab("nt\n(relative to start codon)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16)) -> CDS_start_plot
  
  df[df$region == "CDS" & df$grouping < 0,] %>%
    ggplot(aes(x = grouping, y = average_counts, colour = condition))+
    geom_line(size = 1)+
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    ylim(ylims)+
    CDS_theme+
    xlab("nt\n(relative to stop codon)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16)) -> CDS_end_plot
  
  #3'UTR
  df[df$region == "UTR3" & df$grouping > 0,] %>%
    ggplot(aes(x = grouping, y = average_counts, colour = condition))+
    geom_line(size = 1)+
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    ylim(ylims)+
    UTR3_theme+
    {if(plot_ends)CDS_theme}+
    xlab("nt\n(relative to stop codon)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16)) -> UTR3_start_plot
  
  df[df$region == "UTR3" & df$grouping < 0,] %>%
    ggplot(aes(x = grouping, y = average_counts, colour = condition))+
    geom_line(size = 1)+
    {if(SD)geom_ribbon(aes(ymin = average_counts-sd_counts, ymax = average_counts+sd_counts, fill = condition), alpha = 0.3, colour = NA)}+
    ylim(ylims)+
    UTR3_theme+
    xlab("nt\n(relative to 3' end)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16)) -> UTR3_end_plot
  
  if (plot_ends == T) {
    return(list(UTR5_start_plot, UTR5_end_plot, CDS_start_plot, CDS_end_plot, UTR3_start_plot, UTR3_end_plot))
  } else {
    return(list(UTR5_end_plot, CDS_start_plot, CDS_end_plot, UTR3_start_plot))
  }
}

#the following functions will plot the binned/positional/single_nt delta
#Setting SD to True will add shaded areas to represent 95% confidence intervals of the delta
plot_binned_delta <- function(df, SD = T) {
  
  #calculate axis limits
  lower_delta_ylim <- min(c(df$delta, df$upper))
  upper_delta_ylim <- max(c(df$delta,df$lower))
  ylims <- c(lower_delta_ylim, upper_delta_ylim)
  
  #5'UTR
  df[df$region == "UTR5",] %>%
    ggplot(aes(x = bin, y = delta))+
    geom_col(fill = "grey")+
    ylim(ylims)+
    UTR5_theme+
    {if(SD)geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, colour = NA)} -> UTR5_delta_plot
  
  #CDS
  df[df$region == "CDS",] %>%
    ggplot(aes(x = bin, y = delta))+
    geom_col(fill = "grey")+
    ylim(ylims)+
    CDS_theme+
    {if(SD)geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, colour = NA)} -> CDS_delta_plot
  
  #3'UTR
  df[df$region == "UTR3",] %>%
    ggplot(aes(x = bin, y = delta))+
    geom_col(fill = "grey")+
    ylim(ylims)+
    UTR3_theme+
    {if(SD)geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, colour = NA)} -> UTR3_delta_plot
  
  return(list(UTR5_delta_plot, CDS_delta_plot, UTR3_delta_plot))
}

plot_single_nt_delta <- function(df, SD = T, plot_ends = F) {
  
  #calculate axis limits
  lower_delta_ylim <- min(c(df$delta, df$upper))
  upper_delta_ylim <- max(c(df$delta,df$lower))
  ylims <- c(lower_delta_ylim, upper_delta_ylim)
  
  #5'UTR
  df[df$region == "UTR5" & df$window > 0,] %>%
    ggplot(aes(x = window, y = delta))+
    geom_col(fill = "grey")+
    ylim(ylims)+
    UTR5_theme+
    xlab("nt\n(relative to 5' end)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16))+
    {if(SD)geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, colour = NA)} -> UTR5_start_plot
  
  df[df$region == "UTR5" & df$window < 0,] %>%
    ggplot(aes(x = window, y = delta))+
    geom_col(fill = "grey")+
    ylim(ylims)+
    UTR5_theme+
    {if(plot_ends)CDS_theme}+
    xlab("nt\n(relative to start codon)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16))+
    {if(SD)geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, colour = NA)} -> UTR5_end_plot
  
  #CDS
  df[df$region == "CDS" & df$window > 0,] %>%
    ggplot(aes(x = window, y = delta))+
    geom_col(fill = "grey")+
    ylim(ylims)+
    CDS_theme+
    xlab("nt\n(relative to start codon)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16))+
    {if(SD)geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, colour = NA)} -> CDS_start_plot
  
  df[df$region == "CDS" & df$window < 0,] %>%
    ggplot(aes(x = window, y = delta))+
    geom_col(fill = "grey")+
    ylim(ylims)+
    CDS_theme+
    xlab("nt\n(relative to stop codon)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16))+
    {if(SD)geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, colour = NA)} -> CDS_end_plot
  
  #3'UTR
  df[df$region == "UTR3"  & df$window > 0,] %>%
    ggplot(aes(x = window, y = delta))+
    geom_col(fill = "grey")+
    ylim(ylims)+
    UTR3_theme+
    {if(plot_ends)CDS_theme}+
    xlab("nt\n(relative to stop codon)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16))+
    {if(SD)geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, colour = NA)} -> UTR3_start_plot
  
  df[df$region == "UTR3"  & df$window < 0,] %>%
    ggplot(aes(x = window, y = delta))+
    geom_col(fill = "grey")+
    ylim(ylims)+
    UTR3_theme+
    xlab("nt\n(relative to 3' end)")+
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16))+
    {if(SD)geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, colour = NA)} -> UTR3_stop_plot
  
  if (plot_ends == T) {
    return(list(UTR5_start_plot, UTR5_end_plot, CDS_start_plot, CDS_end_plot, UTR3_start_plot, UTR3_end_plot))
  } else {
    return(list(UTR5_end_plot, CDS_start_plot, CDS_end_plot, UTR3_start_plot))
  }
}

plot_positional_delta <- function(df, SD = T) {
  
  #calculate axis limits
  lower_delta_ylim <- min(c(df$delta, df$upper))
  upper_delta_ylim <- max(c(df$delta,df$lower))
  ylims <- c(lower_delta_ylim, upper_delta_ylim)
  
  df %>%
    ggplot(aes(x = bin, y = delta))+
    geom_col(fill = "grey")+
    ylim(ylims)+
    UTR3_theme+
    theme(axis.text.y = element_text(size = 18))+
    {if(SD)geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, colour = NA)} -> UTR5_delta_plot
  
  return(UTR5_delta_plot)
}

#filter data----
filter_outliers <- function(df, outliers) {
  return(df[!(df$transcript %in% outliers),])
}

#the following function will filter the data to only include the supplied transcript IDs
filter_transcripts <- function(df, transcript_IDs) {
  return(df[df$transcript %in% transcript_IDs,])
}

#plot subsets of data----
plot_subset <- function(IDs, subset, sub_dir,
                        binned_value, single_nt_value, control = control, treatment = treatment,
                        plot_binned = T, plot_single_nt = F, plot_positional = F, plot_delta = T,
                        SD = T, paired_data = T) {
  
  if (!(dir.exists(file.path(parent_dir, "plots/binned_plots", sub_dir)))) {
    dir.create(file.path(parent_dir, "plots/binned_plots", sub_dir))
  }
  
  if (plot_binned == T) {
    #binned
    
    #subset
    subset_binned_list <- lapply(binned_list, filter_transcripts, transcript_IDs = IDs)
    
    #summarise
    summarised_subset_binned_list <- lapply(subset_binned_list, summarise_data, value = binned_value, grouping = "bin")
    do.call("rbind", summarised_subset_binned_list) %>%
      group_by(grouping, condition, region) %>%
      summarise(average_counts = mean(mean_counts),
                sd_counts = sd(mean_counts)) %>%
      ungroup() -> summarised_subset_binned
    
    #plot lines
    subset_binned_line_plots <- plot_binned_lines(summarised_subset_binned, SD = SD, control = control, treatment = treatment)
    
    png(filename = file.path(parent_dir, "plots/binned_plots", sub_dir, paste(treatment, subset, binned_value, "lines.png")), width = 1000, height = 200)
    grid.arrange(subset_binned_line_plots[[1]], subset_binned_line_plots[[2]], subset_binned_line_plots[[3]], nrow = 1, widths = c(1,2,1.5))
    dev.off()
    
    if (plot_delta == T) {
      #calculate and plot delta
      subset_binned_delta_data <- calculate_binned_delta(subset_binned_list, value = binned_value, control = control, treatment = treatment, paired_data = paired_data)
      subset_binned_delta_plots <- plot_binned_delta(subset_binned_delta_data)
      
      png(filename = file.path(parent_dir, "plots/binned_plots", sub_dir, paste(treatment, subset, binned_value, "delta.png")), width = 1000, height = 200)
      grid.arrange(subset_binned_delta_plots[[1]], subset_binned_delta_plots[[2]], subset_binned_delta_plots[[3]], nrow = 1, widths = c(1,2,1))
      dev.off()
    }
  }
  
  if (plot_single_nt == T) {
    
    #single nt
    
    #subset
    subset_single_nt_list <- lapply(single_nt_list, filter_transcripts, transcript_IDs = IDs)
    
    #summarise
    summarised_subset_single_nt_list <- lapply(subset_single_nt_list, summarise_data, value = single_nt_value, grouping = "window")
    do.call("rbind", summarised_subset_single_nt_list) %>%
      group_by(grouping, condition, region) %>%
      summarise(average_counts = mean(mean_counts),
                sd_counts = sd(mean_counts)) %>%
      ungroup() -> summarised_subset_single_nt
    
    #plot lines
    subset_single_nt_line_plots <- plot_single_nt_lines(summarised_subset_single_nt, SD = SD, control = control, treatment = treatment)
    
    png(filename = file.path(parent_dir, "plots/binned_plots", sub_dir, paste(treatment, subset, single_nt_value, "lines.png")), width = 1300, height = 300)
    grid.arrange(subset_single_nt_line_plots[[1]], subset_single_nt_line_plots[[2]], subset_single_nt_line_plots[[3]], subset_single_nt_line_plots[[4]], nrow = 1, widths = c(1,2,2,1.5))
    dev.off()
    
    if (plot_delta == T) {
      #calculate and plot delta
      subset_single_nt_delta_data <- calculate_single_nt_delta(subset_single_nt_list, value = single_nt_value, control = control, treatment = treatment, paired_data = paired_data)
      subset_single_nt_delta_plots <- plot_single_nt_delta(subset_single_nt_delta_data)
      
      png(filename = file.path(parent_dir, "plots/binned_plots", sub_dir, paste(treatment, subset, single_nt_value, "delta.png")), width = 1300, height = 200)
      grid.arrange(subset_single_nt_delta_plots[[1]], subset_single_nt_delta_plots[[2]], subset_single_nt_delta_plots[[3]], subset_single_nt_delta_plots[[4]], nrow = 1, widths = c(1,2,2,1))
      dev.off()
    }
  }
  
  if (plot_positional == T) {
    #positional
    subset_positional_binned_list <- lapply(positional_list, filter_transcripts, transcript_IDs = IDs)
    
    summarised_positional_subset_binned_list <- lapply(subset_positional_binned_list, summarise_data, value = binned_value, grouping = "bin")
    do.call("rbind", summarised_positional_subset_binned_list) %>%
      group_by(grouping, condition) %>%
      summarise(average_counts = mean(mean_counts),
                sd_counts = sd(mean_counts)) %>%
      ungroup() -> summarised_positional_subset_binned
    
    subset_positional_line_plots <- plot_positional_lines(df = summarised_positional_subset_binned, SD = SD, control = control, treatment = treatment)
    
    png(filename = file.path(parent_dir, "plots/binned_plots", sub_dir, paste(treatment, subset, "binned positional lines.png")), width = 500, height = 200)
    print(subset_positional_line_plots)
    dev.off()
    
    if (plot_delta == T) {
      subset_binned_positional_delta <- calculate_positional_delta(subset_positional_binned_list, control = control, treatment = treatment, paired_data = paired_data)
      subset_positional_delta_plots <- plot_positional_delta(subset_binned_positional_delta)
      
      png(filename = file.path(parent_dir, "plots/binned_plots", sub_dir, paste(treatment, subset, "binned positional delta.png")), width = 500, height = 200)
      print(subset_positional_delta_plots)
      dev.off()
    }
  }
}


plot_GSEA_binned <- function(GSEA_set, pathway, subdir,
                             human = T, conversion_table = NULL,
                             binned_value = binned_value, single_nt_value = single_nt_value,
                             plot_binned = T, plot_single_nt = F, plot_positional = F,
                             sub_dir, control = control, treatment = treatment, paired_data = T, SD = T, plot_delta = T) {
  
  gene_list <- GSEA_set[[pathway]]
  
  if (human == T) {
    most_abundant_transcripts %>%
      filter(gene_sym %in% gene_list) %>%
      pull(transcript) -> GSEA_transcript_IDs
  } else {
    most_abundant_transcripts %>%
      inner_join(conversion_table, by = "gene") %>%
      filter(Human_gene_name %in% gene_list) %>%
      pull(transcript) -> GSEA_transcript_IDs
  }
  
  plot_subset(IDs = GSEA_transcript_IDs, subset = pathway, sub_dir= sub_dir,
                            binned_value = binned_value, single_nt_value = single_nt_value, control = control, treatment = treatment,
                            plot_binned = plot_binned, plot_single_nt = plot_single_nt, plot_positional = plot_positional, plot_delta = plot_delta,
                            SD = SD, paired_data = paired_data)
  
}

plot_single_transcripts <- function(gene, dir,
                                    plot_binned = T, plot_single_nt = F, plot_positional = F,
                                    SD = T, plot_replicates = T, plot_delta = T,
                                    control = control, treatment = treatment, paired_data = T,
                                    region_cutoffs = c(0,0,0)) {
  
  if (!(dir.exists(file.path(parent_dir, "plots/binned_plots/single_transcripts", dir)))) {
    dir.create(file.path(parent_dir, "plots/binned_plots/single_transcripts", dir))
  }
  
  transcript <- most_abundant_transcripts$transcript[most_abundant_transcripts$gene_sym == gene]
  
  #extract data from counts list
  filtered_counts_list <- lapply(counts_list, filter_transcripts, transcript_IDs = transcript)
  
  #bin data
  if (plot_binned == T) {
    binned_list <- lapply(filtered_counts_list, bin_data, region_lengths = region_lengths, region_cutoffs = region_cutoffs)
    
    summarised_binned_list <- lapply(binned_list, summarise_data, value = "binned_normalised_cpm", grouping = "bin")
    
    if (plot_replicates == F) {
      
      #plot average across replicates
      
      do.call("rbind", summarised_binned_list) %>%
        group_by(grouping, condition, region) %>%
        summarise(average_counts = mean(mean_counts),
                  sd_counts = sd(mean_counts)) %>%
        ungroup() -> summarised_binned
      
      binned_line_plots <- plot_binned_lines(df = summarised_binned, SD = SD, control = control, treatment = treatment)
      
      png(filename = file.path(parent_dir, "plots/binned_plots/single_transcripts", dir, paste(treatment, gene, "binned lines.png")), width = 1000, height = 200)
      grid.arrange(binned_line_plots[[1]], binned_line_plots[[2]], binned_line_plots[[3]], nrow = 1, widths = c(1,2,1.5))
      dev.off()
      
    } else {
      
      #plot individual replicates
      
      binned_line_plots_all_replicates <- plot_binned_all_replicates(summarised_binned_list, control = control, treatment = treatment)
      
      png(filename = file.path(parent_dir, "plots/binned_plots/single_transcripts", dir, paste(treatment, gene, "binned lines all replicates.png")), width = 1000, height = 200)
      grid.arrange(binned_line_plots_all_replicates[[1]], binned_line_plots_all_replicates[[2]], binned_line_plots_all_replicates[[3]], nrow = 1, widths = c(1,2,1.5))
      dev.off()
    }
    
    #plot delta
    if (plot_delta == T) {
      binned_delta_data <- calculate_binned_delta(binned_list, value = "binned_normalised_cpm", control = control, treatment = treatment, paired_data = paired_data)
      binned_delta_plots <- plot_binned_delta(binned_delta_data)
      
      png(filename = file.path(parent_dir, "plots/binned_plots/single_transcripts", dir, paste(treatment, gene, "binned delta.png")), width = 1000, height = 200)
      grid.arrange(binned_delta_plots[[1]], binned_delta_plots[[2]], binned_delta_plots[[3]],
                   nrow = 1, widths = c(1,2,1))
      dev.off()
    }
  }
  
  #single_nt
  if (plot_single_nt == T) {
    single_nt_list <- lapply(filtered_counts_list, splice_single_nt, region_lengths = region_lengths, region_cutoffs = region_cutoffs)
    
    summarised_single_nt_list <- lapply(single_nt_list, summarise_data, value = "single_nt_normalised_cpm", grouping = "window")
    
    do.call("rbind", summarised_single_nt_list) %>%
      group_by(grouping, condition, region) %>%
      summarise(average_counts = mean(mean_counts),
                sd_counts = sd(mean_counts)) %>%
      ungroup() -> summarised_single_nt
    
    #plot
    single_nt_line_plots <- plot_single_nt_lines(summarised_single_nt, SD=SD, control = control, treatment = treatment)
    
    png(filename = file.path(parent_dir, "plots/binned_plots/single_transcripts", dir, paste(treatment, gene, "single nt lines.png")), width = 1000, height = 200)
    grid.arrange(single_nt_line_plots[[1]], single_nt_line_plots[[2]], single_nt_line_plots[[3]], single_nt_line_plots[[4]], nrow = 1, widths = c(rep(1,3),1.2))
    dev.off()
    
    if (plot_delta == T) {
      #calculate and plot delta
      single_nt_delta_data <- calculate_single_nt_delta(single_nt_list, value = "single_nt_normalised_cpm", control = control, treatment = treatment, paired_data = paired_data)
      single_nt_delta_plots <- plot_single_nt_delta(single_nt_delta_data, SD = SD)
      
      png(filename = file.path(parent_dir, "plots/binned_plots/single_transcripts", dir, paste(treatment, gene, "single nt delta.png")), width = 1000, height = 200)
      grid.arrange(single_nt_delta_plots[[1]], single_nt_delta_plots[[2]], single_nt_delta_plots[[3]], single_nt_delta_plots[[4]], nrow = 1)
      dev.off()
    }
  }
}

plot_binned_heatmaps <- function(gene_names, remove_IDs, col_lims) {
  
  #get transcript IDs from gene list
  transcripts <- most_abundant_transcripts$transcript[most_abundant_transcripts$gene_sym %in% gene_names]
  
  #extract data from counts list
  filtered_counts_list <- lapply(counts_list, filter_transcripts, transcript_IDs = transcripts)
  
  #bin data
  single_transcript_binned_data_list <- lapply(filtered_counts_list, bin_data, region_lengths = region_lengths, region_cutoffs = c(50,300,0))
  single_transcript_binned_data <- do.call("rbind", single_transcript_binned_data_list)
  
  #plot delta
  single_transcript_binned_data %>%
    inner_join(most_abundant_transcripts, by = "transcript") %>%
    spread(key = condition, value = binned_counts) %>%
    mutate(delta = EFT226 - Ctrl) %>%
    group_by(gene_sym, bin, region) %>%
    summarise(mean_delta = mean(delta)) %>%
    filter(!(gene_sym %in% remove_IDs)) %>%
    filter(region == "CDS" | region == "UTR5") %>%
    mutate(plot_bin = case_when(region == "CDS" ~ bin,
                                region == "UTR5" ~ bin-25)) %>%
    ggplot(aes(x = plot_bin, y = gene_sym, fill = mean_delta))+
    geom_tile()+
    scale_fill_viridis(limits = col_lims)+
    geom_vline(xintercept = 0, lty=2)+
    theme_classic()+
    theme(axis.title.y = element_blank(),
          legend.title = element_blank())+
    xlab("Bin (relative to start codon)") -> binned_heatmap
  return(binned_heatmap)
}

plot_single_nt_heatmaps <- function(gene_names, remove_IDs, col_lims) {
  transcripts <- most_abundant_transcripts$transcript[most_abundant_transcripts$gene_sym %in% gene_names]
  UTR5_lens <- region_lengths[region_lengths$transcript %in% transcripts,]
  
  #extract data from counts list
  single_transcript_data_list <- lapply(counts_list, filter_transcripts, transcript_IDs = transcripts)
  single_transcript_data <- do.call("rbind", single_transcript_data_list)
  
  #plot delta
  single_transcript_data %>%
    inner_join(most_abundant_transcripts, by = "transcript") %>%
    inner_join(UTR5_lens, by = "transcript") %>%
    mutate(plot_position = Position - UTR5_len) %>%
    filter(plot_position <= 1000) %>%
    select(gene_sym, condition, replicate, plot_position, normalised_CPM) %>%
    spread(key = condition, value = normalised_CPM) %>%
    mutate(delta = EFT226 - Ctrl) %>%
    group_by(gene_sym, plot_position) %>%
    summarise(mean_delta = mean(delta)) -> summarised_single_transcript_data
  
  summarised_single_transcript_data %>%
    filter(!(gene_sym %in% remove_IDs)) %>%
    ggplot(aes(x = plot_position, y = gene_sym, fill = mean_delta))+
    geom_tile()+
    scale_fill_viridis(limits = col_lims)+
    geom_vline(xintercept = 0, lty=2)+
    theme_classic()+
    theme(axis.title.y = element_blank(),
          legend.title = element_blank())+
    xlab("nt (relative to start codon)") -> single_nt_heatmap
  return(single_nt_heatmap)
}
