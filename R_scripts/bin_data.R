#load packages----
library(tidyverse)

#read in common variables----
source("/home/local/BICR/jwaldron/data/JWALDRON/Ribosome_profiling/Organoids/4AIn/Scripts/R_scripts/common_variables_VM.R")
#source("\\\\data.beatson.gla.ac.uk/data/JWALDRON/Ribosome_profiling/Organoids/4AIn/Scripts/R_scripts/common_variables.R")

#set the threshold for the average CDS counts a transcript has to have across all samples for it to be included
min_counts <- 50

#set the thresholds for region lengths
UTR5_min_len <- 50
CDS_min_len <- 300
UTR3_min_len <- 50

#read in functions----
source("/home/local/BICR/jwaldron/data/JWALDRON/Scripts/R/Functions/binning_RiboSeq_functions.R")
#source("\\\\data.beatson.gla.ac.uk/data/JWALDRON/Scripts/R/Functions/binning_RiboSeq_functions.R")

#read in data----
region_lengths <- read_csv(file = "/home/local/BICR/jwaldron/data/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_region_lengths.csv", col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))
#region_lengths <- read_csv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/FASTAs/mouse/GENCODE/vM27/transcript_info/gencode.vM27.pc_transcripts_region_lengths.csv", col_names = c("transcript", "UTR5_len", "CDS_len", "UTR3_len"))

#read in tpms----
#this reads in the tpms for each transcript and gathers them in tidy format
#It then inner joins the sample info data frame (from common variables) to obtain condition and replicate info
read_csv(file = file.path(parent_dir, "Analysis/DESeq2_output/tpms.csv")) %>%
  select(transcript, Total_sample_names) %>%
  gather(key = sample, value = tpm, all_of(Total_sample_names)) %>%
  inner_join(Total_sample_info, by = "sample") -> tpms

#sanity check----
#sum of all transcript tpms within each sample should equal 1,000,000
tpms %>%
  group_by(condition, replicate) %>%
  summarise(summed_tpm = sum(tpm))

#read in the CDS counts
#the following for loop reads in each final CDS counts file and renames the counts column by the sample name and saves each data frame to a list
data_list <- list()
for (sample in RPF_sample_names) {
  df <- read_csv(file = file.path(parent_dir, "Analysis/CDS_counts", paste0(sample, "_pc_final_CDS_counts_all_frames.csv")), col_names = T)
  colnames(df) <- c("transcript", sample)
  data_list[[sample]] <- df
}

#extract transcript IDs to keep----
#merge all data within the above list using reduce
#remove any transcripts with an NA in any sample (these will have had 0 counts)
#remove any transcripts with a mean count (across all samples) less than the min_counts threshold
#remove any transcripts with UTRs/CDSs shorter than thresholds defined above
#extract the transcript IDs
data_list %>%
  reduce(full_join, by = "transcript") %>%
  column_to_rownames("transcript") %>%
  drop_na() %>%
  filter(rowMeans(.) >= min_counts) %>%
  rownames_to_column("transcript") %>%
  inner_join(region_lengths, by = "transcript") %>%
  filter(UTR5_len >= UTR5_min_len & CDS_len >= CDS_min_len & UTR3_len >= UTR3_min_len) %>%
  pull(transcript) -> filtered_transcripts

#read in counts data----
#read in csvs using parLapply (parallel version of lapply)
counts_list <- list()
for (sample in RPF_sample_names) {
  
  #extract the condition and replicate from the sample name (this may need to be edited depending on how your sample names are structured)
  condition <- RPF_sample_info$condition[RPF_sample_info$sample == sample]
  replicate <- RPF_sample_info$replicate[RPF_sample_info$sample == sample]
  
  #get all the csv file names from the directory and filter
  read_csv(file = file.path(parent_dir, "Counts_files/csv_files", paste0(sample, "_pc_final_counts.csv"))) %>%
    filter(transcript %in% filtered_transcripts) %>%
    mutate(condition = rep(condition),
           replicate = rep(replicate)) %>%
    select(-sample) %>%
    normalise_data(tpms = tpms) -> counts_list[[sample]]
}

#return to parent directory
setwd(parent_dir)

#check counts list looks OK
summary(counts_list[[1]])
head(counts_list[[1]])

#sanity check----
#sum of all transcript CPMs within each sample should equal 1,000,000
sum_sample_counts <- function(df) {
  df %>%
    summarise(summed_counts = sum(CPM)) -> df2
  
  return(df2)
}

lapply(counts_list, sum_sample_counts)

#bin all data----
binned_list <- lapply(counts_list, bin_data, region_lengths = region_lengths, region_cutoffs = c(50,300,50), bins = c(25,50,25))
summary(binned_list[[1]])
head(binned_list[[1]])

#remove outliers----
#calculate the SD for each transcript across all replicates within each condition
do.call("rbind", binned_list) %>%
  group_by(transcript, condition, region, bin) %>%
  summarise(binned_normalised_cpm_sd = sd(binned_normalised_cpm),
            binned_cpm_sd = sd(binned_cpm)) -> transcript_SDs

summary(transcript_SDs)

n_distinct(transcript_SDs$transcript[is.na(transcript_SDs$binned_cpm_sd)])
n_distinct(transcript_SDs$transcript)

sd_quantiles <- data.frame(quantile = seq(0.95, 1, 0.00001),
                           binned_normalised_cpm_sd = quantile(transcript_SDs$binned_normalised_cpm_sd, probs = seq(0.95, 1, 0.00001), na.rm = T),
                           binned_cpm_sd = quantile(transcript_SDs$binned_cpm_sd, probs = seq(0.95, 1, 0.00001), na.rm = T))

#plot all SD quantiles
png(filename = file.path(parent_dir, "plots/binned_plots/normalisation/binned_normalised_cpm_sd_quantiles_pre_filtering.png"))
sd_quantiles %>%
  ggplot(aes(quantile, binned_normalised_cpm_sd))+
  geom_point()+
  theme_classic()
dev.off()

png(filename = file.path(parent_dir, "plots/binned_plots/normalisation/binned_cpm_sd_quantiles_pre_filtering.png"))
sd_quantiles %>%
  ggplot(aes(quantile, binned_cpm_sd))+
  geom_point()+
  theme_classic()
dev.off()

#extract transcripts with a binned SD above a quantile threshold
#set outlier quantile (0.99 will remove the transcripts with the top 1% SD) and make sure this is set appropriatly based on the number of transcripts being removed
perc <- 1 - 0.00003

normalised_cpm_discard_transcripts <- transcript_SDs$transcript[transcript_SDs$binned_normalised_cpm_sd >= quantile(transcript_SDs$binned_normalised_cpm_sd, probs = perc, na.rm = T)]
cpm_discard_transcripts <- transcript_SDs$transcript[transcript_SDs$binned_cpm_sd >= quantile(transcript_SDs$binned_cpm_sd, probs = perc, na.rm = T)]
discard_transcripts <- unique(c(normalised_cpm_discard_transcripts, cpm_discard_transcripts))
paste(length(discard_transcripts), "outliers to be removed")

#remove outliers from counts list and re-normalise
for (i in 1:length(counts_list)) {
  counts_list[[i]] <- counts_list[[i]][!(counts_list[[i]]$transcript %in% discard_transcripts),]
  counts_list[[i]] <- normalise_data(counts_list[[i]], tpms = tpms)
}

summary(counts_list[[1]])
head(counts_list[[1]])

#bin all data again having removed the outliers
binned_list <- lapply(counts_list, bin_data, region_lengths = region_lengths, region_cutoffs = c(50,300,50), bins = c(25,50,25))
summary(binned_list[[1]])
head(binned_list[[1]])

#plot quantiles of SD again to check all outliers have been removed
do.call("rbind", binned_list) %>%
  group_by(transcript, condition, region, bin) %>%
  summarise(binned_normalised_cpm_sd = sd(binned_normalised_cpm),
            binned_cpm_sd = sd(binned_cpm)) -> transcript_SDs

sd_quantiles <- data.frame(quantile = seq(0.95, 1, 0.00001),
                           binned_normalised_cpm_sd = quantile(transcript_SDs$binned_normalised_cpm_sd, probs = seq(0.95, 1, 0.00001), na.rm = T),
                           binned_cpm_sd = quantile(transcript_SDs$binned_cpm_sd, probs = seq(0.95, 1, 0.00001), na.rm = T))

#plot all SD quantiles
png(filename = file.path(parent_dir, "plots/binned_plots/normalisation/binned_normalised_cpm_sd_quantiles_post_filtering.png"))
sd_quantiles %>%
  ggplot(aes(quantile, binned_normalised_cpm_sd))+
  geom_point()+
  theme_classic()
dev.off()

png(filename = file.path(parent_dir, "plots/binned_plots/normalisation/binned_cpm_sd_quantiles_post_filtering.png"))
sd_quantiles %>%
  ggplot(aes(quantile, binned_cpm_sd))+
  geom_point()+
  theme_classic()
dev.off()

#single nt----
single_nt_list <- lapply(counts_list, splice_single_nt, region_lengths = region_lengths)
summary(single_nt_list[[1]])
head(single_nt_list[[1]])

#plot quantiles of SD again to check all outliers have been removed
do.call("rbind", single_nt_list) %>%
  group_by(transcript, condition, region, window) %>%
  summarise(single_nt_normalised_cpm_sd = sd(single_nt_normalised_cpm),
            single_nt_cpm_sd = sd(single_nt_cpm)) -> transcript_SDs

sd_quantiles <- data.frame(quantile = seq(0.95, 1, 0.00001),
                           single_nt_normalised_cpm_sd = quantile(transcript_SDs$single_nt_normalised_cpm_sd, probs = seq(0.95, 1, 0.00001), na.rm = T),
                           single_nt_cpm_sd = quantile(transcript_SDs$single_nt_cpm_sd, probs = seq(0.95, 1, 0.00001), na.rm = T))

#plot all SD quantiles
png(filename = file.path(parent_dir, "plots/binned_plots/normalisation/single_nt_normalised_cpm_sd_quantiles_post_filtering.png"))
sd_quantiles %>%
  ggplot(aes(quantile, single_nt_normalised_cpm_sd))+
  geom_point()+
  theme_classic()
dev.off()

png(filename = file.path(parent_dir, "plots/binned_plots/normalisation/single_nt_cpm_sd_quantiles_post_filtering.png"))
sd_quantiles %>%
  ggplot(aes(quantile, single_nt_cpm_sd))+
  geom_point()+
  theme_classic()
dev.off()

#save lists----
save(file = file.path(parent_dir, "Counts_files/R_objects/binned_list.Rdata"), binned_list)
save(file = file.path(parent_dir, "Counts_files/R_objects/single_nt_list.Rdata"), single_nt_list)
