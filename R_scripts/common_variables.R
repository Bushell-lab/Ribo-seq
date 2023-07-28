#Set the parent directory (this should be the same directory as is set in the common_variables.sh script
parent_dir <- 'Path/to/dir'

#set sample names
RPF_sample_names <- c('Ctrl_RPFs_1', 'Ctrl_RPFs_2', 'Ctrl_RPFs_3', 'Treatment_RPFs_1', 'Treatment_RPFs_2', 'Treatment_RPFs_3')
Total_sample_names <- c('Ctrl_Totals_1', 'Ctrl_Totals_2', 'Ctrl_Totals_3', 'Treatment_Totals_1', 'Treatment_Totals_2', 'Treatment_Totals_3')

RPF_sample_info <- data.frame(sample = RPF_sample_names,
                             condition = c(rep("Ctrl", 3), rep("Treatment", 3)),
                             replicate = factor(rep(c("1", "2", "3"), 2)))

Total_sample_info <- data.frame(sample = Total_sample_names,
                             condition = c(rep("Ctrl", 3), rep("Treatment", 3)),
                             replicate = factor(rep(c("1", "2", "3"), 2)))
