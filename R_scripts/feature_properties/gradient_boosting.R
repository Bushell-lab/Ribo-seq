library(rsample)      # data splitting 
library(xgboost)      # a faster implementation of gbm
library(gbm)          # basic implementation
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # a java-based platform
library(pdp)          # model visualization
library(lime)         # model visualization
library(tidyverse)

#read in data----
merged_DEseq_data <- read_csv("merged_DESeq2.csv")

most_abundant_transcripts <- read_csv("most_abundant_transcripts_IDs.csv")

feature_properties <- read_csv(file = "gencode.v38.pc_transcripts_filtered_feature_properties.csv")

#merge data----
feature_properties %>%
  inner_join(most_abundant_transcripts, by = c("gene", "gene_sym", "transcript")) %>%
  inner_join(merged_DEseq_data[,c("gene", "gene_sym", "TE_log2FC")], by = c("gene", "gene_sym")) %>%
  column_to_rownames("gene") %>%
  select(-c("transcript", "gene_sym")) %>%
  mutate(Signal_seq = factor(Signal_seq),
         uORF = factor(uORF)) -> data_set

summary(data_set)

#Split data into training and test
data_split <- initial_split(data_set, prop = .7) #generates a random % (prop factor) split of your data set; that should include your to be predicted column and all other test-parameters)
data_train <- training(data_split) # will be used to optimise the model
data_test  <- testing(data_split) # will be used to test the model

# train GBM model on a manual basis to test performance
gbm.fit2 <- gbm(
  formula = TE_log2FC ~ ., # this is the column you want to have predicted, so log2FC or classes
  distribution = "gaussian", # this is basically the type of data of your to be predicted column
  data = data_train,
  n.trees = 1000, # the next 4 parameters are for optimising the model
  interaction.depth = 1, #connections between the tree nodes
  shrinkage = 0.01, #learning rate
  cv.folds = 5, #cross validation
  n.cores = NULL, # will use all cores by default
  verbose = FALSE
)  

# finds the minimum error in the series of n trees above
min_MSE <- which.min(gbm.fit2$cv.error)

# computes the RMSE of the model
sqrt(gbm.fit2$cv.error[min_MSE])

# plots error over the interation of the trees, good to check if model fits your data or noise
gbm.perf(gbm.fit2, method = "cv")

# plots the predictor strength of the parameters
par(mar = c(5, 12, 1, 1))
summary(
  gbm.fit2, 
  cBars = 11, #top 10 plotted
  method = relative.influence, # also can use permutation.test.gbm
  las = 2
)

# generates a df with the predictor values vs parameter
model <- tibble::as_tibble(gbm::summary.gbm(gbm.fit2, plotit = FALSE))
model

# optimse the model across all the parameters ####
# randomize data
random_index <- sample(1:nrow(data_train), nrow(data_train))
random_data_train <- data_train[random_index, ]

hyper_grid <- expand.grid(
  shrinkage = c(0.001, 0.01, 0.1),
  interaction.depth = c(1,5),
  n_trees =c(1000,5000),
  n.minobsinnode = c(5,10),
  bag.fraction = c(0.75), 
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)


# computes the model error for each model-parameter above in the matrix, care, this can take some time if too big matrix 
for(i in 1:nrow(hyper_grid)) {
  
  # reproducibility
  set.seed(123)
  
  # train model
  gbm.tune <- gbm(
    formula = TE_log2FC ~ .,
    distribution = "gaussian",
    data = random_data_train,
    n.trees =  hyper_grid$n_trees[i],
    cv.folds = 5,
    interaction.depth = hyper_grid$interaction.depth[i],
    shrinkage = hyper_grid$shrinkage[i],
    n.minobsinnode = hyper_grid$n.minobsinnode[i],
    n.cores = NULL, # will use all cores by default
    verbose = FALSE
  )
  
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- gbm.perf(gbm.tune, method = "cv") #saves the tree position for min error 
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$cv.error)) # saves the min error
  print(i)
}

hyper_grid %>% 
  dplyr::arrange(min_RMSE) %>%
  head(10)


# generate the final model using the optimised parameters from aboves matrix search ####

# for reproducibility
set.seed(123)

# train GBM model
gbm.fit.final <- gbm(
  formula = TE_log2FC ~ .,
  distribution = "gaussian",
  data = data_train,
  n.trees = 765,
  interaction.depth = 5,
  shrinkage = 0.01,
  cv.folds=5,
  n.minobsinnode = 10,
  bag.fraction = 0.75, 
  train.fraction = 0.7,
  n.cores = NULL, # will use all cores by default
  verbose = FALSE
) 

par(mar = c(5, 12, 1, 1))
summary(
  gbm.fit.final, 
  cBars = 10,
  method = relative.influence, # also can use permutation.test.gbm
  las = 2
)


# get error of the model
min_MSE <- which.min(gbm.fit.final$cv.error)
sqrt(gbm.fit.final$cv.error[min_MSE])


# plots error over the interation of the trees, good to check if model fits your data or noise
gbm.perf(gbm.fit2, method = "cv")

# get the list of predictor values 
model <- tibble::as_tibble(gbm::summary.gbm(gbm.fit.final, 
                                            plotit = FALSE))
# plot the results in a custom plot ####
plot <- ggplot(model, aes(fct_reorder(var,rel.inf), fill=fct_reorder(var,rel.inf)))+
  geom_bar(aes(weight=rel.inf), show.legend = F)+
  theme_classic()+
  theme(axis.text = element_text(size=18), 
        axis.title = element_text(size=20, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  scale_fill_viridis_d()+
  coord_flip(xlim=c(0,(nrow(model))))+
  labs(title="", x="feature", y="relative influence")

png(filename = file.path(parent_dir, "plots/feature_properties", "gradient_boosting.png"), width = 500, height = 300)
print(plot)
dev.off()

