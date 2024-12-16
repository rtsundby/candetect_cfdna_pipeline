#!/usr/bin/env Rscript

#Comparison individual F-profiles ability to distinguish two cohorts (one versus one or one versus all). 10 repeats of 5-fold CV

#Parse arguments
args = commandArgs(trailingOnly = TRUE)
fprofiles_tsv = args[1]
libraries_tsv = args[2]
plot_out = args[3]
csv_out = args[4]

# Load necessary packages
library(readr)
library(dplyr)
library(pROC)
library(ggplot2)
library(knitr)
library(tidyr)
library(cowplot)
library(caret)

# Load data frames, join data frames
fprofiles <- read_tsv(fprofiles_tsv)
libraries <- read_tsv(libraries_tsv)

colnames(fprofiles)[1]<-"fprofile"
colnames(fprofiles)[-1] <- sub("_motifs$", "", colnames(fprofiles)[-1]) #Note the fprofiles_tsv appends _motifs to library names. Pre-processing to remove
fprofiles_long <- pivot_longer(fprofiles, cols = -1, names_to = "library", values_to = "pct")
data <- fprofiles_long %>%
  left_join(libraries, by = "library") %>% na.omit(.)

#Define fprofile names and cohorts
fprofiles_list <- unique(fprofiles$fprofile)

cohorts <- unique(data$cohort)
if (length(cohorts) != 2) {
  stop("This script is designed for the comparison of two cohorts")
}

# Initialize an empty data frame for the results table and a ggplot object for overlapping ROC
roc_results <- data.frame(fprofile = character(), AUC = numeric(), J_cutpoint = numeric(), stringsAsFactors = FALSE)

roc_plot <- ggplot()

# Loop through each fprofile to calculate and plot ROC
for (fprofile_name in fprofiles_list) {
  
  subset_data <- data %>%
    filter(fprofile == !!fprofile_name) %>%  
    select(pct, cohort)
  
  # Set up cross-validation: 10 repeats of 5-fold CV
  train_control <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 10,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  # Perform CV and collect AUC values
  auc_values <- numeric()
  
  tryCatch({
    for (i in 1:10) {
      folds <- createFolds(subset_data$cohort, k = 5, returnTrain = TRUE)
      for (fold in folds) {
        train <- subset_data[fold, ]
        test <- subset_data[-fold, ]
        
        roc_obj <- roc(test$cohort, test$pct, levels = cohorts)
        auc_values <- c(auc_values, auc(roc_obj))
      }
    }
    
    mean_auc <- mean(auc_values)
    sd_auc <- sd(auc_values)
    
    roc_results <- rbind(roc_results, data.frame(fprofile = fprofile_name, AUC = mean_auc, AUC_sd = sd_auc))
    
    # Prepare for plotting
    roc_plot <- roc_plot + 
      geom_smooth(aes(x = 1 - roc_obj$specificities, y = roc_obj$sensitivities, color = fprofile_name), 
                  method = "loess", span = 0.2, size = 1.5, se = FALSE, linetype = "solid")
  }, error = function(e) {
    message(paste("Error in ROC computation for fprofile", fprofile_name, ":", e$message))
  })
}


# Add diagonal line and labels
roc_plot <- roc_plot + 
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  theme_half_open() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(name = "F-profile", values = rainbow(length(fprofiles_list)), 
                     labels = sapply(1:length(fprofiles_list), function(i) {
                       fprofile <- fprofiles_list[i]
                       auc_val <- round(roc_results$AUC[roc_results$fprofile == fprofile], 2)
                       paste(fprofile, "AUC:", auc_val)
                     }))

# Print results table
kable(roc_results, caption = "AUC and J-cutpoint for Each Fprofile")

print(roc_plot)

write.csv(roc_results, csv_out, row.names = FALSE)
ggsave(plot_out, plot = roc_plot, width = 10, height = 8)
