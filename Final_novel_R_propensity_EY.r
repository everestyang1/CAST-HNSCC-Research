rm(list=ls())

library(caret)
library(car)
library(boot)
library(pscl)
library(tidyverse)
library(Hmisc)
library(psych)
library(survival)
library(survminer)
library(pdp)
library(glmnet)
library(grf)
library(devEMF)
library(matrixStats)
library(vioplot)
library(fastshap)
library(MASS)
library(corrplot)
library(GGally)
library(dplyr)

plot_effects_distribution <- function(results, method, horizon) {
  effects_df <- data.frame(effects = results$predictions)
  
  if (var(effects_df$effects) == 0) {
    effects_df$effects <- jitter(effects_df$effects, amount = 0.01)
  }
  
  # Calculate the number of bins using the Freedman-Diaconis rule
  iqr <- IQR(effects_df$effects)
  n <- length(effects_df$effects)
  bin_width <- 2 * iqr / (n^(1/3))
  num_bins <- ceiling((max(effects_df$effects) - min(effects_df$effects)) / bin_width)
  
  print(paste("ATE:", results$ate[1]))
  print(paste("Mean of predictions:", mean(effects_df$effects)))
  
  p <- ggplot(effects_df, aes(x = effects)) +
    geom_histogram(bins = num_bins, fill = "lightblue", color = "black") +

# here ATE does not always correspond to mean of sample predictions because ATE is not just the 
#  mean of sample estimates, it is calculated by doubly robust methods
#  To avoid confusion I deleted the red line
#    geom_vline(xintercept = results$ate[1], color = "red", linetype = "dashed") +
    theme_bw() +
    labs(title = sprintf("Distribution of Treatment Effects at %d month Horizon (%s)", 
                         horizon, method),
         x = "Individual treatment effect",
         y = "Count")
  
  ggsave(sprintf("treatment_effects_%s_%dm.png", method, horizon), p, 
         width = 7, height = 7)
}


plot_treatment_effects_by_covariate <- function(predictions, X, var_name, 
                                              keep_indices = NULL,
                                              title = "Treatment Effects by Covariate") {
    if (!is.null(keep_indices)) {
        X <- X[keep_indices, ]
    }
    
    effects_df <- data.frame(
        effect = predictions$predictions,
        covariate = X[, var_name]
    )
    
    p <- ggplot(effects_df, aes(x = covariate, y = effect)) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "loess", se = TRUE) +
        theme_bw() +
        labs(title = paste(title, "-", var_name),
             x = var_name,
             y = "Estimated Treatment Effect")
    
    ggsave(paste0("effects_by_", tolower(var_name), ".png"), p, width = 8, height = 6)
    return(p)
}

# data simulation function with simple treatment effect
simulate_simple_data <- function(n_samples, true_effect, seed = 42) {
  set.seed(seed)
  
  data <- data.frame(
    Age = rnorm(n_samples, mean = 60, sd = 10),
    Sex = rbinom(n_samples, 1, 0.5),
    Stage = sample(1:4, n_samples, replace = TRUE),
    Chemotherapy = rbinom(n_samples, 1, 0.5),
    HPV = rbinom(n_samples, 1, 0.6),
    Smoking = rbinom(n_samples, 1, 0.25),
    TumorSite_OralCavity = as.numeric(sample(1:4, n_samples, replace = TRUE) == 1),
    TumorSite_Oropharynx = as.numeric(sample(1:4, n_samples, replace = TRUE) == 2),
    TumorSite_Larynx = as.numeric(sample(1:4, n_samples, replace = TRUE) == 3),
    TumorSite_Hypopharynx = as.numeric(sample(1:4, n_samples, replace = TRUE) == 4),
    Comorbidities_0 = as.numeric(sample(0:2, n_samples, replace = TRUE, prob = c(0.5, 0.3, 0.2)) == 0),
    Comorbidities_1 = as.numeric(sample(0:2, n_samples, replace = TRUE, prob = c(0.5, 0.3, 0.2)) == 1),
    Comorbidities_2plus = as.numeric(sample(0:2, n_samples, replace = TRUE, prob = c(0.5, 0.3, 0.2)) == 2),
    NeckDissection = rbinom(n_samples, 1, 0.4),
    FeedingTube = rbinom(n_samples, 1, 0.2),
    Tracheostomy = rbinom(n_samples, 1, 0.1),
    AlcoholUse_None = as.numeric(sample(1:3, n_samples, replace = TRUE) == 1),
    AlcoholUse_Moderate = as.numeric(sample(1:3, n_samples, replace = TRUE) == 2),
    AlcoholUse_Heavy = as.numeric(sample(1:3, n_samples, replace = TRUE) == 3),
    Education_HighSchool = as.numeric(sample(1:3, n_samples, replace = TRUE) == 1),
    Education_College = as.numeric(sample(1:3, n_samples, replace = TRUE) == 2),
    Education_Graduate = as.numeric(sample(1:3, n_samples, replace = TRUE) == 3),
    MaritalStatus_Single = as.numeric(sample(1:3, n_samples, replace = TRUE) == 1),
    MaritalStatus_Married = as.numeric(sample(1:3, n_samples, replace = TRUE) == 2),
    MaritalStatus_Divorced = as.numeric(sample(1:3, n_samples, replace = TRUE) == 3),
    BMI = rnorm(n_samples, mean = 25, sd = 4),
    PerformanceStatus_0 = as.numeric(sample(0:4, n_samples, replace = TRUE) == 0),
    PerformanceStatus_1 = as.numeric(sample(0:4, n_samples, replace = TRUE) == 1),
    PerformanceStatus_2 = as.numeric(sample(0:4, n_samples, replace = TRUE) == 2),
    PerformanceStatus_3 = as.numeric(sample(0:4,n_samples ,replace=TRUE)==3), 
    PerformanceStatus_4 = as.numeric(sample(0:4,n_samples ,replace=TRUE)==4)
  ) 
  
  # Even better balanced treatment assignment
  prob_Cause <- plogis(-7 -
      (0.3 * scale(data$Age)) + 
      (3 * data$Stage) - 
      (2 * data$HPV)
  )
  
  data$Cause <- rbinom(n_samples, 1, prob_Cause)
  
  # Survival time with true effect and some other effects
  Survival_time <- round(pmax(rnorm(n_samples, mean=15+true_effect*data$Cause
            +5*data$Chemotherapy -2*data$Stage + 5*data$HPV, sd=8),1),0)
  
  # Survival status with censoring
  censoring_prob <- function(time) {
    pmin(0.01 + 0.5*(time / max(time)), 0.99)
  }
  
  Survival_status <- rbinom(n_samples, 1, 1-censoring_prob(Survival_time))
  
  data$time <- Survival_time
  data$event <- Survival_status
  
  # Print detailed diagnostics
  event_rate <- mean(data$event)
  treatment_rate <- mean(data$Cause)
  cat(sprintf("Event rate: %.1f%%\n", event_rate * 100))
  cat(sprintf("Treatment rate: %.1f%%\n", treatment_rate * 100))
  
  # Calculate true survival probabilities at multiple timepoints
  for(t in c(12*1, 12*2, 12*3, 12*4)) {
    true_surv_control <- mean(data$time[data$Cause == 0] > t)
    true_surv_treated <- mean(data$time[data$Cause == 1] > t)
    cat(sprintf("True %d-month survival - Control: %.1f%% , Treated: %.1f%% (Diff: %.1f%%)\n", 
                t ,true_surv_control *100 ,true_surv_treated *100 ,
                (true_surv_treated - true_surv_control)*100))
  }
  
  # Calculate and print median survival times for treated and control groups
  median_control <- median(data$time[data$Cause == 0])
  median_treated <- median(data$time[data$Cause == 1])
  cat(sprintf("Median survival time - Control: %.2f months , Treated: %.2f months\n", 
              median_control ,median_treated))
  
  # Check balance in covariates
  cat("\nCovariate balance:\n")
  for(var in c("Age" ,"Stage" ,"HPV", "Sex")) {
    mean_control <- mean(data[[var]][data$Cause == 0])
    mean_treated <- mean(data[[var]][data$Cause == 1])
    cat(sprintf("%s - Control: %.2f , Treated: %.2f\n", var ,mean_control ,mean_treated))
  }
  
  return(data)
}


# call the simulation function

simulated_data <- simulate_simple_data(n_samples = 3000, true_effect = 10)
str(simulated_data)
write.csv(simulated_data, file = file.path("simulated_data.csv"), row.names = FALSE)

# some correlations
summary(simulated_data$Cause)
cor(simulated_data$Cause, simulated_data$Stage)
cor(simulated_data$Cause, simulated_data$Age)
cor(simulated_data$Cause, simulated_data$HPV)

########

# function to split data into train and test parts, separate causes, outcomes and covariates
prepare_cancer_data <- function(data_set, seed = 99) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Create training and test sets
  index <- caret::createDataPartition(data_set$event, p=0.75, list=FALSE)
  trainSet <- data_set[index,]
  testSet <- data_set[-index,]
  
  # Print dataset sizes and survival status distribution
  cat("Training set size:", nrow(trainSet), "\n")
  cat("Test set size:", nrow(testSet), "\n")
  cat("Training set survival status:\n")
  print(table(trainSet$event))
  cat("Test set survival status:\n")
  print(table(testSet$event))
  
  # Create the necessary matrices for the causal forest
  trainSet_X <- as.data.frame(subset(trainSet, select=-c(time, event, Cause)))
  trainSet_W <- trainSet$Cause
  trainSet_times <- trainSet$time
  trainSet_events <- trainSet$event
  
  testSet_X <- as.data.frame(subset(testSet, select=-c(time, event, Cause)))
  testSet_W <- testSet$Cause
  testSet_times <- testSet$time
  testSet_events <- testSet$event
  
  # Fit logistic regression with elastic net penalty to get propensity scores for Cause
  X_matrix_train <- model.matrix(~ . - Cause - time - event, data = trainSet)[,-1]
  Y_vector_train <- trainSet$Cause
  X_matrix_test <- model.matrix(~ . - Cause - time - event, data = testSet)[,-1]
  
  return(list(
    trainSet = trainSet,
    testSet = testSet,
    trainSet_X = trainSet_X,
    trainSet_W = trainSet_W,
    trainSet_times = trainSet_times,
    trainSet_events = trainSet_events,
    testSet_X = testSet_X,
    testSet_W = testSet_W,
    testSet_times = testSet_times,
    testSet_events = testSet_events,
    X_matrix_train = X_matrix_train,
    Y_vector_train = Y_vector_train,
    X_matrix_test = X_matrix_test
  ))
}

# Call the function on simulated data
#prepared_data <- prepare_cancer_data(simulated_data)

# or read real data instead
real_data_raw <- read.csv("selected_data_with_sites_processed.csv")
str(real_data_raw)

# calculate BEDs on this data set
#   calculate BED_DI and BED_DD variants


#   parameters

# cell killing linear term (1/Gy)
alpha <- 0.2

# cell killing alpha/beta ratio (Gy)
r <- 10

# accelerated repopulation rate (exponential, 1/day) for DI model
L_DI <- 0.2

# accelerated maximum repopulation rate (exponential, 1/day) for DD model
L_DD <- 0.5

# time after treatment begins when accelerated repopulation starts, DI model (days)
Tk <- 28

# number of dose fractions
m <- real_data_raw$Fx

# dose/fraction (Gy)
d <- real_data_raw$d_Frac

# total treatment time (days)
T <- real_data_raw$Total_days_RT

# dose-independent model BED
real_data_raw$BED_DI <- m*d*(1 + d/r) - L_DI*log(1 + exp(3*T - 3*Tk))/(3*alpha)

# dose-dependent model BED
real_data_raw$BED_DD <- (r*L_DD*(exp(-m*alpha*d*(d + r)/(T*r)) - 1)*log(1 + exp(3*(47*m*d^2 + 47*d*m*r 
          + 28*r*log((1 + exp(141 - 3*Tk))^(-L_DI/((exp(-(84*alpha)/47) - 1)*L_DD)) - 1) 
          - 3948*r)*T/(47*m*d*(d + r)))) + 3*m*alpha*d*(d + r))/(3*r*alpha)

# rename data
real_data <- real_data_raw

# remove redundant features which were combined into BEDs
real_data <- real_data[, !(names(real_data) %in% c("Fx", "d_Frac", "Total_days_RT"))]

summary(real_data)


write.csv(real_data, file = "real_data.csv", row.names = F)


# use real data
prepared_data <- prepare_cancer_data(real_data)

trainSet <- prepared_data$trainSet
testSet <- prepared_data$testSet
trainSet_X <- prepared_data$trainSet_X
trainSet_W <- prepared_data$trainSet_W
trainSet_times <- prepared_data$trainSet_times
trainSet_events <- prepared_data$trainSet_events
testSet_X <- prepared_data$testSet_X
testSet_W <- prepared_data$testSet_W
testSet_times <- prepared_data$testSet_times
testSet_events <- prepared_data$testSet_events
X_matrix_train <- prepared_data$X_matrix_train
Y_vector_train <- prepared_data$Y_vector_train
X_matrix_test <- prepared_data$X_matrix_test


##########
# visualize KM plots
# Prepare the data
data_KM <- data.frame(trainSet_X, trainSet_times, trainSet_W, trainSet_events)
# Create a Surv object
surv_object_KM <- Surv(time = trainSet_times, event = trainSet_events)
# Fit the Kaplan-Meier survival curves
fit_KM <- survfit(surv_object_KM ~ trainSet_W, data = data_KM)
# Create the Kaplan-Meier plot with ggsurvplot
p_KM <- ggsurvplot(fit_KM, data = data_KM,
                   pval = TRUE, conf.int = TRUE,
                   risk.table = TRUE,
                   ggtheme = theme_bw(),
                   xlab = "Time",
                   ylab = "Survival Probability",
                   title = "Kaplan-Meier Survival Curves")
# Customize the plot to make fonts bigger and bold
p_KM$plot <- p_KM$plot + theme(
  plot.title = element_text(size = 20, face = "bold"),
  axis.title.x = element_text(size = 16, face = "bold"),
  axis.title.y = element_text(size = 16, face = "bold"),
  axis.text = element_text(size = 14, face = "bold"),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12, face = "bold"),
  strip.text = element_text(size = 14, face = "bold")
)
p_KM$table <- p_KM$table + theme(
  
  axis.title.x = element_text(size = 16, face = "bold"),
  axis.title.y = element_text(size = 16, face = "bold"),
  axis.text = element_text(size = 14, face = "bold")
)
# Save the plot
ggsave(filename = "Kaplan_Meier_Survival_Curves.png", plot = p_KM$plot, dpi = 300, width = 10,
       height = 8)
# Save the risk table
ggsave(filename = "Kaplan_Meier_Survival_Curves_risk_table.png", plot = p_KM$table, dpi = 300,
       width = 10, height = 4)

#######
# propensity score model

fit_propensity_model <- function(trainSet, testSet, seed = 1234, alpha_range = seq(0.01, 0.99, by = 0.01), output_dir = ".") {
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE)
  
  # Fit logistic regression with elastic net penalty to get propensity scores for Cause
  X_matrix_train <- model.matrix(~ . - Cause - time - event, data = trainSet)[,-1]
  Y_vector_train <- trainSet$Cause
  X_matrix_test <- model.matrix(~ . - Cause - time - event, data = testSet)[,-1]
  
  set.seed(seed)
  cv_models <- list()
  
  for (alpha_val in alpha_range) {
    cv_model <- cv.glmnet(X_matrix_train, Y_vector_train,
                          family = "binomial",
                          type.measure = "class",
                          alpha = alpha_val,
                          nfolds = 10)
    cv_models[[as.character(alpha_val)]] <- cv_model
  }
  
  optimal_lambdas <- sapply(cv_models, function(model) model$lambda.min)
  best_alpha <- alpha_range[which.min(optimal_lambdas)]
  optimal_lambda <- cv_models[[as.character(best_alpha)]]$lambda.min
  final_model <- glmnet(X_matrix_train, Y_vector_train,
                        family = "binomial",
                        alpha = best_alpha,
                        lambda = optimal_lambda)
  cv_glmnet_model <- final_model
  
  # model coefficients
  cv_glmnet_model_coefs <- as.data.frame(as.matrix(coef(cv_glmnet_model)))
  write.csv(cv_glmnet_model_coefs, file = file.path(output_dir, "cv_glmnet_model_coefs.csv"), row.names = TRUE)
  
  # Predict propensity scores
  W.hat_train <- predict(cv_glmnet_model, newx = X_matrix_train, type = "response", s = optimal_lambda)
  W.hat_test <- predict(cv_glmnet_model, newx = X_matrix_test, type = "response", s = optimal_lambda)
  
  # propensity scores to datasets
  trainSet_with_scores <- as.data.frame(cbind(trainSet, W.hat_train))
  colnames(trainSet_with_scores) <- c(colnames(trainSet), "W_hat")
  
  testSet_with_scores <- as.data.frame(cbind(testSet, W.hat_test))
  colnames(testSet_with_scores) <- c(colnames(testSet), "W_hat")
  
  # datasets with scores
  write.csv(trainSet_with_scores, file.path(output_dir, "trainSet_with_scores.csv"), row.names = FALSE)
  write.csv(testSet_with_scores, file.path(output_dir, "testSet_with_scores.csv"), row.names = FALSE)
  
  # correlation matrix
  varscor_trainSet_with_scores_pearson <- corr.test(trainSet_with_scores, method = "pearson", adjust = "bonf", alpha = .05, ci = FALSE)
  varscor_trainSet_with_scores_pearson_p <- varscor_trainSet_with_scores_pearson$p
  
  write.csv(varscor_trainSet_with_scores_pearson_p, file = file.path(output_dir, "varscor_trainSet_with_scores_pearson_p.csv"), row.names = TRUE)
  write.csv(varscor_trainSet_with_scores_pearson$r, file = file.path(output_dir, "varscor_trainSet_with_scores_pearson_r.csv"), row.names = TRUE)
  
#  emf(file.path(output_dir, "correlation_matrix_pearson_propensity_scores_train.emf"), width = 7, height = 7, bg = "transparent",
#      fg = "black", pointsize = 8, family = "Arial", coordDPI = 600)
    png(file.path(output_dir, "correlation_matrix_pearson_propensity_scores_train.png"), 
        width = 7, height = 7, units = "in", res = 600, bg = "white", pointsize = 8)
  
    corrplot(varscor_trainSet_with_scores_pearson$r, p.mat = varscor_trainSet_with_scores_pearson$p, method = 'circle', 
           tl.col = "black", type = "upper", sig.level = 0.05, pch.cex = 0.6, cl.cex = 1, tl.cex = 1, insig = 'pch', 
           pch = 19, pch.col = "white", diag = FALSE, font = 1)

    dev.off()
  
    png(file.path(output_dir, "propensity_scores_histogram_train.png"), 
        width = 7, height = 7, units = "in", res = 600, bg = "white", pointsize = 8)
    
    hist(trainSet_with_scores$W_hat, breaks=30, xlab="Popensity score")
    
    dev.off()
    
  propensity_score_cutoff <- 0.1    
    
  trainSet_with_scores_filtered <- subset(trainSet_with_scores, W_hat >= propensity_score_cutoff & W_hat <= (1-propensity_score_cutoff))
  testSet_with_scores_filtered <- subset(testSet_with_scores, W_hat >= propensity_score_cutoff & W_hat <= (1-propensity_score_cutoff))
  
  trainSet_X <- as.data.frame(subset(trainSet_with_scores_filtered, select = -c(W_hat, time, event, Cause)))
  trainSet_W <- trainSet_with_scores_filtered$Cause
  trainSet_times <- trainSet_with_scores_filtered$time
  trainSet_events <- trainSet_with_scores_filtered$event
  W_hat_train_adj2 <- trainSet_with_scores_filtered$W_hat
  
  testSet_X <- as.data.frame(subset(testSet_with_scores_filtered, select = -c(W_hat, time, event, Cause)))
  testSet_W <- testSet_with_scores_filtered$Cause
  testSet_times <- testSet_with_scores_filtered$time
  testSet_events <- testSet_with_scores_filtered$event
  W_hat_test_adj2 <- testSet_with_scores_filtered$W_hat
  
  write.csv(trainSet_with_scores_filtered, file.path(output_dir, "trainSet_with_scores_filtered.csv"), row.names = FALSE)
  write.csv(testSet_with_scores_filtered, file.path(output_dir, "testSet_with_scores_filtered.csv"), row.names = FALSE)
  
  return(list(
    best_alpha = best_alpha,
    optimal_lambda = optimal_lambda,
    cv_glmnet_model = cv_glmnet_model,
    trainSet_with_scores = trainSet_with_scores,
    testSet_with_scores = testSet_with_scores,
    trainSet_with_scores_filtered = trainSet_with_scores_filtered,
    testSet_with_scores_filtered = testSet_with_scores_filtered,
    trainSet_X = trainSet_X,
    trainSet_W = trainSet_W,
    trainSet_times = trainSet_times,
    trainSet_events = trainSet_events,
    W_hat_train_adj2 = W_hat_train_adj2,
    testSet_X = testSet_X,
    testSet_W = testSet_W,
    testSet_times = testSet_times,
    testSet_events = testSet_events,
    W_hat_test_adj2 = W_hat_test_adj2
  ))
}

propensity_results <- fit_propensity_model(trainSet, testSet, output_dir = "propensity_score_output")

best_alpha <- propensity_results$best_alpha
optimal_lambda <- propensity_results$optimal_lambda
cv_glmnet_model <- propensity_results$cv_glmnet_model
trainSet_with_scores <- propensity_results$trainSet_with_scores
testSet_with_scores <- propensity_results$testSet_with_scores
trainSet_with_scores_filtered <- propensity_results$trainSet_with_scores_filtered
testSet_with_scores_filtered <- propensity_results$testSet_with_scores_filtered
trainSet_X <- propensity_results$trainSet_X
trainSet_W <- propensity_results$trainSet_W
trainSet_times <- propensity_results$trainSet_times
trainSet_events <- propensity_results$trainSet_events
W_hat_train_adj2 <- propensity_results$W_hat_train_adj2
testSet_X <- propensity_results$testSet_X
testSet_W <- propensity_results$testSet_W
testSet_times <- propensity_results$testSet_times
testSet_events <- propensity_results$testSet_events
W_hat_test_adj2 <- propensity_results$W_hat_test_adj2

print(paste("Best alpha:", best_alpha))
print(paste("Optimal lambda:", optimal_lambda))
print(paste("Original training set size:", nrow(trainSet_with_scores)))
print(paste("Filtered training set size:", nrow(trainSet_with_scores_filtered)))
print(paste("Original test set size:", nrow(testSet_with_scores)))
print(paste("Filtered test set size:", nrow(testSet_with_scores_filtered)))


########
# causal survival forests implementation


m_trainSet_X <- as.matrix(trainSet_X)
m_trainSet_W <- as.matrix(trainSet_W)
m_trainSet_times <- as.matrix(trainSet_times) 
m_trainSet_events <-as.matrix(trainSet_events)
m_testSet_X <-as.matrix(testSet_X)
m_testSet_W <- as.matrix(testSet_W)
m_testSet_times <-as.matrix(testSet_times)
m_testSet_events <- as.matrix(testSet_events)
m_W_hat_train_adj2 <- as.matrix(W_hat_train_adj2)
m_W_hat_test_adj2 <- as.matrix(W_hat_test_adj2)

implement_causal_forests <- function(m_trainSet_X, m_trainSet_W, m_trainSet_times, m_trainSet_events,
                                     m_testSet_X, m_testSet_W, m_testSet_times, m_testSet_events,
                                     m_W_hat_train_adj2, m_W_hat_test_adj2,
                                     n_trees_val = 5000, 
                                     horizons = seq(12, 120, by=12),
                                     output_dir = ".") {
  
  # Input validation
  if (!is.matrix(m_trainSet_X) || !is.matrix(m_testSet_X)) 
    stop("m_trainSet_X and m_testSet_X must be matrices")
  
  if (any(W_hat_train_adj2 < 0) || any(W_hat_train_adj2 > 1))
    stop("W_hat probabilities must be between 0 and 1")
  
  if (any(horizons <= 0))
    stop("Horizons must be positive")
  
  # Create output directory safely
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Pre-allocate results data frames
  results_SP <- data.frame(
    horizon_sel = integer(length(horizons)),
    ATE_estimate_train_SP = numeric(length(horizons)),
    ATE_se_train_SP = numeric(length(horizons)),
    ATE_estimate_test_SP = numeric(length(horizons)),
    ATE_se_test_SP = numeric(length(horizons))
  )
  
  results_RMST <- data.frame(
    horizon_sel = integer(length(horizons)),
    ATE_estimate_train_RMST = numeric(length(horizons)),
    ATE_se_train_RMST = numeric(length(horizons)),
    ATE_estimate_test_RMST = numeric(length(horizons)),
    ATE_se_test_RMST = numeric(length(horizons))
  )
  
  # Set seed for reproducibility
  set.seed(123)
  
  for (i in seq_along(horizons)) {
    horizon <- horizons[i]
    tryCatch({
      # Survival Probability model
      csf_model_SP <- causal_survival_forest(
        X = m_trainSet_X, 
        Y = m_trainSet_times, 
        W = m_trainSet_W,
        D = m_trainSet_events, 
        W.hat = as.vector(W_hat_train_adj2), 
        num.trees = n_trees_val, 
        target = "survival.probability", 
        horizon = horizon, 
#        honesty = TRUE,
#        min.node.size = 5,
#        alpha = 0.05,
#        imbalance.penalty = 0.1,
#        stabilize.splits = TRUE,
        tune.parameters = "all", 
        seed = 123
      )
      
      ate_train_SP <- average_treatment_effect(csf_model_SP)
      csf_pred_test_SP <- predict(csf_model_SP, m_testSet_X, 
                                  W.hat = as.vector(W_hat_test_adj2), 
                                  estimate.variance = TRUE)

      plot_effects_distribution(
        list(
            predictions = csf_pred_test_SP$predictions,
            ate = ate_train_SP
        ),
        "SP",
        horizon
    )
      
      # Store results
      results_SP[i, ] <- c(
        horizon,
        ate_train_SP[1],
        ate_train_SP[2],
        mean(csf_pred_test_SP$predictions),
        mean(sqrt(csf_pred_test_SP$variance.estimates))
      )
      
      # Similar block for RMST model
      csf_model_RMST <- causal_survival_forest(
        X = m_trainSet_X, 
        Y = m_trainSet_times, 
        W = m_trainSet_W,
        D = m_trainSet_events, 
        W.hat = as.vector(W_hat_train_adj2), 
        num.trees = n_trees_val, 
        target = "RMST", 
        horizon = horizon, 
#        honesty = TRUE,
#        min.node.size = 5,
#        alpha = 0.05,
#        imbalance.penalty = 0.1,
#        stabilize.splits = TRUE,
        tune.parameters = "all", 
        seed = 123
      )
      
      ate_train_RMST <- average_treatment_effect(csf_model_RMST)
      csf_pred_test_RMST <- predict(csf_model_RMST, m_testSet_X, 
                                  W.hat = as.vector(W_hat_test_adj2), 
                                  estimate.variance = TRUE)
      
      plot_effects_distribution(
        list(
            predictions = csf_pred_test_RMST$predictions,
            ate = ate_train_RMST
        ),
        "RMST",
        horizon
    )

      # Store results
      results_RMST[i, ] <- c(
        horizon,
        ate_train_RMST[1],
        ate_train_RMST[2],
        mean(csf_pred_test_RMST$predictions),
        mean(sqrt(csf_pred_test_RMST$variance.estimates))
      )
      
      
    }, error = function(e) {
      warning(sprintf("Error in horizon %d: %s", horizon, e$message))
    })
  }
  
  # Save results safely
  tryCatch({
    write.csv(results_SP, 
              file.path(output_dir, "train_and_test_ATE_SP.csv"), 
              row.names = FALSE)
    write.csv(results_RMST, 
              file.path(output_dir, "train_and_test_ATE_RMST.csv"), 
              row.names = FALSE)
  }, error = function(e) {
    warning("Error saving results: ", e$message)
  })
  
  # Create and save plots
  # [plotting code remains the same]
  
  return(list(results_SP = results_SP, results_RMST = results_RMST))
}


#########
# Novel method: Time-Varying Treatment Effects Framework
#########

# This function models the trajectory of treatment effects over time
# It captures the rising and diminishing patterns observed in radiotherapy outcomes

model_temporal_treatment_effects <- function(results_df, 
                                           method = c("parametric", "nonparametric"),
                                           output_dir = ".",
                                           result_type = "SP") {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE)
  
  # Add diagnostic information
  cat("Input dataframe dimensions:", dim(results_df), "\n")
  cat("Column names:", paste(colnames(results_df), collapse=", "), "\n")
  
  # Dynamically set the column names based on result_type
  estimate_col <- paste0("ATE_estimate_test_", result_type)
  se_col <- paste0("ATE_se_test_", result_type)
  
  # Check if the columns exist in the dataset
  if(!estimate_col %in% colnames(results_df)) {
    stop(paste("Column", estimate_col, "not found in results data frame."))
  }
  if(!se_col %in% colnames(results_df)) {
    stop(paste("Column", se_col, "not found in results data frame."))
  }
  
  # Print NA diagnostics
  cat("NA values in horizon_sel:", sum(is.na(results_df$horizon_sel)), "\n")
  cat("NA values in", estimate_col, ":", sum(is.na(results_df[[estimate_col]])), "\n")
  cat("NA values in", se_col, ":", sum(is.na(results_df[[se_col]])), "\n")
  
  # Filter out rows with NA values
  valid_df <- results_df[!is.na(results_df$horizon_sel) & 
                          !is.na(results_df[[estimate_col]]) & 
                          !is.na(results_df[[se_col]]), ]
  
  # Check if we still have data after filtering
  if(nrow(valid_df) == 0) {
    stop("No valid data remains after filtering NA values.")
  }
  
  cat("Filtered dataframe dimensions:", dim(valid_df), "\n")
  print(head(valid_df))
  
  # Extract horizons and estimates from valid data
  horizons <- valid_df$horizon_sel
  estimates <- valid_df[[estimate_col]]
  se_values <- valid_df[[se_col]]
  
  # Basic visualization of raw estimates
  raw_plot <- ggplot(valid_df, aes(x = horizon_sel, y = .data[[estimate_col]])) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = .data[[estimate_col]] - .data[[se_col]], 
                      ymax = .data[[estimate_col]] + .data[[se_col]]), width = 2) +
    labs(title = paste("Raw Treatment Effect Estimates Over Time -", result_type),
         x = "Months after treatment",
         y = paste("Treatment effect (", result_type, ")", sep = "")) +
    theme_bw() +
    theme(text = element_text(size = 12))
  
  ggsave(file.path(output_dir, paste0("raw_temporal_effects_", result_type, ".png")), raw_plot, width = 8, height = 6)
  
  # Method 1: Parametric approach - Fit quadratic model to capture rise and fall
  if(method[1] == "parametric" || length(method) > 1) {
    # Quadratic model to capture rise and fall with valid data
    quad_model <- lm(estimates ~ horizons + I(horizons^2), weights = 1/se_values^2)
    
    # Calculate peak timing
    peak_coefs <- coef(quad_model)
    
    # Make sure we have valid coefficients
    if(length(peak_coefs) < 3 || is.na(peak_coefs[2]) || is.na(peak_coefs[3]) || peak_coefs[3] == 0) {
      warning("Cannot calculate peak time: insufficient or invalid model coefficients")
      peak_time <- NA
      peak_effect <- NA
    } else {
      peak_time <- -peak_coefs[2] / (2 * peak_coefs[3])
      peak_effect <- predict(quad_model, newdata = data.frame(horizons = peak_time))
      
      # Check if peak_time is within range
      if(peak_time < min(horizons) || peak_time > max(horizons)) {
        warning("Calculated peak time is outside the range of observed horizons")
      }
    }
    
    # Calculate half-life with error handling
    if(!is.na(peak_time) && !is.na(peak_effect)) {
      predict_fn <- function(x) predict(quad_model, newdata = data.frame(horizons = x)) - (peak_effect/2)
      half_life <- try(uniroot(predict_fn, c(peak_time, max(horizons) + 50))$root - peak_time, silent = TRUE)
      
      if(inherits(half_life, "try-error")) {
        half_life <- NA
        warning("Could not calculate half-life, effect may not decline to half of peak within range")
      }
    } else {
      half_life <- NA
    }
    
    # Generate predictions across full time range
    pred_horizons <- seq(min(horizons), max(horizons), length.out = 100)
    pred_effects <- predict(quad_model, newdata = data.frame(horizons = pred_horizons), 
                          se.fit = TRUE)
    
    # Create prediction dataframe
    pred_df <- data.frame(
      horizons = pred_horizons,
      effect = pred_effects$fit,
      lower = pred_effects$fit - 1.96 * pred_effects$se.fit,
      upper = pred_effects$fit + 1.96 * pred_effects$se.fit
    )
    
    # Plot parametric model with error handling for NA values
    param_plot <- ggplot() +
      geom_ribbon(data = pred_df, aes(x = horizons, ymin = lower, ymax = upper), 
                 alpha = 0.2, fill = "blue") +
      geom_line(data = pred_df, aes(x = horizons, y = effect), color = "blue", linewidth = 1.2) +
      geom_point(data = valid_df, aes(x = horizon_sel, y = .data[[estimate_col]]), 
                size = 3) +
      geom_errorbar(data = valid_df, 
                   aes(x = horizon_sel, 
                       ymin = .data[[estimate_col]] - .data[[se_col]], 
                       ymax = .data[[estimate_col]] + .data[[se_col]]), 
                   width = 2)
    
    # Add peak line and annotation if peak is valid
    if(!is.na(peak_time) && peak_time >= min(horizons) && peak_time <= max(horizons)) {
      param_plot <- param_plot +
        geom_vline(xintercept = peak_time, linetype = "dashed", color = "red") +
        annotate("text", x = peak_time + 5, y = max(pred_df$upper), 
                label = paste("Peak at", round(peak_time, 1), "months"), 
                hjust = 0, color = "red")
    }
    
    # Add title and labels
    param_plot <- param_plot +
      labs(title = paste("Parametric Time-Varying Treatment Effect Model -", result_type),
           subtitle = if(!is.na(peak_time)) {
             paste("Peak effect at", round(peak_time, 1), "months, effect half-life:", 
                  ifelse(is.na(half_life), "Not reached", paste(round(half_life, 1), "months")))
           } else {
             "Peak effect could not be determined"
           },
           x = "Months after treatment",
           y = paste("Treatment effect (", result_type, ")", sep = "")) +
      theme_bw() +
      theme(text = element_text(size = 12))
    
    ggsave(file.path(output_dir, paste0("parametric_temporal_effects_", result_type, ".png")), 
           param_plot, width = 10, height = 7)
    
    # Save parametric model results
    param_results <- list(
      model = quad_model,
      summary = summary(quad_model),
      peak_time = peak_time,
      peak_effect = peak_effect,
      half_life = half_life,
      predictions = pred_df
    )
    
    saveRDS(param_results, file.path(output_dir, paste0("parametric_model_results_", result_type, ".rds")))
    write.csv(pred_df, file.path(output_dir, paste0("parametric_predictions_", result_type, ".csv")), row.names = FALSE)
  }
  
  # Method 2: Nonparametric approach using smoothing splines
  if(method[1] == "nonparametric" || length(method) > 1) {
    # Use weights based on standard errors
    weights <- 1/se_values^2
    
    # Fit smoothing spline with cross-validation for optimal smoothing using valid data
    tryCatch({
      spline_model <- smooth.spline(horizons, estimates, w = weights, cv = TRUE)
      
      # Generate predictions
      spline_pred <- predict(spline_model, seq(min(horizons), max(horizons), length.out = 100))
      spline_df <- data.frame(horizons = spline_pred$x, effect = spline_pred$y)
      
      # Calculate derivative to find peaks and inflection points
      spline_deriv <- predict(spline_model, spline_df$horizons, deriv = 1)
      spline_deriv2 <- predict(spline_model, spline_df$horizons, deriv = 2)
      
      spline_df$derivative <- spline_deriv$y
      spline_df$derivative2 <- spline_deriv2$y
      
      # Find where derivative crosses zero (peaks/troughs)
      sign_changes <- diff(sign(spline_df$derivative))
      peak_indices <- which(sign_changes < 0) + 1
      trough_indices <- which(sign_changes > 0) + 1
      
      # Only keep peaks/troughs with sufficient magnitude
      if(length(peak_indices) > 0) {
        peak_time_np <- spline_df$horizons[peak_indices[which.max(spline_df$effect[peak_indices])]]
        peak_effect_np <- max(spline_df$effect[peak_indices])
      } else {
        peak_time_np <- NA
        peak_effect_np <- NA
        warning("No peaks detected in nonparametric model")
      }
      
      # Calculate acceleration/deceleration phases
      inflection_indices <- which(diff(sign(spline_df$derivative2)) != 0) + 1
      
      # Plot nonparametric model
      nonparam_plot <- ggplot() +
        geom_line(data = spline_df, aes(x = horizons, y = effect), color = "blue", linewidth = 1.2) +
        geom_point(data = valid_df, aes(x = horizon_sel, y = .data[[estimate_col]]), size = 3) +
        geom_errorbar(data = valid_df, 
                     aes(x = horizon_sel, 
                         ymin = .data[[estimate_col]] - .data[[se_col]], 
                         ymax = .data[[estimate_col]] + .data[[se_col]]), width = 2)
      
      # Add peak point if found
      if(!is.na(peak_time_np)) {
        nonparam_plot <- nonparam_plot +
          geom_vline(xintercept = peak_time_np, linetype = "dashed", color = "red") +
          annotate("text", x = peak_time_np + 5, y = max(valid_df[[estimate_col]]), 
                  label = paste("Peak at", round(peak_time_np, 1), "months"), 
                  hjust = 0, color = "red")
      }
      
      # Add inflection points if found
      if(length(inflection_indices) > 0) {
        nonparam_plot <- nonparam_plot +
          geom_vline(xintercept = spline_df$horizons[inflection_indices], 
                    linetype = "dotted", color = "darkgreen", alpha = 0.7)
      }
      
      nonparam_plot <- nonparam_plot +
        labs(title = paste("Nonparametric Time-Varying Treatment Effect Model -", result_type),
             subtitle = ifelse(!is.na(peak_time_np), 
                              paste("Peak effect at", round(peak_time_np, 1), "months"), 
                              "No clear peak detected"),
             x = "Months after treatment",
             y = paste("Treatment effect (", result_type, ")", sep = "")) +
        theme_bw() +
        theme(text = element_text(size = 12))
      
      ggsave(file.path(output_dir, paste0("nonparametric_temporal_effects_", result_type, ".png")), 
             nonparam_plot, width = 10, height = 7)
      
      # Save nonparametric model results
      nonparam_results <- list(
        model = spline_model,
        peak_time = peak_time_np,
        peak_effect = peak_effect_np,
        inflection_points = if(exists("inflection_indices")) spline_df$horizons[inflection_indices] else NULL,
        predictions = spline_df
      )
      
      saveRDS(nonparam_results, file.path(output_dir, paste0("nonparametric_model_results_", result_type, ".rds")))
      write.csv(spline_df, file.path(output_dir, paste0("nonparametric_predictions_", result_type, ".csv")), row.names = FALSE)
      
    }, error = function(e) {
      warning("Error in nonparametric modeling: ", e$message)
      nonparam_results <- list(error = e$message)
      saveRDS(nonparam_results, file.path(output_dir, paste0("nonparametric_model_results_error_", result_type, ".rds")))
    })
  }
  
  # Create combined visualization only if both models were successful
  if(length(method) > 1 && exists("spline_df") && exists("pred_df") && 
     !is.na(peak_time) && !is.na(peak_time_np)) {
    
    combined_plot <- ggplot() +
      # Add nonparametric model
      geom_line(data = spline_df, aes(x = horizons, y = effect, color = "Nonparametric"), 
               linewidth = 1.2) +
      # Add parametric model
      geom_line(data = pred_df, aes(x = horizons, y = effect, color = "Parametric"), 
               linewidth = 1.2) +
      geom_ribbon(data = pred_df, aes(x = horizons, ymin = lower, ymax = upper), 
                 alpha = 0.1, fill = "blue") +
      # Add raw data points
      geom_point(data = valid_df, aes(x = horizon_sel, y = .data[[estimate_col]]), 
                size = 3, color = "black") +
      geom_errorbar(data = valid_df, 
                   aes(x = horizon_sel, 
                       ymin = .data[[estimate_col]] - .data[[se_col]], 
                       ymax = .data[[estimate_col]] + .data[[se_col]]), width = 2, color = "black")
    
    # Add peak lines if they exist and are within valid range
    if(!is.na(peak_time) && peak_time >= min(horizons) && peak_time <= max(horizons)) {
      combined_plot <- combined_plot +
        geom_vline(xintercept = peak_time, linetype = "dashed", color = "blue")
    }
    
    if(!is.na(peak_time_np) && peak_time_np >= min(horizons) && peak_time_np <= max(horizons)) {
      combined_plot <- combined_plot +
        geom_vline(xintercept = peak_time_np, linetype = "dashed", color = "red")
    }
    
    # Add labels and styling
    combined_plot <- combined_plot +
      scale_color_manual(values = c("Parametric" = "blue", "Nonparametric" = "red"),
                        name = "Model Type") +
      labs(title = paste("Comparison of Time-Varying Treatment Effect Models -", result_type),
           subtitle = paste("Parametric peak at", round(peak_time, 1), 
                           "months, Nonparametric peak at", 
                           round(peak_time_np, 1), "months"),
           x = "Months after treatment",
           y = paste("Treatment effect (", result_type, ")", sep = "")) +
      theme_bw() +
      theme(text = element_text(size = 12),
            legend.position = "bottom")
    
    ggsave(file.path(output_dir, paste0("combined_temporal_effects_", result_type, ".png")), 
           combined_plot, width = 12, height = 8)
  }
  
  # Return comprehensive results
  results <- list(
    raw_data = valid_df
  )
  
  if(method[1] == "parametric" || length(method) > 1) {
    if(exists("param_results")) {
      results$parametric <- param_results
    } else {
      results$parametric <- list(error = "Parametric modeling failed")
    }
  }
  
  if(method[1] == "nonparametric" || length(method) > 1) {
    if(exists("nonparam_results")) {
      results$nonparametric <- nonparam_results
    } else {
      results$nonparametric <- list(error = "Nonparametric modeling failed")
    }
  }
  
  return(results)
}


CSF_results <- implement_causal_forests(
  m_trainSet_X, m_trainSet_W, m_trainSet_times, m_trainSet_events,
  m_testSet_X, m_testSet_W, m_testSet_times, m_testSet_events,
  m_W_hat_train_adj2, m_W_hat_test_adj2,
  output_dir = "causal_forest_results"
)

print(CSF_results$results_SP)
print(CSF_results$results_RMST)

time_varying_results_SP <- model_temporal_treatment_effects(
  CSF_results$results_SP,
  method = c("parametric", "nonparametric"),
  output_dir = "time_varying_effects_results"
)

# Create summary of the temporal patterns
temporal_patterns_summary <- data.frame(
  Metric = "SP",
  Parametric_Peak_Time = time_varying_results_SP$parametric$peak_time,
  Parametric_Peak_Effect = time_varying_results_SP$parametric$peak_effect,
  Parametric_Half_Life = time_varying_results_SP$parametric$half_life,
  Nonparametric_Peak_Time = time_varying_results_SP$nonparametric$peak_time,
  Nonparametric_Peak_Effect = time_varying_results_SP$nonparametric$peak_effect
)

write.csv(temporal_patterns_summary, "time_varying_effects_summary.csv", row.names = FALSE)

# Compare to standard approach without time-varying effects modeling
print("Standard approach vs. Time-varying effects framework:")
print(paste("Standard approach identifies effects at discrete timepoints but cannot model the trajectory"))
print(paste("Time-varying framework identifies peak effect at", 
            round(time_varying_results_SP$parametric$peak_time, 1), 
            "months with half-life of", 
            round(time_varying_results_SP$parametric$half_life, 1), 
            "months"))


#########
# dummy outcome test

perform_dummy_outcome_tests <- function(trainSet_X, trainSet_W, trainSet_times, trainSet_events, 
                                        W_hat_train_adj2,  
                                        num_repetitions = 20, n_trees_val = 5000, 
                                        seed = 1234, output_dir = ".") {
  
  set.seed(seed)
  dir.create(output_dir, showWarnings = FALSE)
  
  # Define horizons inside the function
  horizons <- seq(12, 120, by=12)
  
  # Function to perform dummy test for a specific target
  perform_dummy_test <- function(target) {
    dummy_results <- data.frame(
      Repetition = integer(),
      Horizon = integer(),
      CATE_Estimate = numeric(),
      Standard_Error = numeric()
    )
    
    for (rep in 1:num_repetitions) {
      trainSet_W_Dummy <- sample(trainSet_W, length(trainSet_W), replace=TRUE)
      train_DummyOutcomeTimes <- sample(trainSet_times, length(trainSet_times), replace=TRUE)
      
      for (horizon in horizons) {
        forest_dummy <- causal_survival_forest(X = trainSet_X, 
                                               Y = train_DummyOutcomeTimes, 
                                               W = trainSet_W_Dummy, 
                                               W.hat = as.vector(W_hat_train_adj2),
                                               D = trainSet_events, 
                                               num.trees = n_trees_val, 
                                               target = target, 
                                               horizon = horizon, 
                                               seed = seed)
        ate_dummy <- average_treatment_effect(forest_dummy)
        
        if (is.list(ate_dummy)) {
          cate_estimate <- ate_dummy$estimate 
          standard_error <- ate_dummy$std.err
        } else if (is.vector(ate_dummy)) {
          cate_estimate <- ate_dummy[1] 
          standard_error <- ate_dummy[2]
        }
        
        dummy_results <- rbind(dummy_results, data.frame(
          Repetition = rep,
          Horizon = horizon,
          CATE_Estimate = cate_estimate,
          Standard_Error = standard_error
        ))
      }
    }
    
    return(dummy_results)
  }
  
  # tests for both targets
  dummy_results_SP <- perform_dummy_test("survival.probability")
  dummy_results_RMST <- perform_dummy_test("RMST")
  
  write.csv(dummy_results_SP, file.path(output_dir, "Dummy_Outcome_SP_Results.csv"), row.names = FALSE)
  write.csv(dummy_results_RMST, file.path(output_dir, "Dummy_Outcome_RMST_Results.csv"), row.names = FALSE)
  
  create_boxplot <- function(data, y_label, filename) {
    data$Horizon <- as.factor(data$Horizon)
    p <- ggplot(data, aes(x = Horizon, y = CATE_Estimate, fill = Horizon)) +
      geom_boxplot() +
      labs(title = "Boxplot of CATE Estimates with Dummy Outcome",
           x = "Survival time (months)",
           y = y_label) +
      theme_bw() +
      scale_fill_brewer(palette = "Set3") +
      theme(legend.position = "none",
            axis.text = element_text(face = "bold", size = 12),
            axis.title = element_text(face = "bold", size = 14),
            plot.title = element_text(face = "bold", size = 16))
    
    ggsave(file.path(output_dir, filename), p, width = 10, height = 8, dpi = 600)
  }
  
  create_boxplot(dummy_results_SP, "Causal effect (SP)", "Dummy_Outcome_SP_Boxplot.png")
  create_boxplot(dummy_results_RMST, "Causal effect (RMST, months)", "Dummy_Outcome_RMST_Boxplot.png")
  
  return(list(SP_results = dummy_results_SP, RMST_results = dummy_results_RMST))
}

dummy_results <- perform_dummy_outcome_tests(
  trainSet_X = trainSet_X,
  trainSet_W = trainSet_W,
  trainSet_times = trainSet_times,
  trainSet_events = trainSet_events,
  W_hat_train_adj2 = W_hat_train_adj2,
  output_dir = "dummy_outcome_results"
)

########
# refutation tests with fake confounders

refutation_fake_confounder_tests <- function(trainSet_X, trainSet_W, trainSet_times, trainSet_events,
                                          W_hat_train_adj2, testSet_X, testSet_W, testSet_times, 
                                          testSet_events, W_hat_test_adj2,
                                          num_repetitions = 20, n_trees_val = 5000,
                                          confounder_strength = c(0.1, 0.3, 0.5),
                                          seed = 1234, output_dir = ".") {
  
  set.seed(seed)
  dir.create(output_dir, showWarnings = FALSE)
  
  horizons <- seq(12, 120, by=12)
  
  generate_fake_confounder <- function(X, strength) {
    n <- nrow(X)
    # Create confounder correlated with covariates
    Z <- scale(as.matrix(X)) %*% rnorm(ncol(X)) * strength + rnorm(n) * (1 - strength)
    return(scale(Z)) # Standardize confounder
  }
  
  results <- data.frame(
    Repetition = integer(),
    Horizon = integer(),
    Strength = numeric(),
    Target = character(),
    Original_ATE = numeric(),
    Confounded_ATE = numeric(),
    ATE_Difference = numeric(),
    Original_SE = numeric(),
    Confounded_SE = numeric()
  )
  
  for(strength in confounder_strength) {
    for(rep in 1:num_repetitions) {
      Z_train <- generate_fake_confounder(trainSet_X, strength)
      Z_test <- generate_fake_confounder(testSet_X, strength)
      
      # fake confounder to feature matrices
      trainSet_X_augmented <- cbind(trainSet_X, fake_confounder = Z_train)
      testSet_X_augmented <- cbind(testSet_X, fake_confounder = Z_test)
      
      for(horizon in horizons) {
        # survival probability and RMST targets
        for(target in c("survival.probability", "RMST")) {

          original_forest <- causal_survival_forest(
            X = trainSet_X,
            Y = trainSet_times,
            W = trainSet_W,
            D = trainSet_events,
            W.hat = as.vector(W_hat_train_adj2),
            num.trees = n_trees_val,
            target = target,
            horizon = horizon,
            seed = seed + rep
          )
          
          confounded_forest <- causal_survival_forest(
            X = trainSet_X_augmented,
            Y = trainSet_times,
            W = trainSet_W,
            D = trainSet_events,
            W.hat = as.vector(W_hat_train_adj2),
            num.trees = n_trees_val,
            target = target,
            horizon = horizon,
            seed = seed + rep
          )
          
          pred_original <- predict(original_forest, testSet_X)
          pred_confounded <- predict(confounded_forest, testSet_X_augmented)
          
          ate_original <- mean(pred_original$predictions)
          ate_confounded <- mean(pred_confounded$predictions)
          
          se_original <- sqrt(mean(pred_original$variance.estimates))
          se_confounded <- sqrt(mean(pred_confounded$variance.estimates))
          
          results <- rbind(results, data.frame(
            Repetition = rep,
            Horizon = horizon,
            Strength = strength,
            Target = target,
            Original_ATE = ate_original,
            Confounded_ATE = ate_confounded,
            ATE_Difference = abs(ate_original - ate_confounded),
            Original_SE = se_original,
            Confounded_SE = se_confounded
          ))
        }
      }
      
      cat(sprintf("Completed repetition %d for strength %.1f\n", rep, strength))
    }
  }
  
  write.csv(results, file.path(output_dir, "fake_confounder_results.csv"), row.names = FALSE)
  
  for(target_type in c("survival.probability", "RMST")) {
    target_data <- subset(results, Target == target_type)
    
    p <- ggplot(target_data, aes(x = as.factor(Horizon), y = ATE_Difference, fill = as.factor(Strength))) +
      geom_boxplot() +
      facet_wrap(~Strength, labeller = label_both) +
      labs(title = paste("ATE Differences with Fake Confounders -", target_type),
           x = "Horizon",
           y = "Absolute ATE Difference",
           fill = "Confounder Strength") +
      theme_bw() +
      theme(
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 14),
        plot.title = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(face = "bold", size = 10)
      )
    
    ggsave(
      file.path(output_dir, paste0("fake_confounder_effects_", tolower(target_type), ".png")),
      p,
      width = 12,
      height = 8,
      dpi = 600
    )
  }
  
  summary_stats <- results %>%
    group_by(Strength, Target) %>%
    summarise(
      Mean_ATE_Diff = mean(ATE_Difference),
      SD_ATE_Diff = sd(ATE_Difference),
      Max_ATE_Diff = max(ATE_Difference),
      Mean_Original_SE = mean(Original_SE),
      Mean_Confounded_SE = mean(Confounded_SE)
    )
  
  write.csv(summary_stats, file.path(output_dir, "fake_confounder_summary.csv"), row.names = FALSE)
  
  return(list(
    detailed_results = results,
    summary_stats = summary_stats
  ))
}

# refutation tests with fake confounders
refutation_results <- refutation_fake_confounder_tests(
  trainSet_X = trainSet_X,
  trainSet_W = trainSet_W,
  trainSet_times = trainSet_times,
  trainSet_events = trainSet_events,
  W_hat_train_adj2 = W_hat_train_adj2,
  testSet_X = testSet_X,
  testSet_W = testSet_W,
  testSet_times = testSet_times,
  testSet_events = testSet_events,
  W_hat_test_adj2 = W_hat_test_adj2,
  # You can reduce/add if necessary b/c right now it is taking too long to run
  num_repetitions = 20,
  output_dir = "refutation_test_results"
)

print("Summary of refutation test results:")
print(refutation_results$summary_stats)

sp_results <- dummy_results$SP_results
rmst_results <- dummy_results$RMST_results

print(summary(sp_results$CATE_Estimate))
print(summary(rmst_results$CATE_Estimate))

###########

horizons <- seq(12, 120, by=12)
n_trees_val <- 5000

# Negative control with an irrelevant causal variable: causal effect should be zero at all times
set.seed(123)
binary_vector <- sample(c(0, 1), size = nrow(trainSet_X), replace = TRUE)

trainSet_W_neg_c <- binary_vector

trainSet_X_neg_c <- trainSet_X

HN_cate_results_RMST_train_neg_c <- data.frame(Horizon = integer(), Estimate = numeric(), Standard_Error = numeric())

for (h in horizons) {
  forest_h_RMST_train_neg_c <- causal_survival_forest(X = trainSet_X_neg_c, Y = trainSet_times, W = trainSet_W_neg_c, W.hat = as.vector(W_hat_train_adj2),
                                                      D = trainSet_events, num.trees = n_trees_val, target = "RMST", horizon = h, tune.parameters = "all", seed = 123)
  ate_h_RMST_train_neg_c <- average_treatment_effect(forest_h_RMST_train_neg_c)
  HN_cate_results_RMST_train_neg_c <- rbind(HN_cate_results_RMST_train_neg_c, data.frame(Horizon = h, Estimate = ate_h_RMST_train_neg_c[1], Standard_Error = ate_h_RMST_train_neg_c[2]))
}

print(HN_cate_results_RMST_train_neg_c)
write.csv(HN_cate_results_RMST_train_neg_c, file = "HN_cate_results_RMST_train_neg_c.csv", row.names = TRUE)

# Same for survival probabilities with negative control
HN_cate_results_RMST_train_neg_c_surv <- data.frame(Horizon = integer(), Estimate = numeric(), Standard_Error = numeric())

for (h in horizons) {
  forest_h_RMST_train_neg_c_surv <- causal_survival_forest(X = trainSet_X_neg_c, Y = trainSet_times, W = trainSet_W_neg_c, W.hat = as.vector(W_hat_train_adj2),
                                                           D = trainSet_events, num.trees = n_trees_val, target = "survival.probability", horizon = h, 
                                                           tune.parameters = "all",  seed = 123)
  ate_h_RMST_train_neg_c_surv <- average_treatment_effect(forest_h_RMST_train_neg_c_surv)
  HN_cate_results_RMST_train_neg_c_surv <- rbind(HN_cate_results_RMST_train_neg_c_surv, data.frame(Horizon = h, Estimate = ate_h_RMST_train_neg_c_surv[1], Standard_Error = ate_h_RMST_train_neg_c_surv[2]))
}

print(HN_cate_results_RMST_train_neg_c_surv)
write.csv(HN_cate_results_RMST_train_neg_c_surv, file = "HN_cate_results_RMST_train_neg_c_surv.csv", row.names = TRUE)


#######
# Generating several noise variables as another test: they should not affect the causal effect much
trainSet_X_noise <- trainSet_X
for (i in 1:5) {
  trainSet_X_noise[paste("NoiseVar_", i)] <- rnorm(nrow(trainSet_X_noise))
}

# for RMST
results_noise <- matrix(0, length(horizons), 2)

for (i in 1:length(horizons)) {
  forest_with_noise <- causal_survival_forest(X = trainSet_X_noise, Y = trainSet_times, W = trainSet_W, W.hat = as.vector(W_hat_train_adj2),
                                              D = trainSet_events, num.trees = n_trees_val, target = "RMST", horizon = horizons[i], seed = 123)
  ate_with_noise <- average_treatment_effect(forest_with_noise)
  
  if (is.list(ate_with_noise)) {
    results_noise[i, 1] <- ate_with_noise$estimate
    results_noise[i, 2] <- ate_with_noise$std.err
  } else if (is.vector(ate_with_noise)) {
    results_noise[i, 1] <- ate_with_noise[1]
    results_noise[i, 2] <- ate_with_noise[2]
  }
}

results_noise_df <- data.frame(Horizon = horizons, CATE_Estimate_with_noise = results_noise[, 1], Standard_Error = results_noise[, 2])
write.csv(results_noise_df, file = "results_df_noise.csv", row.names = TRUE)


# for SP
results_noise_SP <- matrix(0, length(horizons), 2)

for (i in 1:length(horizons)) {
  forest_with_noise_SP <- causal_survival_forest(X = trainSet_X_noise, Y = trainSet_times, W = trainSet_W, W.hat = as.vector(W_hat_train_adj2),
                                              D = trainSet_events, num.trees = n_trees_val, target = "survival.probability", horizon = horizons[i], seed = 123)
  ate_with_noise_SP <- average_treatment_effect(forest_with_noise_SP)
  
  if (is.list(ate_with_noise_SP)) {
    results_noise_SP[i, 1] <- ate_with_noise_SP$estimate
    results_noise_SP[i, 2] <- ate_with_noise_SP$std.err
  } else if (is.vector(ate_with_noise_SP)) {
    results_noise_SP[i, 1] <- ate_with_noise_SP[1]
    results_noise_SP[i, 2] <- ate_with_noise_SP[2]
  }
}

results_noise_SP_df <- data.frame(Horizon = horizons, CATE_Estimate_with_noise_SP = results_noise_SP[, 1], Standard_Error = results_noise_SP[, 2])
write.csv(results_noise_SP_df, file = "results_df_noise_SP.csv", row.names = TRUE)

# noise should not matter much

#####

horizon_sel <- 3*12
forest_sel <- causal_survival_forest(X = trainSet_X, Y = trainSet_times, W = trainSet_W, W.hat = as.vector(W_hat_train_adj2), D = trainSet_events, num.trees = n_trees_val, target = "survival.probability", horizon = horizon_sel, tune.parameters = "all", seed=1234)
forest_sel_preds_train <- predict(forest_sel, trainSet_X, estimate.variance = TRUE)
forest_sel_preds_test <- predict(forest_sel, testSet_X, estimate.variance = TRUE)

forest_sel_preds_train_df <- cbind(trainSet_X, forest_sel_preds_train)
forest_sel_preds_test_df <- cbind(testSet_X, forest_sel_preds_test)

write.csv(forest_sel_preds_train_df, file = "HN_forest_sel_preds_train_df.csv", row.names = FALSE)
write.csv(forest_sel_preds_test_df, file = "HN_forest_sel_preds_test_df.csv", row.names = FALSE)


pfun <- function(object, newdata) {
  predict(object, newdata = newdata, estimate.variance = TRUE)$predictions
}

nsim_shap <- 1000

shap <- fastshap::explain(forest_sel, X = trainSet_X, pred_wrapper = pfun, nsim = nsim_shap)
colnames(shap) <- paste0(colnames(shap), "_SHAP")
shap_vals <- cbind(trainSet_X, shap, forest_sel_preds_train)

average_prediction <- mean(shap_vals$predictions)

# normalize the SHAP values so the Monte Carlo random error is removed as much as possible
# Step 1: Sum all the columns ending in _SHAP
shap_cols <- grep("_SHAP$", names(shap_vals), value = TRUE)
shap_vals$sum_SHAP <- rowSums(shap_vals[, shap_cols])

# Step 2: Normalize each _SHAP column
for (col in shap_cols) {
  norm_col <- paste0(col, "_norm_SHAP")
  shap_vals[[norm_col]] <- (shap_vals[[col]] / shap_vals$sum_SHAP) * (shap_vals$predictions - average_prediction)
}

# Step 3: sum of the normalized SHAP columns
norm_shap_cols <- grep("_norm_SHAP$", names(shap_vals), value = TRUE)
shap_vals$sum_norm_SHAP <- rowSums(shap_vals[, norm_shap_cols])

# Check if sum of normalized SHAP columns equals the difference between predictions and average prediction
all.equal(shap_vals$sum_norm_SHAP, shap_vals$predictions - average_prediction)

write.csv(shap_vals, file = "HN_shap_vals_raw.csv", row.names = FALSE)

cols_to_remove <- grep("_SHAP$", names(shap_vals), value = TRUE)
cols_to_keep <- grep("_SHAP_norm_SHAP$", names(shap_vals), value = TRUE)
cols_to_remove <- setdiff(cols_to_remove, cols_to_keep)

shap_vals <- shap_vals[, !(names(shap_vals) %in% cols_to_remove)]

cols_to_rename <- grep("_SHAP_norm_SHAP$", names(shap_vals), value = TRUE)
new_names <- sub("_SHAP_norm_SHAP$", "_SHAP", cols_to_rename)
names(shap_vals)[names(shap_vals) %in% cols_to_rename] <- new_names
str(shap_vals)

shap_cols <- grep("_SHAP$", names(shap_vals), value = TRUE)

shap_sum <- rowSums(shap_vals[, shap_cols])
is_equal <- all.equal(shap_sum, shap_vals$predictions, tolerance = 1e-6)

if (isTRUE(is_equal)) {
  print("The sum of all _SHAP columns equals the predictions column.")
} else {
  print("The sum of all _SHAP columns does NOT equal the predictions column.")
  
  max_diff <- max(abs(shap_sum - shap_vals$predictions))
  print(paste("Maximum difference:", max_diff))
  
  diff_summary <- summary(shap_sum - shap_vals$predictions)
  print("Summary of differences:")
  print(diff_summary)
}


write.csv(shap_vals, file = "HN_shap_vals.csv", row.names = FALSE)

# plot age effects for interest

#emf("Age_SHAP_plot.emf", width = 7, height = 7, bg = "white", fg = "black", pointsize = 12, family = "Arial", coordDPI = 600)
png("Age_SHAP_plot.png", width = 7, height = 7, units = "in", res = 500, bg = "white", pointsize = 12)

plot(shap_vals$Age, shap_vals$Age_SHAP, xlab="Age (years)", ylab="SHAP value")

dev.off()

ordered_shap_vals <- shap_vals[order(shap_vals$predictions), ]
plot_data <- data.frame(
  Patients = 1:nrow(ordered_shap_vals),
  predictions = ordered_shap_vals$predictions,
  lower = ordered_shap_vals$predictions - 1.96 * sqrt(ordered_shap_vals$variance.estimates),
  upper = ordered_shap_vals$predictions + 1.96 * sqrt(ordered_shap_vals$variance.estimates)
)

#emf("causal_effects_plot.emf", width = 7, height = 7, bg = "white", fg = "black", pointsize = 12, family = "Arial", coordDPI = 600)

png("causal_effects_plot.png", width = 7, height = 7, units = "in", res = 500, bg = "white", pointsize = 12)
ggplot(plot_data, aes(x = Patients, y = predictions)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(x = "Patients", y = "Causal effect (SP)")

dev.off()

#emf("causal_effects_histogram.emf", width = 7, height = 7, bg = "white", fg = "black", pointsize = 12, family = "Arial", coordDPI = 600)

png("causal_effects_histogram.png", width = 7, height = 7, units = "in", res = 500, bg = "white", pointsize = 12)
hist(ordered_shap_vals$predictions, xlab = "Causal effect (SP)")

dev.off()

shap_cols <- grep("_SHAP$", names(shap_vals), value = TRUE)
medians <- apply(shap_vals[shap_cols], 2, median)
abs_medians <- apply(abs(shap_vals[shap_cols]), 2, median)

shap_vals_medians <- data.frame(medians, abs_medians)
write.csv(shap_vals_medians, file = "HN_shap_vals_medians.csv", row.names = TRUE)

shap_vals_medians <- shap_vals_medians[order(-shap_vals_medians$abs_medians), ]
shap_vals_medians_top_10 <- head(shap_vals_medians, 10)
write.csv(shap_vals_medians_top_10, file = "HN_shap_vals_medians_top_10.csv", row.names = TRUE)

std_devs <- apply(shap_vals, 2, sd)
cols <- names(std_devs)[std_devs > 0]
selected_shap_vals <- shap_vals[, cols]

varscor_selected_shap_vals <- corr.test(selected_shap_vals, method = "spearman", adjust = "bonf", alpha = .05, ci = FALSE)
varscor_selected_shap_vals_p <- varscor_selected_shap_vals$p

write.csv(varscor_selected_shap_vals_p, file = "varscor_selected_shap_vals_spearman_p.csv", row.names = TRUE)
write.csv(varscor_selected_shap_vals$r, file = "varscor_selected_shap_vals_spearman_r.csv", row.names = TRUE)

# Spearman correlation matrix
#emf("correlation_matrix_spearman_SHAP.emf", width = 7, height = 7, bg = "white", fg = "black", pointsize = 8, family = "Arial", coordDPI = 600)

png("correlation_matrix_spearman_SHAP.png", width = 7, height = 7, units = "in", res = 500, bg = "white", pointsize = 8)
corrplot(varscor_selected_shap_vals$r, p.mat = varscor_selected_shap_vals$p, method = 'circle', tl.col = "black", type = "upper", sig.level = 0.05, pch.cex = 0.6, cl.cex = 1, tl.cex = 1, insig = 'pch', pch = 19, pch.col = "white", diag = FALSE, font = 1)

dev.off()

# Pearson correlation matrix
varscor_selected_shap_vals_pearson <- corr.test(selected_shap_vals, method = "pearson", adjust = "bonf", alpha = .05, ci = FALSE)
varscor_selected_shap_vals_pearson_p <- varscor_selected_shap_vals_pearson$p

write.csv(varscor_selected_shap_vals_pearson_p, file = "varscor_selected_shap_vals_pearson_p.csv", row.names = TRUE)
write.csv(varscor_selected_shap_vals_pearson$r, file = "varscor_selected_shap_vals_pearson_r.csv", row.names = TRUE)

#emf("correlation_matrix_pearson_SHAP.emf", width = 7, height = 7, bg = "white", fg = "black", pointsize = 8, family = "Arial", coordDPI = 600)

png("correlation_matrix_pearson_SHAP.png", width = 7, height = 7, units = "in", res = 500, bg = "white", pointsize = 8)
corrplot(varscor_selected_shap_vals_pearson$r, p.mat = varscor_selected_shap_vals_pearson$p, method = 'circle', tl.col = "black", type = "upper", sig.level = 0.05, pch.cex = 0.6, cl.cex = 1, tl.cex = 1, insig = 'pch', pch = 19, pch.col = "white", diag = FALSE, font = 1)

dev.off()

SHAP_data_1 <- as.data.frame(subset(shap_vals,
                                    select=-c(predictions, variance.estimates)))

features <- names(SHAP_data_1)[!grepl("_SHAP$", names(SHAP_data_1))]
causal_effects <- shap_vals$predictions

for (feature in features) {

  # 5th and 95th percentiles of SHAP values for this feature
  shap_col <- paste0(feature, "_SHAP")
  y_min <- quantile(SHAP_data_1[[shap_col]], 0.05)
  y_max <- quantile(SHAP_data_1[[shap_col]], 0.95)
  
 p <- ggplot(SHAP_data_1, aes(x = .data[[feature]], y = .data[[paste0(feature, "_SHAP")]])) +
    geom_point(aes(color = causal_effects, shape = factor(Sex)), size = 3) +
    scale_shape_manual(values = c(16, 15)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(causal_effects)) +
    coord_cartesian(ylim = c(y_min, y_max)) +  # Set y-axis limits using coord_cartesian
    labs(title = paste("SHAP values for", feature),
         subtitle = "Y-axis limited to 5th-95th percentiles",
         x = feature,
         y = paste0(feature, "_SHAP"),
         color = "Causal Effects",
         shape = "Sex") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "darkgray"),
      axis.title = element_text(face = "bold", size = 12),
      legend.title = element_text(face = "bold", size = 10)
    )
  
  # emf_file <- paste0(make.names(feature), "_SHAP.emf")
  # emf(file = emf_file)

  png_file <- paste0(make.names(feature), "_SHAP.png")
  png(file = png_file, width = 7, height = 7, units = "in", res = 500, bg = "white", pointsize = 12)
  print(p)
  dev.off()
  
  # Print the percentages of points included/excluded for documentation
  total_points <- nrow(SHAP_data_1)
  points_in_range <- sum(SHAP_data_1[[shap_col]] >= y_min & SHAP_data_1[[shap_col]] <= y_max)
  percent_included <- (points_in_range/total_points) * 100
  cat(sprintf("Feature %s: %.1f%% of points shown (y-axis limits: %.3f to %.3f)\n", 
              feature, percent_included, y_min, y_max))
}

