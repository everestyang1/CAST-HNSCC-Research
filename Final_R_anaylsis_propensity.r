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
    
    p <- ggplot(effects_df, aes(x = effects)) +
        geom_histogram(bins = 30, fill = "lightblue", color = "black") +
        geom_vline(xintercept = results$ate[1], color = "red", linetype = "dashed") +
        theme_bw() +
        labs(title = sprintf("Distribution of Treatment Effects at %d Year Horizon (%s)", 
                            horizon, method),
             x = "Effect on Survival Probability",
             y = "Count")
    
    ggsave(sprintf("treatment_effects_%s_%dy.png", method, horizon), p, 
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
  
  # Generate covariates with more realistic distributions
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
  for(var in c("Age" ,"Stage" ,"HPV", "Sex", "Chemotherapy")) {
    mean_control <- mean(data[[var]][data$Cause == 0])
    mean_treated <- mean(data[[var]][data$Cause == 1])
    cat(sprintf("%s - Control: %.2f , Treated: %.2f\n", var ,mean_control ,mean_treated))
  }
  
  return(data)
}


# call the simulation function

simulated_data <- simulate_simple_data(n_samples = 3000, true_effect = 6)
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
prepared_data <- prepare_cancer_data(simulated_data)
# or read real data instead

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
    
  trainSet_with_scores_filtered <- subset(trainSet_with_scores, W_hat >= 0.05 & W_hat <= 0.95)
  testSet_with_scores_filtered <- subset(testSet_with_scores, W_hat >= 0.05 & W_hat <= 0.95)
  
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

implement_causal_forests <- function(trainSet_X, trainSet_W, trainSet_times, trainSet_events,
                                     testSet_X, testSet_W, testSet_times, testSet_events,
                                     W_hat_train_adj2, W_hat_test_adj2,
                                     n_trees_val = 3000, horizons = seq(10, 50, by=5),
                                     output_dir = ".") {
  
  # Ensure output directory exists
  dir.create(output_dir, showWarnings = FALSE)
  
  results_SP <- data.frame(horizon_sel = integer(), 
                           ATE_estimate_train_SP = numeric(), ATE_se_train_SP = numeric(),
                           ATE_estimate_test_SP = numeric(), ATE_se_test_SP = numeric())
  
  results_RMST <- data.frame(horizon_sel = integer(), 
                             ATE_estimate_train_RMST = numeric(), ATE_se_train_RMST = numeric(),
                             ATE_estimate_test_RMST = numeric(), ATE_se_test_RMST = numeric())
  
  for (horizon in horizons) {
    # Survival Probability
    csf_model_SP <- causal_survival_forest(X = trainSet_X, Y = trainSet_times, W = trainSet_W,
                                           D = trainSet_events, W.hat = as.vector(W_hat_train_adj2), 
                                           num.trees = n_trees_val, target = "survival.probability", 
                                           horizon = horizon, honesty = TRUE,
                                           min.node.size = 5,
                                           alpha = 0.05,
                                           imbalance.penalty = 0.1,  # Added imbalance penalty
                                           stabilize.splits = TRUE,  # Added split stabilization
                                           seed = 123)
    
    ate_train_SP <- average_treatment_effect(csf_model_SP)
    csf_pred_test_SP <- predict(csf_model_SP, testSet_X, W.hat = as.vector(W_hat_test_adj2), estimate.variance = TRUE)
    ate_h_SP_test <- mean(csf_pred_test_SP$predictions)
    ate_h_SP_test_sd <- mean(sqrt(csf_pred_test_SP$variance.estimates))
    
    plot_effects_distribution(csf_model_SP, "propensity_SP", horizon)
    
    for(var in c("Age", "Stage", "HPV")) {
      plot_treatment_effects_by_covariate(
        csf_pred_test_SP, 
        testSet_X, 
        var,
        title = paste("SP - Horizon", horizon)
      )
    }

    results_SP <- rbind(results_SP, data.frame(horizon_sel = horizon, 
                                               ATE_estimate_train_SP = ate_train_SP[1], ATE_se_train_SP = ate_train_SP[2],
                                               ATE_estimate_test_SP = ate_h_SP_test, ATE_se_test_SP = ate_h_SP_test_sd))
    
    # RMST
    csf_model_RMST <- causal_survival_forest(X = trainSet_X, Y = trainSet_times, W = trainSet_W,
                                             D = trainSet_events, W.hat = as.vector(W_hat_train_adj2), 
                                             num.trees = n_trees_val, target = "RMST", 
                                             horizon = horizon, honesty = TRUE,
                                             min.node.size = 5,
                                             alpha = 0.05,
                                             imbalance.penalty = 0.1,  # Added imbalance penalty
                                             stabilize.splits = TRUE,  # Added split stabilization
                                             seed = 123)
    
    ate_train_RMST <- average_treatment_effect(csf_model_RMST)
    csf_pred_test_RMST <- predict(csf_model_RMST, testSet_X, W.hat = as.vector(W_hat_test_adj2), estimate.variance = TRUE)
    ate_h_RMST_test <- mean(csf_pred_test_RMST$predictions)
    ate_h_RMST_test_sd <- mean(sqrt(csf_pred_test_RMST$variance.estimates))

    plot_effects_distribution(csf_model_RMST, "propensity_RMST", horizon)
    
    for(var in c("Age", "Stage", "HPV")) {
      plot_treatment_effects_by_covariate(
        csf_pred_test_RMST, 
        testSet_X, 
        var,
        title = paste("RMST - Horizon", horizon)
      )
    }
    
    results_RMST <- rbind(results_RMST, data.frame(horizon_sel = horizon, 
                                                   ATE_estimate_train_RMST = ate_train_RMST[1], ATE_se_train_RMST = ate_train_RMST[2],
                                                   ATE_estimate_test_RMST = ate_h_RMST_test, ATE_se_test_RMST = ate_h_RMST_test_sd))
  }
  
  write.csv(results_SP, file.path(output_dir, "train_and_test_ATE_SP.csv"), row.names = FALSE)
  write.csv(results_RMST, file.path(output_dir, "train_and_test_ATE_RMST.csv"), row.names = FALSE)
  
  plot_SP <- ggplot(results_SP, aes(x = horizon_sel)) +
    geom_point(aes(y = ATE_estimate_train_SP, color = "Train"), size = 5, position = position_nudge(x = -0.2)) +
    geom_errorbar(aes(ymin = ATE_estimate_train_SP - 1.96 * ATE_se_train_SP, 
                      ymax = ATE_estimate_train_SP + 1.96 * ATE_se_train_SP, color = "Train"), 
                  width = 0.5, size = 1, position = position_nudge(x = -0.2)) +
    geom_point(aes(y = ATE_estimate_test_SP, color = "Test"), size = 5, position = position_nudge(x = 0.2)) +
    geom_errorbar(aes(ymin = ATE_estimate_test_SP - 1.96 * ATE_se_test_SP, 
                      ymax = ATE_estimate_test_SP + 1.96 * ATE_se_test_SP, color = "Test"), 
                  width = 0.5, size = 1, position = position_nudge(x = 0.2)) +
    labs(title = "ATE for Survival Probability", x = "Horizon", y = "ATE") +
    theme_bw(base_size = 16) +
    theme(
      axis.text = element_text(color = "black", face = "bold", size = 14),
      axis.title = element_text(color = "black", face = "bold", size = 16),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 14)
    ) +
    scale_color_manual(values = c("Train" = "blue", "Test" = "red"))
  
  plot_RMST <- ggplot(results_RMST, aes(x = horizon_sel)) +
    geom_point(aes(y = ATE_estimate_train_RMST, color = "Train"), size = 5, position = position_nudge(x = -0.2)) +
    geom_errorbar(aes(ymin = ATE_estimate_train_RMST - 1.96 * ATE_se_train_RMST, 
                      ymax = ATE_estimate_train_RMST + 1.96 * ATE_se_train_RMST, color = "Train"), 
                  width = 0.5, size = 1, position = position_nudge(x = -0.2)) +
    geom_point(aes(y = ATE_estimate_test_RMST, color = "Test"), size = 5, position = position_nudge(x = 0.2)) +
    geom_errorbar(aes(ymin = ATE_estimate_test_RMST - 1.96 * ATE_se_test_RMST, 
                      ymax = ATE_estimate_test_RMST + 1.96 * ATE_se_test_RMST, color = "Test"), 
                  width = 0.5, size = 1, position = position_nudge(x = 0.2)) +
    labs(title = "ATE for RMST", x = "Horizon", y = "ATE") +
    theme_bw(base_size = 16) +
    theme(
      axis.text = element_text(color = "black", face = "bold", size = 14),
      axis.title = element_text(color = "black", face = "bold", size = 16),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 14)
    ) +
    scale_color_manual(values = c("Train" = "blue", "Test" = "red"))
  
  # high-resolution plots
  ggsave(file.path(output_dir, "ATE_SP_plot.png"), plot_SP, width = 12, height = 8, dpi = 600, units = "in")
  ggsave(file.path(output_dir, "ATE_RMST_plot.png"), plot_RMST, width = 12, height = 8, dpi = 600, units = "in")
  
  
  return(list(results_SP = results_SP, results_RMST = results_RMST))
}

CSF_results <- implement_causal_forests(
  trainSet_X, trainSet_W, trainSet_times, trainSet_events,
  testSet_X, testSet_W, testSet_times, testSet_events,
  W_hat_train_adj2, W_hat_test_adj2,
  output_dir = "causal_forest_results"
)

print(CSF_results$results_SP)
print(CSF_results$results_RMST)

#########
# dummy outcome test

perform_dummy_outcome_tests <- function(trainSet_X, trainSet_W, trainSet_times, trainSet_events, 
                                        W_hat_train_adj2,  
                                        num_repetitions = 20, n_trees_val = 3000, 
                                        seed = 1234, output_dir = ".") {
  
  set.seed(seed)
  dir.create(output_dir, showWarnings = FALSE)
  
  # Define horizons inside the function
  horizons <- seq(10, 50, by=5)
  
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
           x = "Survival time (years)",
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
  create_boxplot(dummy_results_RMST, "Causal effect (RMST, years)", "Dummy_Outcome_RMST_Boxplot.png")
  
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

sp_results <- dummy_results$SP_results
rmst_results <- dummy_results$RMST_results

print(summary(sp_results$CATE_Estimate))
print(summary(rmst_results$CATE_Estimate))

