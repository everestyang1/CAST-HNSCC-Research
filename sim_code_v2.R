
# Libraries
rm(list=ls()) #will remove ALL objects 

library(caret)
library(car)
library(boot)
library(pscl)
library(tidyverse)
library(Hmisc)
library(psych)
library(missForest)
library(pdp)
library(glmnet)
library(grf)
library(survival)
library(devEMF)
library(matrixStats)
library(vioplot)
library(fastshap)
library(MASS)
library(corrplot)
library(GGally)
library(dplyr)


# # simulate data

# n_samples <- 5000
# set.seed(199)

# data <- data.frame(
#   Age = rnorm(n_samples, mean = 60, sd = 10),
#   Stage = sample(1:4, n_samples, replace = TRUE),
#   Chemotherapy = rbinom(n_samples, 1, 0.5),
#   HPV = rbinom(n_samples, 1, 0.6),
#   Smoking = rbinom(n_samples, 1, 0.25),
#   TumorSite_OralCavity = as.numeric(sample(1:4, n_samples, replace = TRUE) == 1),
#   TumorSite_Oropharynx = as.numeric(sample(1:4, n_samples, replace = TRUE) == 2),
#   TumorSite_Larynx = as.numeric(sample(1:4, n_samples, replace = TRUE) == 3),
#   TumorSite_Hypopharynx = as.numeric(sample(1:4, n_samples, replace = TRUE) == 4),
#   Comorbidities_0 = as.numeric(sample(0:2, n_samples, replace = TRUE, prob = c(0.5, 0.3, 0.2)) == 0),
#   Comorbidities_1 = as.numeric(sample(0:2, n_samples, replace = TRUE, prob = c(0.5, 0.3, 0.2)) == 1),
#   Comorbidities_2plus = as.numeric(sample(0:2, n_samples, replace = TRUE, prob = c(0.5, 0.3, 0.2)) == 2),
#   NeckDissection = rbinom(n_samples, 1, 0.4),
#   FeedingTube = rbinom(n_samples, 1, 0.2),
#   Tracheostomy = rbinom(n_samples, 1, 0.1),
#   AlcoholUse_None = as.numeric(sample(1:3, n_samples, replace = TRUE) == 1),
#   AlcoholUse_Moderate = as.numeric(sample(1:3, n_samples, replace = TRUE) == 2),
#   AlcoholUse_Heavy = as.numeric(sample(1:3, n_samples, replace = TRUE) == 3),
#   Education_HighSchool = as.numeric(sample(1:3, n_samples, replace = TRUE) == 1),
#   Education_College = as.numeric(sample(1:3, n_samples, replace = TRUE) == 2),
#   Education_Graduate = as.numeric(sample(1:3, n_samples, replace = TRUE) == 3),
#   MaritalStatus_Single = as.numeric(sample(1:3, n_samples, replace = TRUE) == 1),
#   MaritalStatus_Married = as.numeric(sample(1:3, n_samples, replace = TRUE) == 2),
#   MaritalStatus_Divorced = as.numeric(sample(1:3, n_samples, replace = TRUE) == 3),
#   BMI = rnorm(n_samples, mean = 25, sd = 4),
#   PerformanceStatus_0 = as.numeric(sample(0:4, n_samples, replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.1, 0.1)) == 0),
#   PerformanceStatus_1 = as.numeric(sample(0:4, n_samples, replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.1, 0.1)) == 1),
#   PerformanceStatus_2 = as.numeric(sample(0:4, n_samples, replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.1, 0.1)) == 2),
#   PerformanceStatus_3 = as.numeric(sample(0:4, n_samples, replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.1, 0.1)) == 3),
#   PerformanceStatus_4 = as.numeric(sample(0:4, n_samples, replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.1, 0.1)) == 4)
# ) 

# str(data)

# # causal (treatment) feature
# Cause <- rbinom(n_samples, 1, 0.5)

# # survival time with true effect

# true_effect <- 12

# Survival_time <- 6+round(pmax(rnorm(n_samples, mean = 10 + true_effect*Cause, sd = 7), 1), 0)

# # Survival status with censoring
# #Survival_status <- rbinom(n_samples, 1, 0.5)

# censoring_prob <- function(time) {
#   pmin(0.01 + 0.2 * (time / max(time)), 0.95)  
# }

# # Generate survival status with increasing censoring probability
# Survival_status <- rbinom(n_samples, 1, 1 - censoring_prob(Survival_time))

# cor(Survival_time, Survival_status)


# # compile dataset
# cancer_data_encoded <- cbind(data, Cause, Survival_time, Survival_status)
# str(cancer_data_encoded)

# ######

# # check censoring

# check_censoring <- function(Y, D, horizon) {
#   km_cens <- survfit(Surv(Y, 1-D) ~ 1)
#   times <- sort(unique(Y[Y <= horizon]))
#   probs <- summary(km_cens, times = times)$surv
  
#   n_total <- length(Y)
#   valid_horizons <- times[which(
#     probs > 0.05 &
#       sapply(times, function(t) sum(Y >= t)) > 0.1 * n_total
#   )]
  
#   if (length(valid_horizons) == 0) return(NULL)
  
#   valid_horizon <- max(valid_horizons[valid_horizons <= horizon])
#   if (is.null(valid_horizon) || is.infinite(valid_horizon)) {
#     valid_horizon <- max(valid_horizons)
#   }
  
#   # Calculate fraction of non-censored samples at the valid horizon
#   fraction_non_censored <- sum(Y >= valid_horizon & D == 1) / n_total
  
#   return(list(valid_horizon = valid_horizon, fraction_non_censored = fraction_non_censored))
# }

# # set horizons
# horizons <- c(12*1, 12*2, 12*3, 12*4, 12*5)

# for(h in horizons) {
#   cat(sprintf("\nTesting horizon: %d months", h))
  
#   result <- check_censoring(cancer_data_encoded$Survival_time, cancer_data_encoded$Survival_status, h)
#   if(is.null(result)) {
#     cat("\nInvalid horizon due to censoring\n")
#     next
#   }
  
#   cat(sprintf("\nValid horizon: %d months", result$valid_horizon))
#   cat(sprintf("\nFraction of non-censored samples: %.2f%%\n", result$fraction_non_censored * 100))
# }



# #######
# # summarize continuous variables
# library(reshape2)
# continuous_vars <- sapply(cancer_data_encoded, function(x) length(unique(x)) > 2)
# continuous_data <- cancer_data_encoded[, continuous_vars]

# summary_data <- summary(continuous_data)
# summary_df <- as.data.frame(summary_data)

# # Write the dataframe to a CSV file
# write.csv(summary_df, file = "summary_continuous_data.csv", row.names = FALSE)

# ##########

cancer_data_encoded <- read.csv("cancer_data_chemo_cause.csv")
# Prepare data for ML analysis
# Create training and test sets
set.seed(99)


index <- createDataPartition(cancer_data_encoded$Survival_status, p=0.75, list=FALSE)
trainSet <- cancer_data_encoded[index,]
testSet <- cancer_data_encoded[-index,]
nrow(trainSet)
nrow(testSet)
table(trainSet$Survival_status)
table(testSet$Survival_status)

# Create the necessary matrices for the causal forest
trainSet_X <- as.data.frame(subset(trainSet, select=-c(Survival_time, Survival_status, Cause)))
trainSet_W <- trainSet$Cause
trainSet_times <- trainSet$Survival_time
trainSet_events <- trainSet$Survival_status

testSet_X <- as.data.frame(subset(testSet, select=-c(Survival_time, Survival_status, Cause)))
testSet_W <- testSet$Cause
testSet_times <- testSet$Survival_time
testSet_events <- testSet$Survival_status

# Fit logistic regression with elastic net penalty to get propensity scores for Cause
X_matrix_train <- model.matrix(~ . - Cause - Survival_time - Survival_status, data = trainSet)[,-1]
Y_vector_train <- trainSet$Cause
X_matrix_test <- model.matrix(~ . - Cause - Survival_time - Survival_status, data = testSet)[,-1]

set.seed(1234)
alpha_values <- seq(0.01, 0.99, by = 0.01)
cv_models <- list()

for (alpha_val in alpha_values) {
  cv_model <- cv.glmnet(X_matrix_train, Y_vector_train,
                        family = "binomial",
                        type.measure = "class",
                        alpha = alpha_val,
                        nfolds = 10)
  cv_models[[as.character(alpha_val)]] <- cv_model
}

optimal_lambdas <- sapply(cv_models, function(model) model$lambda.min)
best_alpha <- alpha_values[which.min(optimal_lambdas)]
optimal_lambda <- cv_models[[as.character(best_alpha)]]$lambda.min
final_model <- glmnet(X_matrix_train, Y_vector_train,
                      family = "binomial",
                      alpha = best_alpha,
                      lambda = optimal_lambda)
cv_glmnet_model <- final_model

best_alpha
optimal_lambda

cv_glmnet_model_coefs <- as.data.frame(as.matrix(coef(cv_glmnet_model)))
write.csv(cv_glmnet_model_coefs, file = "cv_glmnet_model_coefs.csv", row.names = TRUE)


# Predict propensity scores for training and testing sets
W.hat_train <- predict(cv_glmnet_model, newx = X_matrix_train, type = "response", s = optimal_lambda)
W.hat_test <- predict(cv_glmnet_model, newx = X_matrix_test, type = "response", s = optimal_lambda)

# Examine propensity scores
summary(W.hat_train)
hist(W.hat_train)

summary(W.hat_test)
hist(W.hat_test)

# Identify patients with propensity scores very close to 0 or 1
trainSet_with_scores <- as.data.frame(cbind(trainSet, W.hat_train))
colnames(trainSet_with_scores) <- c(colnames(trainSet), "W_hat")
head(trainSet_with_scores)

testSet_with_scores <- as.data.frame(cbind(testSet, W.hat_test))
colnames(testSet_with_scores) <- c(colnames(testSet), "W_hat")
head(testSet_with_scores)

write.csv(trainSet_with_scores, "trainSet_with_scores.csv", row.names = FALSE)
write.csv(testSet_with_scores, "testSet_with_scores.csv", row.names = FALSE)

# Plot Pearson correlation matrix
varscor_trainSet_with_scores_pearson <- corr.test(trainSet_with_scores, method = "pearson", adjust = "bonf", alpha = .05, ci = FALSE)
varscor_trainSet_with_scores_pearson_p <- varscor_trainSet_with_scores_pearson$p

write.csv(varscor_trainSet_with_scores_pearson_p, file = "varscor_trainSet_with_scores_pearson_p.csv", row.names = TRUE)
write.csv(varscor_trainSet_with_scores_pearson$r, file = "varscor_trainSet_with_scores_pearson_r.csv", row.names = TRUE)

emf("correlation_matrix_pearson_propensity_scores_train.emf", width = 7, height = 7, bg = "transparent",
    fg = "black", pointsize = 8, family = "Arial", coordDPI = 600)
corrplot(varscor_trainSet_with_scores_pearson$r, p.mat = varscor_trainSet_with_scores_pearson$p, method = 'circle', 
         tl.col = "black", type = "upper", sig.level = 0.05, pch.cex = 0.6, cl.cex = 1, tl.cex = 1, insig = 'pch', 
         pch = 19, pch.col = "white", diag = FALSE, font = 1)
dev.off()

# Filter training and test sets by propensity score range
trainSet_with_scores_filtered <- subset(trainSet_with_scores, W_hat >= 0.1 & W_hat <= 0.9)
testSet_with_scores_filtered <- subset(testSet_with_scores, W_hat >= 0.1 & W_hat <= 0.9)

trainSet_X <- as.data.frame(subset(trainSet_with_scores_filtered, select = -c(W_hat, Survival_time, Survival_status, Cause)))
trainSet_W <- trainSet_with_scores_filtered$Cause
trainSet_times <- trainSet_with_scores_filtered$Survival_time
trainSet_events <- trainSet_with_scores_filtered$Survival_status
W_hat_train_adj2 <- trainSet_with_scores_filtered$W_hat

testSet_X <- as.data.frame(subset(testSet_with_scores_filtered, select = -c(W_hat, Survival_time, Survival_status, Cause)))
testSet_W <- testSet_with_scores_filtered$Cause
testSet_times <- testSet_with_scores_filtered$Survival_time
testSet_events <- testSet_with_scores_filtered$Survival_status
W_hat_test_adj2 <- testSet_with_scores_filtered$W_hat


write.csv(trainSet_with_scores_filtered, "trainSet_with_scores_filtered.csv", row.names = FALSE)
write.csv(testSet_with_scores_filtered, "testSet_with_scores_filtered.csv", row.names = FALSE)


nrow(trainSet_with_scores)
nrow(trainSet_with_scores_filtered)


nrow(testSet_with_scores)
nrow(testSet_with_scores_filtered)


# visualize KM plots
library(survival)
library(survminer)

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
                   ggtheme = theme_minimal(),
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
ggsave(filename = "Kaplan_Meier_Survival_Curves.png", plot = p_KM$plot, dpi = 300, width = 10, height = 8)

# Save the risk table
ggsave(filename = "Kaplan_Meier_Survival_Curves_risk_table.png", plot = p_KM$table, dpi = 300, width = 10, height = 4)


# Implement causal forests
n_trees_val <- 5000
horizons <- c(12*1, 12*2, 12*3, 12*4, 12*5)
results_SP <- data.frame(horizon_sel = integer(), ATE_estimate_train_SP = numeric(), ATE_se_train_SP = numeric(),
                      ATE_estimate_test_SP = numeric(), ATE_se_test_SP = numeric())

for (horizon in horizons) {
  csf_model_SP <- causal_survival_forest(X = trainSet_X, Y = trainSet_times, W = trainSet_W,
                                      D = trainSet_events, W.hat = as.vector(W_hat_train_adj2), 
                                      num.trees = n_trees_val, target = "survival.probability", 
                                      horizon = horizon, seed = 123)
  
  ate_train_SP <- average_treatment_effect(csf_model_SP)
  csf_pred_test_SP <- predict(csf_model_SP, testSet_X, W.hat = as.vector(W_hat_test_adj2), estimate.variance = TRUE)
  ate_h_SP_test <- mean(csf_pred_test_SP$predictions)
  ate_h_SP_test_sd <- mean(sqrt(csf_pred_test_SP$variance.estimates))
  
  results_SP <- rbind(results_SP, data.frame(horizon_sel = horizon, 
                                       ATE_estimate_train_SP = ate_train_SP[1], ATE_se_train_SP = ate_train_SP[2],
                                       ATE_estimate_test_SP = ate_h_SP_test, ATE_se_test_SP = ate_h_SP_test_sd))
}

print(results_SP)
write.csv(results_SP, "train_and_test_ATE_SP.csv", row.names = FALSE)


########

# same for RMST
results_RMST <- data.frame(horizon_sel = integer(), ATE_estimate_train_RMST = numeric(), ATE_se_train_RMST = numeric(),
                         ATE_estimate_test_RMST = numeric(), ATE_se_test_RMST = numeric())

for (horizon in horizons) {
  csf_model_RMST <- causal_survival_forest(X = trainSet_X, Y = trainSet_times, W = trainSet_W,
                                         D = trainSet_events, W.hat = as.vector(W_hat_train_adj2), 
                                         num.trees = n_trees_val, target = "RMST", 
                                         horizon = horizon, seed = 123)
  
  ate_train_RMST <- average_treatment_effect(csf_model_RMST)
  csf_pred_test_RMST <- predict(csf_model_RMST, testSet_X, W.hat = as.vector(W_hat_test_adj2), estimate.variance = TRUE)
  ate_h_RMST_test <- mean(csf_pred_test_RMST$predictions)
  ate_h_RMST_test_sd <- mean(sqrt(csf_pred_test_RMST$variance.estimates))
  
  results_RMST <- rbind(results_RMST, data.frame(horizon_sel = horizon, 
                                             ATE_estimate_train_RMST = ate_train_RMST[1], ATE_se_train_RMST = ate_train_RMST[2],
                                             ATE_estimate_test_RMST = ate_h_RMST_test, ATE_se_test_RMST = ate_h_RMST_test_sd))
}

print(results_RMST)
write.csv(results_RMST, "train_and_test_ATE_RMST.csv", row.names = FALSE)

#######



##############
# Fit a causal survival forest at one specific horizon 
horizon_sel <- 2 * 12
csf_model_try1 <- causal_survival_forest(X = trainSet_X, Y = trainSet_times, W = trainSet_W,
                                         D = trainSet_events, W.hat = as.vector(W_hat_train_adj2), num.trees = n_trees_val, 
                                         target = "survival.probability", horizon = horizon_sel, seed = 123)

# Estimate ATE
average_treatment_effect(csf_model_try1)

###

###{r}

# Best linear projection over features of interest
# best_linear_projection(csf_model_try1, subset(trainSet_X, select = c("Grade", "Age", "Sex")))
blp_results <- best_linear_projection(csf_model_try1, trainSet_X)

# Extract the coefficients and their statistics
blp_matrix <- as.matrix(blp_results)

# Ensure the matrix has the correct dimensions and column names
colnames(blp_matrix) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

# Create a data frame from the matrix
blp_df <- data.frame(
  Term = rownames(blp_matrix),
  Estimate = blp_matrix[, 1],
  `Std. Error` = blp_matrix[, 2],
  `t value` = blp_matrix[, 3],
  `Pr(>|t|)` = blp_matrix[, 4]
)

# Export the data frame to a CSV file
write.csv(blp_df, "best_linear_projection.csv", row.names = FALSE)

#############


# Fake confounder sensitivity test
# Create a fake confounder variable which affects both the treatment and outcome with some 
# strength K_val. Check sensitivity of the causal effects to this variable as function of K_val

# Define K_val values to iterate over
K_vals <- seq(-2.4, 2.4, by = 0.2)

# Store ATE results
ate_results <- data.frame(K_val = numeric(), ATE = numeric())

# Original ATE for comparison
original_ate <- average_treatment_effect(csf_model_try1)

# Iterate over K_val values
for (K_val in K_vals) {
  
  set.seed(123)  # For reproducibility
  
  # Generate fake confounder
  FakeConfounder <- exp(rnorm(nrow(trainSet_X)))
  
  # Modify times based on K_val
  modified_times <- trainSet_times * exp(-(K_val / 100) * FakeConfounder)
  
  # Modify treatment
  original_binary_vector <- trainSet_W
  
  # Define the probability of flipping
  flip_probability <- min(0.5, max(0, abs(K_val) / 100))
  
  # Flip values based on the defined probability
  flipped_vector <- ifelse(runif(length(original_binary_vector)) < flip_probability,
                           1 - original_binary_vector,  # Flip 0 to 1 and 1 to 0
                           original_binary_vector)      # Keep original value
  
  # Fit causal survival forest model
  csf_model_Fake <- causal_survival_forest(X = trainSet_X, Y = modified_times,
                                           W = flipped_vector,
                                           D = trainSet_events,
                                           W.hat = as.vector(W_hat_train_adj2),
                                           num.trees = n_trees_val,
                                           target = "survival.probability", 
                                           horizon = horizon_sel,
                                           seed = 123)
  
  # Estimate ATE for the current K_val
  ate_value <- average_treatment_effect(csf_model_Fake)
  
  # Store results in the data frame
  ate_results <- rbind(ate_results, data.frame(K_val = K_val, ATE = ate_value))
}

# Print results
print(ate_results)

# Compare with original ATE
print(paste("Original ATE:", original_ate))

# Create a new DataFrame for structured results
structured_results <- data.frame(K_val = numeric(),
                                 ATE_estimate = numeric(),
                                 ATE_std_err = numeric())

# Extract estimates and standard errors
for (i in seq_along(K_vals)) {
  K_val <- K_vals[i]
  
  # Extract the estimate and standard error for the current K_val
  estimate <- ate_results[i * 2 - 1, 2]  # Odd rows contain estimates
  std_err <- ate_results[i * 2, 2]       # Even rows contain standard errors
  
  # Append to the new DataFrame
  structured_results <- rbind(structured_results,
                              data.frame(K_val = K_val,
                                         ATE_estimate = estimate,
                                         ATE_std_err = std_err))
}

# Save the structured results as a CSV file
write.csv(structured_results, "ate_results.csv", row.names = FALSE)

# Print the structured results to verify
print(structured_results)

# Create the plot
ate_plot <- ggplot(structured_results, aes(x = K_val, y = ATE_estimate)) +
  geom_point(size = 6, color = "blue") +  # Points for ATE estimates
  geom_errorbar(aes(ymin = ATE_estimate - ATE_std_err, 
                    ymax = ATE_estimate + ATE_std_err), 
                width = 0.2, color = "black") +  # Error bars
  labs(title = "Sensitivity to fake confounder",
       x = "Strength of fake confounder effect",
       y = "ATE for survival probability") +
  theme_minimal() +  # Clean theme
  theme(text = element_text(size = 14))  # Increase text size for better readability

# Save the plot as a high-resolution PNG file
ggsave("ate_plot.png", plot = ate_plot, width = 10, height = 6, dpi = 500)

# Print the plot to the console (optional)
print(ate_plot)

######


########

# Dummy outcome test on causal survival forest: causal effect should be zero at all times
num_repetitions <- 30
set.seed(1234)
dummy_results <- data.frame(
  Repetition = integer(),
  Horizon = integer(),
  CATE_Estimate = numeric(),
  Standard_Error = numeric(),
  cor_values = numeric()
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
                                           target = "survival.probability", 
                                           horizon = horizon, 
                                           seed=123)
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

write.csv(dummy_results, "Dummy_Outcome_Normal_Dist_Results.csv", row.names = FALSE)

# Dummy Outcome Boxplot
dummy_results$Horizon <- as.factor(dummy_results$Horizon)

emf("Dummy_Outcome_Boxplot.emf", width = 7, height = 7,
    bg = "transparent", fg = "black", pointsize = 12,
    family = "Arial", coordDPI = 600)

ggplot(dummy_results, aes(x = Horizon, y = CATE_Estimate, fill = Horizon)) +
  geom_boxplot() +
  labs(title = "Boxplot of CATE Estimates with Dummy Outcome",
       x = "Time after treatment (years)",
       y = "Causal effect (SP)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "none")

dev.off()

###

# Now dummy for RMST

set.seed(1234)
dummy_RMST_results <- data.frame(
  Repetition = integer(),
  Horizon = integer(),
  CATE_Estimate = numeric(),
  Standard_Error = numeric(),
  cor_values = numeric()
)

for (rep in 1:num_repetitions) {
  trainSet_W_Dummy_RMST <- sample(trainSet_W, length(trainSet_W), replace=TRUE)
  train_Dummy_RMSTOutcomeTimes <- sample(trainSet_times, length(trainSet_times), replace=TRUE)
  
  for (horizon in horizons) {
    forest_dummy_RMST <- causal_survival_forest(X = trainSet_X, 
                                                Y = train_Dummy_RMSTOutcomeTimes, 
                                                W = trainSet_W_Dummy_RMST, 
                                                W.hat = as.vector(W_hat_train_adj2),
                                                D = trainSet_events, 
                                                num.trees = n_trees_val, 
                                                target = "RMST", 
                                                horizon = horizon, 
                                                seed=123)
    ate_dummy_RMST <- average_treatment_effect(forest_dummy_RMST)
    
    if (is.list(ate_dummy_RMST)) {
      cate_estimate <- ate_dummy_RMST$estimate 
      standard_error <- ate_dummy_RMST$std.err
    } else if (is.vector(ate_dummy_RMST)) {
      cate_estimate <- ate_dummy_RMST[1] 
      standard_error <- ate_dummy_RMST[2]
    }
    
    dummy_RMST_results <- rbind(dummy_RMST_results, data.frame(
      Repetition = rep,
      Horizon = horizon,
      CATE_Estimate = cate_estimate,
      Standard_Error = standard_error
    ))
  }
}

write.csv(dummy_RMST_results, "Dummy_RMST_Outcome_Normal_Dist_Results.csv", row.names = FALSE)

# Dummy_RMST Outcome Boxplot
dummy_RMST_results$Horizon <- as.factor(dummy_RMST_results$Horizon)

emf("Dummy_RMST_Outcome_Boxplot.emf", width = 7, height = 7,
    bg = "transparent", fg = "black", pointsize = 12,
    family = "Arial", coordDPI = 600)

ggplot(dummy_RMST_results, aes(x = Horizon, y = CATE_Estimate, fill = Horizon)) +
  geom_boxplot() +
  labs(title = "Boxplot of CATE Estimates with Dummy Outcome",
       x = "Time after treatment (years)",
       y = "Causal effect (RMST, years)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "none")

dev.off()


###{r}

# Cross-validation to assess stability of causal forest on real data
numFolds <- trainControl(method = "cv", number = 10)
folds <- createFolds(trainSet_times, k = 10)
HN_cate_cv_results_RMST_train <- data.frame(Horizon = integer(), Estimate = numeric(), Standard_Error = numeric())

for(i in 1:10){
  testIndexes <- folds[[i]]
  trainSet_X_cv <- trainSet_X[-testIndexes,]
  trainSet_times_cv <- trainSet_times[-testIndexes]
  trainSet_W_cv <- trainSet_W[-testIndexes]
  trainSet_events_cv <- trainSet_events[-testIndexes]
  W_hat_train_adj2_cv <- W_hat_train_adj2[-testIndexes]
  
  for (h in horizons) {
    forest_h_RMST_train <- causal_survival_forest(X = trainSet_X_cv, Y = trainSet_times_cv, W = trainSet_W_cv, W.hat=as.vector(W_hat_train_adj2_cv),
                                                  D = trainSet_events_cv, num.trees = n_trees_val, target = "RMST",
                                                  horizon = h, seed=123)
    ate_h_RMST_train <- average_treatment_effect(forest_h_RMST_train)
    HN_cate_cv_results_RMST_train <- rbind(HN_cate_cv_results_RMST_train, data.frame(Horizon = h, Estimate = ate_h_RMST_train[1], Standard_Error = ate_h_RMST_train[2]))
  }
}

print(HN_cate_cv_results_RMST_train)
write.csv(HN_cate_cv_results_RMST_train, file = "HN_cate_cv_results_RMST_train.csv", row.names = TRUE)

# Create boxplot for cross-validation results
HN_cate_cv_results_RMST_train$Horizon <- as.factor(HN_cate_cv_results_RMST_train$Horizon)

emf("CV_Outcome_Boxplot.emf", width = 7, height = 7, bg = "transparent", fg = "black", pointsize = 12, family = "Arial", coordDPI = 600)

ggplot(HN_cate_cv_results_RMST_train, aes(x = Horizon, y = Estimate, fill = Horizon)) +
  geom_boxplot() +
  labs(title = "Boxplot of RMST CATE Estimates with CV", x = "Time after treatment (years)", y = "Causal effect (RMST, years)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "none")

dev.off()

# Another way to plot
HN_cate_cv_results_RMST_train_pl <- HN_cate_cv_results_RMST_train


# Convert Horizon to numeric
HN_cate_cv_results_RMST_train_pl$Horizon_num <- as.numeric(as.character(HN_cate_cv_results_RMST_train_pl$Horizon))

# Define the width of the jitter
jitter_width <- 3

# Add jitter to Horizon and save it as a new variable
HN_cate_cv_results_RMST_train_pl$Horizon_jittered <- jitter(HN_cate_cv_results_RMST_train_pl$Horizon_num, amount=jitter_width)

emf("CV_Outcome_errorplot.emf", width = 7, height = 7, bg = "transparent", fg = "black", pointsize = 12, family = "Arial", coordDPI = 600)

ggplot(HN_cate_cv_results_RMST_train_pl, aes(x = Horizon_jittered, y = Estimate)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = Estimate - Standard_Error, ymax = Estimate + Standard_Error), 
                width = 0.2) +
  #  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, color = "red") +
  labs(title = "Plot of RMST CATE Estimates with CV", x = "Time after treatment (years)", y = "Causal effect (RMST, years)") +
  theme_minimal() +
  theme(legend.position = "none") + xlim(0, 170)

dev.off()



# Predict causal forest outputs on testing data for RMST
test_results_RMST <- data.frame(Horizon = integer(), Estimate = numeric(), Standard_Error = numeric())

for (h in horizons) {
  forest_h_RMST_test <- causal_survival_forest(X = trainSet_X, Y = trainSet_times, W = trainSet_W, W.hat = as.vector(W_hat_train_adj2),
                                               D = trainSet_events, num.trees = n_trees_val, target = "RMST", horizon = h, seed = 123)
  c.pred_RMST_test <- predict(forest_h_RMST_test, testSet_X, W.hat = as.vector(W_hat_test_adj2), estimate.variance = TRUE)
  ate_h_RMST_test <- mean(c.pred_RMST_test$predictions)
  ate_h_RMST_test_sd <- mean(sqrt(c.pred_RMST_test$variance.estimates))
  
  test_results_RMST <- rbind(test_results_RMST, data.frame(Horizon = h, Estimate = ate_h_RMST_test, Standard_Error = ate_h_RMST_test_sd))
}

print(test_results_RMST)
write.csv(test_results_RMST, file = "test_results_RMST_test.csv", row.names = TRUE)

# Same for survival probability
test_results_survival_probability <- data.frame(Horizon = integer(), Estimate = numeric(), Standard_Error = numeric())

for (h in horizons) {
  forest_h_survival_probability_test <- causal_survival_forest(X = trainSet_X, Y = trainSet_times, W = trainSet_W, W.hat = as.vector(W_hat_train_adj2),
                                                               D = trainSet_events, num.trees = n_trees_val, target = "survival.probability", horizon = h, seed = 123)
  c.pred_survival_probability_test <- predict(forest_h_survival_probability_test, testSet_X, W.hat = as.vector(W_hat_test_adj2), estimate.variance = TRUE)
  ate_h_survival_probability_test <- mean(c.pred_survival_probability_test$predictions)
  ate_h_survival_probability_test_sd <- mean(sqrt(c.pred_survival_probability_test$variance.estimates))
  
  test_results_survival_probability <- rbind(test_results_survival_probability, data.frame(Horizon = h, Estimate = ate_h_survival_probability_test, Standard_Error = ate_h_survival_probability_test_sd))
}

print(test_results_survival_probability)
write.csv(test_results_survival_probability, file = "test_results_survival_probability_test.csv", row.names = TRUE)

###

# same for survival probability


# Cross-validation to assess stability on real data
numFolds <- trainControl(method = "cv", number = 10)
folds <- createFolds(trainSet_times, k = 10)
HN_cate_cv_results_SP_train <- data.frame(Horizon = integer(), Estimate = numeric(), Standard_Error = numeric())

for(i in 1:10){
  testIndexes <- folds[[i]]
  trainSet_X_cv <- trainSet_X[-testIndexes,]
  trainSet_times_cv <- trainSet_times[-testIndexes]
  trainSet_W_cv <- trainSet_W[-testIndexes]
  trainSet_events_cv <- trainSet_events[-testIndexes]
  W_hat_train_adj2_cv <- W_hat_train_adj2[-testIndexes]
  
  for (h in horizons) {
    forest_h_SP_train <- causal_survival_forest(X = trainSet_X_cv, Y = trainSet_times_cv, W = trainSet_W_cv, W.hat=as.vector(W_hat_train_adj2_cv),
                                                D = trainSet_events_cv, num.trees = n_trees_val, target = "survival.probability",
                                                horizon = h, seed=123)
    ate_h_SP_train <- average_treatment_effect(forest_h_SP_train)
    HN_cate_cv_results_SP_train <- rbind(HN_cate_cv_results_SP_train, data.frame(Horizon = h, Estimate = ate_h_SP_train[1], Standard_Error = ate_h_SP_train[2]))
  }
}

print(HN_cate_cv_results_SP_train)
write.csv(HN_cate_cv_results_SP_train, file = "HN_cate_cv_results_SP_train.csv", row.names = TRUE)

# Create boxplot for cross-validation results
HN_cate_cv_results_SP_train$Horizon <- as.factor(HN_cate_cv_results_SP_train$Horizon)

emf("CV_Outcome_Boxplot_SP.emf", width = 7, height = 7, bg = "transparent", fg = "black", pointsize = 12, family = "Arial", coordDPI = 600)

ggplot(HN_cate_cv_results_SP_train, aes(x = Horizon, y = Estimate, fill = Horizon)) +
  geom_boxplot() +
  labs(title = "Boxplot of SP CATE Estimates with CV", x = "Time after treatment (years)", y = "Causal effect (SP)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "none")

dev.off()

# Another way to plot
HN_cate_cv_results_SP_train_pl <- HN_cate_cv_results_SP_train


# Convert Horizon to numeric
HN_cate_cv_results_SP_train_pl$Horizon_num <- as.numeric(as.character(HN_cate_cv_results_SP_train_pl$Horizon))

# Define the width of the jitter
jitter_width <- 3

# Add jitter to Horizon and save it as a new variable
HN_cate_cv_results_SP_train_pl$Horizon_jittered <- jitter(HN_cate_cv_results_SP_train_pl$Horizon_num, amount=jitter_width)

emf("CV_Outcome_errorplot_SP.emf", width = 7, height = 7, bg = "transparent", fg = "black", pointsize = 12, family = "Arial", coordDPI = 600)

ggplot(HN_cate_cv_results_SP_train_pl, aes(x = Horizon_jittered, y = Estimate)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = Estimate - Standard_Error, ymax = Estimate + Standard_Error), 
                width = 0.2) +
  #  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, color = "red") +
  labs(title = "Plot of SP CATE Estimates with CV", x = "Time after treatment (years)", y = "Causal effect (SP)") +
  theme_minimal() +
  theme(legend.position = "none") + xlim(0, 170)

dev.off()

# Predict on testing data for SP
test_results_SP <- data.frame(Horizon = integer(), Estimate = numeric(), Standard_Error = numeric())

for (h in horizons) {
  forest_h_SP_test <- causal_survival_forest(X = trainSet_X, Y = trainSet_times, W = trainSet_W, W.hat = as.vector(W_hat_train_adj2),
                                             D = trainSet_events, num.trees = n_trees_val, target = "survival.probability", horizon = h, seed = 123)
  c.pred_SP_test <- predict(forest_h_SP_test, testSet_X, W.hat = as.vector(W_hat_test_adj2), estimate.variance = TRUE)
  ate_h_SP_test <- mean(c.pred_SP_test$predictions)
  ate_h_SP_test_sd <- mean(sqrt(c.pred_SP_test$variance.estimates))
  
  test_results_SP <- rbind(test_results_SP, data.frame(Horizon = h, Estimate = ate_h_SP_test, Standard_Error = ate_h_SP_test_sd))
}

print(test_results_SP)
write.csv(test_results_SP, file = "test_results_SP_test.csv", row.names = TRUE)

# Same for survival probability
test_results_survival_probability <- data.frame(Horizon = integer(), Estimate = numeric(), Standard_Error = numeric())

for (h in horizons) {
  forest_h_survival_probability_test <- causal_survival_forest(X = trainSet_X, Y = trainSet_times, W = trainSet_W, W.hat = as.vector(W_hat_train_adj2),
                                                               D = trainSet_events, num.trees = n_trees_val, target = "survival.probability", horizon = h, seed = 123)
  c.pred_survival_probability_test <- predict(forest_h_survival_probability_test, testSet_X, W.hat = as.vector(W_hat_test_adj2), estimate.variance = TRUE)
  ate_h_survival_probability_test <- mean(c.pred_survival_probability_test$predictions)
  ate_h_survival_probability_test_sd <- mean(sqrt(c.pred_survival_probability_test$variance.estimates))
  
  test_results_survival_probability <- rbind(test_results_survival_probability, data.frame(Horizon = h, Estimate = ate_h_survival_probability_test, Standard_Error = ate_h_survival_probability_test_sd))
}

print(test_results_survival_probability)
write.csv(test_results_survival_probability, file = "test_results_survival_probability_test.csv", row.names = TRUE)

###


###{r}

# Negative control with an irrelevant causal variable: causal effect should be zero at all times
set.seed(123)
binary_vector <- sample(c(0, 1), size = nrow(trainSet_X), replace = TRUE)

trainSet_W_neg_c <- binary_vector

trainSet_X_neg_c <- trainSet_X

HN_cate_results_RMST_train_neg_c <- data.frame(Horizon = integer(), Estimate = numeric(), Standard_Error = numeric())

for (h in horizons) {
  forest_h_RMST_train_neg_c <- causal_survival_forest(X = trainSet_X_neg_c, Y = trainSet_times, W = trainSet_W_neg_c, W.hat = as.vector(W_hat_train_adj2),
                                                      D = trainSet_events, num.trees = n_trees_val, target = "RMST", horizon = h, tune.parameters = "all")
  ate_h_RMST_train_neg_c <- average_treatment_effect(forest_h_RMST_train_neg_c)
  HN_cate_results_RMST_train_neg_c <- rbind(HN_cate_results_RMST_train_neg_c, data.frame(Horizon = h, Estimate = ate_h_RMST_train_neg_c[1], Standard_Error = ate_h_RMST_train_neg_c[2]))
}

print(HN_cate_results_RMST_train_neg_c)
write.csv(HN_cate_results_RMST_train_neg_c, file = "HN_cate_results_RMST_train_neg_c.csv", row.names = TRUE)

# Same for survival probabilities with negative control
HN_cate_results_RMST_train_neg_c_surv <- data.frame(Horizon = integer(), Estimate = numeric(), Standard_Error = numeric())

for (h in horizons) {
  forest_h_RMST_train_neg_c_surv <- causal_survival_forest(X = trainSet_X_neg_c, Y = trainSet_times, W = trainSet_W_neg_c, W.hat = as.vector(W_hat_train_adj2),
                                                           D = trainSet_events, num.trees = n_trees_val, target = "survival.probability", horizon = h, tune.parameters = "all")
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

# noise should not matter much

#####





###########


# SHAP values for causal forest at selected horizon
horizon_sel <- 5*12
forest_sel <- causal_survival_forest(X = trainSet_X, Y = trainSet_times, W = trainSet_W, W.hat = as.vector(W_hat_train_adj2), D = trainSet_events, num.trees = n_trees_val, target = "survival.probability", horizon = horizon_sel, tune.parameters = "all", seed=1234)
forest_sel_preds_train <- predict(forest_sel, trainSet_X, estimate.variance = TRUE)
forest_sel_preds_test <- predict(forest_sel, testSet_X, estimate.variance = TRUE)

forest_sel_preds_train_df <- cbind(trainSet_X, forest_sel_preds_train)
forest_sel_preds_test_df <- cbind(testSet_X, forest_sel_preds_test)

write.csv(forest_sel_preds_train_df, file = "HN_forest_sel_preds_train_df.csv", row.names = FALSE)
write.csv(forest_sel_preds_test_df, file = "HN_forest_sel_preds_test_df.csv", row.names = FALSE)



# Compute SHAP values
pfun <- function(object, newdata) {
  predict(object, newdata = newdata, estimate.variance = TRUE)$predictions
}

# these are approximate SHAP values calculated using a Monte Carlo method

nsim_shap <- 100

shap <- fastshap::explain(forest_sel, X = trainSet_X, pred_wrapper = pfun, nsim = nsim_shap)
colnames(shap) <- paste0(colnames(shap), "_SHAP")
shap_vals <- cbind(trainSet_X, shap, forest_sel_preds_train)

# Calculate the average prediction for all samples
average_prediction <- mean(shap_vals$predictions)

# Normalize the SHAP values so the Monte Carlo random error is removed as much as possible
# Step 1: Sum all the columns ending in _SHAP
shap_cols <- grep("_SHAP$", names(shap_vals), value = TRUE)
shap_vals$sum_SHAP <- rowSums(shap_vals[, shap_cols])

# Step 2: Normalize each _SHAP column
for (col in shap_cols) {
  norm_col <- paste0(col, "_norm_SHAP")
  shap_vals[[norm_col]] <- (shap_vals[[col]] / shap_vals$sum_SHAP) * (shap_vals$predictions - average_prediction)
}

# Step 3: Create a sum of the normalized SHAP columns
norm_shap_cols <- grep("_norm_SHAP$", names(shap_vals), value = TRUE)
shap_vals$sum_norm_SHAP <- rowSums(shap_vals[, norm_shap_cols])

# Check if sum of normalized SHAP columns equals the difference between predictions and average prediction
all.equal(shap_vals$sum_norm_SHAP, shap_vals$predictions - average_prediction)

write.csv(shap_vals, file = "HN_shap_vals_raw.csv", row.names = FALSE)

# Identify columns to remove
cols_to_remove <- grep("_SHAP$", names(shap_vals), value = TRUE)

# Identify columns to keep (those ending with _SHAP_norm_SHAP)
cols_to_keep <- grep("_SHAP_norm_SHAP$", names(shap_vals), value = TRUE)

# Remove the _SHAP_norm_SHAP columns from the list of columns to remove
cols_to_remove <- setdiff(cols_to_remove, cols_to_keep)

# Remove the identified columns
shap_vals <- shap_vals[, !(names(shap_vals) %in% cols_to_remove)]

# Identify columns ending with _SHAP_norm_SHAP
cols_to_rename <- grep("_SHAP_norm_SHAP$", names(shap_vals), value = TRUE)

# Create new names
new_names <- sub("_SHAP_norm_SHAP$", "_SHAP", cols_to_rename)

# Rename the columns
names(shap_vals)[names(shap_vals) %in% cols_to_rename] <- new_names

str(shap_vals)

# Identify all columns ending with _SHAP
shap_cols <- grep("_SHAP$", names(shap_vals), value = TRUE)

# Calculate the sum of all _SHAP columns
shap_sum <- rowSums(shap_vals[, shap_cols])

# Compare with the predictions column: theoretically they should correspond
is_equal <- all.equal(shap_sum, shap_vals$predictions, tolerance = 1e-6)

if (isTRUE(is_equal)) {
  print("The sum of all _SHAP columns equals the predictions column.")
} else {
  print("The sum of all _SHAP columns does NOT equal the predictions column.")
  
  # Calculate and print the maximum difference
  max_diff <- max(abs(shap_sum - shap_vals$predictions))
  print(paste("Maximum difference:", max_diff))
  
  # Optionally, you can print summary statistics of the differences
  diff_summary <- summary(shap_sum - shap_vals$predictions)
  print("Summary of differences:")
  print(diff_summary)
}




write.csv(shap_vals, file = "HN_shap_vals.csv", row.names = FALSE)

# plot age effects for interest

emf("Age_SHAP_plot.emf", width = 7, height = 7, bg = "transparent", fg = "black", pointsize = 12, family = "Arial", coordDPI = 600)

plot(shap_vals$Age, shap_vals$Age_SHAP, xlab="Age (years)", ylab="SHAP value")

dev.off()


# Visualizing SHAP values
ordered_shap_vals <- shap_vals[order(shap_vals$predictions), ]
plot_data <- data.frame(
  Patients = 1:nrow(ordered_shap_vals),
  predictions = ordered_shap_vals$predictions,
  lower = ordered_shap_vals$predictions - 1.96 * sqrt(ordered_shap_vals$variance.estimates),
  upper = ordered_shap_vals$predictions + 1.96 * sqrt(ordered_shap_vals$variance.estimates)
)

emf("causal_effects_plot.emf", width = 7, height = 7, bg = "transparent", fg = "black", pointsize = 12, family = "Arial", coordDPI = 600)
ggplot(plot_data, aes(x = Patients, y = predictions)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(x = "Patients", y = "Causal effect (SP)")

dev.off()

emf("causal_effects_histogram.emf", width = 7, height = 7, bg = "transparent", fg = "black", pointsize = 12, family = "Arial", coordDPI = 600)
hist(ordered_shap_vals$predictions, xlab = "Causal effect (SP)", main = "Causal effects histogram")

# Add a dashed vertical red line at 0
abline(v = 0, col = "red", lty = 2)

# Add a title to the histogram
title(main = "Causal effects histogram")
box()

dev.off()

# make another version with fraction of effects <0 vs >=0
# Calculate the fractions
fraction_less_than_0 <- sum(ordered_shap_vals$predictions < 0) / length(ordered_shap_vals$predictions)
fraction_greater_equal_0 <- sum(ordered_shap_vals$predictions >= 0) / length(ordered_shap_vals$predictions)

emf("causal_effects_histogram_frac.emf", width = 7, height = 7, bg = "transparent", fg = "black", pointsize = 12, family = "Arial", coordDPI = 600)
# Create the histogram
# Create the histogram
hist(ordered_shap_vals$predictions, xlab = "Causal effect (SP)", main = "Causal effects histogram")

# Add a dashed vertical red line at 0
abline(v = 0, col = "red", lty = 2)

# Add the fractions to the plot in the upper right corner
text(x = max(ordered_shap_vals$predictions) * 0.7, y = max(hist(ordered_shap_vals$predictions, plot = FALSE)$counts) * 0.9, 
     labels = paste("Fraction < 0: ", round(fraction_less_than_0, 2)), col = "blue", pos = 4)
text(x = max(ordered_shap_vals$predictions) * 0.7, y = max(hist(ordered_shap_vals$predictions, plot = FALSE)$counts) * 0.8, 
     labels = paste("Fraction >= 0: ", round(fraction_greater_equal_0, 2)), col = "blue", pos = 4)

# Add a box around the plot
box()
dev.off()


# Select columns ending with _SHAP
shap_cols <- grep("_SHAP$", names(shap_vals), value = TRUE)

# Calculate median values for these _SHAP columns
medians <- apply(shap_vals[shap_cols], 2, median)

# Calculate median of absolute values for these _SHAP columns
abs_medians <- apply(abs(shap_vals[shap_cols]), 2, median)

# Combine results into a data frame
shap_vals_medians <- data.frame(medians, abs_medians)
write.csv(shap_vals_medians, file = "HN_shap_vals_medians.csv", row.names = TRUE)

# Sort shap_vals_medians by abs_medians in descending order
shap_vals_medians <- shap_vals_medians[order(-shap_vals_medians$abs_medians), ]

# Save the top 10
shap_vals_medians_top_10 <- head(shap_vals_medians, 10)
write.csv(shap_vals_medians_top_10, file = "HN_shap_vals_medians_top_10.csv", row.names = TRUE)

# Plot correlations between top SHAP values
std_devs <- apply(shap_vals, 2, sd)
cols <- names(std_devs)[std_devs > 0]
selected_shap_vals <- shap_vals[, cols]

varscor_selected_shap_vals <- corr.test(selected_shap_vals, method = "spearman", adjust = "bonf", alpha = .05, ci = FALSE)
varscor_selected_shap_vals_p <- varscor_selected_shap_vals$p

write.csv(varscor_selected_shap_vals_p, file = "varscor_selected_shap_vals_spearman_p.csv", row.names = TRUE)
write.csv(varscor_selected_shap_vals$r, file = "varscor_selected_shap_vals_spearman_r.csv", row.names = TRUE)

# Plot Spearman correlation matrix
emf("correlation_matrix_spearman_SHAP.emf", width = 7, height = 7, bg = "transparent", fg = "black", pointsize = 8, family = "Arial", coordDPI = 600)

corrplot(varscor_selected_shap_vals$r, p.mat = varscor_selected_shap_vals$p, method = 'circle', tl.col = "black", type = "upper", sig.level = 0.05, pch.cex = 0.6, cl.cex = 1, tl.cex = 1, insig = 'pch', pch = 19, pch.col = "white", diag = FALSE, font = 1)

dev.off()

# Plot Pearson correlation matrix
varscor_selected_shap_vals_pearson <- corr.test(selected_shap_vals, method = "pearson", adjust = "bonf", alpha = .05, ci = FALSE)
varscor_selected_shap_vals_pearson_p <- varscor_selected_shap_vals_pearson$p

write.csv(varscor_selected_shap_vals_pearson_p, file = "varscor_selected_shap_vals_pearson_p.csv", row.names = TRUE)
write.csv(varscor_selected_shap_vals_pearson$r, file = "varscor_selected_shap_vals_pearson_r.csv", row.names = TRUE)

emf("correlation_matrix_pearson_SHAP.emf", width = 7, height = 7, bg = "transparent", fg = "black", pointsize = 8, family = "Arial", coordDPI = 600)

corrplot(varscor_selected_shap_vals_pearson$r, p.mat = varscor_selected_shap_vals_pearson$p, method = 'circle', tl.col = "black", type = "upper", sig.level = 0.05, pch.cex = 0.6, cl.cex = 1, tl.cex = 1, insig = 'pch', pch = 19, pch.col = "white", diag = FALSE, font = 1)

dev.off()


# graph the SHAPs


# remove causal effects for now - only keep features and corresponding SHAP values
SHAP_data_1 <- as.data.frame(subset(shap_vals,
                                    select=-c(predictions, variance.estimates)))

# identify feature columns
features <- names(SHAP_data_1)[!grepl("_SHAP$", names(SHAP_data_1))]

# causal effects vector
causal_effects <- shap_vals$predictions

for (feature in features) {
  # Create the plot
  p <- ggplot(SHAP_data_1, aes(x = .data[[feature]], y = .data[[paste0(feature, "_SHAP")]])) +
    geom_point(aes(color = causal_effects, shape = factor(Sex)), size = 3) +
    scale_shape_manual(values = c(16, 15)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(causal_effects)) +
    labs(title = paste("SHAP values for", feature),
         x = feature,
         y = paste0(feature, "_SHAP"),
         color = "Causal Effects",
         shape = "Sex") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed
  
  # Save the plot as an EMF file
  emf_file <- paste0(make.names(feature), "_SHAP.emf")
  emf(file = emf_file)
  print(p)
  dev.off()
}


dev.off()

