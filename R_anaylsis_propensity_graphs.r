
if (!require("grf")) install.packages("grf")
if (!require("survival")) install.packages("survival")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("survminer")) install.packages("survminer")
if (!require("glmnet")) install.packages("glmnet")

library(grf)
library(survival)
library(dplyr)
library(ggplot2)
library(survminer)
library(glmnet)

# Updated simulation function with stronger true effect
simulate_simple_data <- function(n_samples = 2000, true_effect = 0.1, seed = 42) {
  set.seed(seed)
  
  # Generate covariates with more realistic distributions
  data <- data.frame(
    Age = rnorm(n_samples, mean = 60, sd = 8),
    Stage = sample(1:4, n_samples, replace = TRUE, 
                  prob = c(0.35, 0.35, 0.2, 0.1)),  # More early stage cases
    HPV = rbinom(n_samples, 1, 0.6)
  )
  
  # Even better balanced treatment assignment
  prob_high_dose <- plogis(0 +  # Perfectly balanced at baseline
                          0.15 * scale(data$Age) +
                          0.1 * data$Stage +
                          0.2 * data$HPV)
  
  data$high_dose <- rbinom(n_samples, 1, prob_high_dose)
  
  # Modified baseline risk for higher event rate
  baseline_risk <- exp(-2.5 +  # Increased from -3.0
                      0.1 * scale(data$Age) +
                      0.5 * (data$Stage/2))  # Stronger stage effect
  
  # Strengthened treatment effect
  # Using a stronger modification of the hazard
  hazard_multiplier <- exp(2 * log(1-true_effect) * data$high_dose)  # Doubled the effect
  
  # Generate survival times
  data$true_time <- rexp(n_samples, rate = baseline_risk * hazard_multiplier)
  
  # More aggressive censoring
  admin_censor <- rep(3, n_samples)  # Reduced from 4 to 3 years
  random_censor <- rexp(n_samples, rate = 0.03)  # Reduced from 0.05 for fewer early censorings
  
  data$time <- pmin(data$true_time, random_censor, admin_censor)
  data$event <- as.numeric(data$true_time <= pmin(random_censor, admin_censor))
  
  # Print detailed diagnostics
  event_rate <- mean(data$event)
  treatment_rate <- mean(data$high_dose)
  cat(sprintf("Event rate: %.1f%%\n", event_rate * 100))
  cat(sprintf("Treatment rate: %.1f%%\n", treatment_rate * 100))
  
  # Calculate true survival probabilities at multiple timepoints
  for(t in c(1,2,3)) {
    true_surv_control <- mean(data$true_time[data$high_dose == 0] > t)
    true_surv_treated <- mean(data$true_time[data$high_dose == 1] > t)
    cat(sprintf("True %d-year survival - Control: %.1f%%, Treated: %.1f%% (Diff: %.1f%%)\n", 
                t, true_surv_control * 100, true_surv_treated * 100,
                (true_surv_treated - true_surv_control) * 100))
  }
  
  # Check balance in covariates
  cat("\nCovariate balance:\n")
  for(var in c("Age", "Stage", "HPV")) {
    mean_control <- mean(data[[var]][data$high_dose == 0])
    mean_treated <- mean(data[[var]][data$high_dose == 1])
    cat(sprintf("%s - Control: %.2f, Treated: %.2f\n", 
                var, mean_control, mean_treated))
  }
  
  return(data)
}

# Propensity score functions
calculate_propensity_scores <- function(X, W, alpha = 0.5) {
  # Ensure X is a matrix and standardize numeric columns
  X <- as.matrix(X)
  numeric_cols <- which(apply(X, 2, function(x) length(unique(x)) > 10))
  if(length(numeric_cols) > 0) {
    X[,numeric_cols] <- scale(X[,numeric_cols])
  }
  
  # Convert to sparse matrix for glmnet
  x_matrix <- Matrix::Matrix(X, sparse = TRUE)
  
  # Fit elastic net with cross-validation
  cv_fit <- cv.glmnet(x_matrix, W, 
                      family = "binomial", 
                      alpha = alpha,
                      nfolds = 5)
  
  # Get propensity scores
  prop_scores <- predict(cv_fit, newx = x_matrix, 
                        s = "lambda.min", 
                        type = "response")
  
  # Print model coefficients for debugging
  coef_matrix <- as.matrix(coef(cv_fit, s = "lambda.min"))
  nonzero_coefs <- coef_matrix[coef_matrix != 0, , drop = FALSE]
  cat("\nNon-zero coefficients in propensity score model:\n")
  print(nonzero_coefs)
  
  return(as.vector(prop_scores))
}

trim_extreme_propensity_scores <- function(prop_scores, lower_quantile = 0.025, 
                                         upper_quantile = 0.975) {
  lower_bound <- quantile(prop_scores, lower_quantile)
  upper_bound <- quantile(prop_scores, upper_quantile)
  
  # Print trimming info
  cat(sprintf("\nTrimming bounds: [%.3f, %.3f]\n", lower_bound, upper_bound))
  
  keep_indices <- prop_scores >= lower_bound & prop_scores <= upper_bound
  cat(sprintf("Trimming removed %d observations (%.1f%%)\n", 
              sum(!keep_indices), mean(!keep_indices) * 100))
  
  return(keep_indices)
}

# Analysis functions
check_censoring <- function(Y, D, horizon) {
  km_cens <- survfit(Surv(Y, 1-D) ~ 1)
  times <- sort(unique(Y[Y <= horizon]))
  probs <- summary(km_cens, times = times)$surv
  
  n_total <- length(Y)
  valid_horizons <- times[which(
    probs > 0.05 &
      sapply(times, function(t) sum(Y >= t)) > 0.1 * n_total
  )]
  
  if (length(valid_horizons) == 0) return(NULL)
  
  valid_horizon <- max(valid_horizons[valid_horizons <= horizon])
  if (is.null(valid_horizon) || is.infinite(valid_horizon)) {
    valid_horizon <- max(valid_horizons)
  }
  
  return(valid_horizon)
}

summarize_propensity <- function(scores, label = "") {
  cat(sprintf("\n%s Propensity Score Summary:\n", label))
  cat("  Quantiles:\n")
  quants <- quantile(scores, probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
  for(q in names(quants)) {
    cat(sprintf("    %s: %.3f\n", q, quants[q]))
  }
  cat(sprintf("  Range: [%.3f, %.3f]\n", min(scores), max(scores)))
  cat(sprintf("  Mean: %.3f\n", mean(scores)))
  cat(sprintf("  SD: %.3f\n", sd(scores)))
}

run_analysis <- function(X, W, Y, D, horizon, method = "propensity") {
  valid_horizon <- check_censoring(Y, D, horizon)
  if (is.null(valid_horizon)) {
    stop(sprintf("Horizon %d years does not meet censoring requirements", horizon))
  }
  
  keep_indices <- rep(TRUE, length(W))
  
  if (method == "propensity") {
    prop_scores <- calculate_propensity_scores(X, W)
    summarize_propensity(prop_scores, "Before Trimming")
    
    keep_indices <- trim_extreme_propensity_scores(prop_scores)
    summarize_propensity(prop_scores[keep_indices], "After Trimming")
    
    X_trimmed <- X[keep_indices, ]
    W_trimmed <- W[keep_indices]
    Y_trimmed <- Y[keep_indices]
    D_trimmed <- D[keep_indices]
    
    cs_forest <- causal_survival_forest(
      X = X_trimmed,
      W = W_trimmed,
      Y = Y_trimmed,
      D = D_trimmed,
      horizon = valid_horizon,
      num.trees = 6000,
      honesty = TRUE,
      min.node.size = 5,
      alpha = 0.05,
      imbalance.penalty = 0.1,
      stabilize.splits = TRUE,
      seed = 42
    )
    
  } else {
    cs_forest <- causal_survival_forest(
      X = X,
      W = W,
      Y = Y,
      D = D,
      horizon = valid_horizon,
      num.trees = 6000,
      honesty = TRUE,
      min.node.size = 5,
      alpha = 0.05,
      imbalance.penalty = 0.1,
      stabilize.splits = TRUE,
      seed = 42
    )
  }
  
  ate <- average_treatment_effect(cs_forest)
  predictions <- predict(cs_forest)
  
  return(list(
    forest = cs_forest,
    predictions = predictions,
    ate = ate,
    horizon = valid_horizon,
    n_excluded = sum(!keep_indices),
    keep_indices = keep_indices
  ))
}

# Plotting functions
plot_effects_distribution <- function(results, method, horizon) {
  effects_df <- data.frame(effects = results$predictions$predictions)
  
  p <- ggplot(effects_df, aes(x = effects)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    geom_vline(xintercept = results$ate[1], color = "red", linetype = "dashed") +
    labs(title = sprintf("Distribution of Treatment Effects at %d Year Horizon (%s)", 
                        horizon, method),
         x = "Effect on Survival Probability",
         y = "Count") +
    theme_bw()
  
  ggsave(sprintf("treatment_effects_%s_%dy.png", method, horizon), p, 
         width = 7, height = 7)
}

plot_propensity_distributions <- function(scores, W, keep_indices = NULL, title = "Propensity Score Distributions by Treatment Group") {

  if (!is.null(keep_indices)) {
    scores <- scores[keep_indices]
    W <- W[keep_indices]
  }
  
  W <- as.numeric(W)
  
  scores_df <- data.frame(
    propensity = scores,
    treatment = factor(W, levels = c(0, 1), labels = c("Low Dose", "High Dose"))
  )
  
  p <- ggplot(scores_df, aes(x = propensity, fill = treatment)) +
    geom_density(alpha = 0.5) +
    labs(title = title,
         x = "Propensity Score",
         y = "Density",
         fill = "Treatment Group") +
    theme_bw()
  
  ggsave(paste0("propensity_dist_", gsub(" ", "_", tolower(title)), ".png"), p, 
         width = 8, height = 6)
  
  return(p)
}

plot_survival_curves <- function(Y, D, W, title = "Kaplan-Meier Survival Curves by Treatment Group") {
  surv_data <- data.frame(
    time = Y,
    status = D,
    treatment = factor(W, labels = c("Low Dose", "High Dose"))
  )
  
  fit <- survfit(Surv(time, status) ~ treatment, data = surv_data)
  
  p <- ggsurvplot(fit,
             data = surv_data,
             title = title,
             xlab = "Time (years)",
             ylab = "Survival Probability",
             conf.int = TRUE,
             risk.table = TRUE,
             pval = TRUE)
  
  # Save the plot
  ggsave(paste0("survival_curves_", gsub(" ", "_", tolower(title)), ".png"), 
         p$plot, width = 10, height = 8)
}

plot_covariate_balance <- function(X, W, varnames = NULL) {
  if(is.null(varnames)) varnames <- colnames(X)
  
  balance_stats <- data.frame(
    variable = varnames,
    diff = sapply(1:ncol(X), function(i) 
      mean(X[W==1,i]) - mean(X[W==0,i])),
    se = sapply(1:ncol(X), function(i) 
      sqrt(var(X[W==1,i])/sum(W==1) + var(X[W==0,i])/sum(W==0)))
  )
  
  p <- ggplot(balance_stats, aes(x = variable, y = diff)) +
    geom_point() +
    geom_errorbar(aes(ymin = diff - 1.96*se, ymax = diff + 1.96*se), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(title = "Standardized Mean Differences",
         x = "Covariate",
         y = "Standardized Mean Difference") +
    theme_bw()
  
  ggsave("covariate_balance.png", p, width = 8, height = 6)
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
    labs(title = paste(title, "-", var_name),
         x = var_name,
         y = "Estimated Treatment Effect") +
    theme_bw()
  
  ggsave(paste0("effects_by_", tolower(var_name), ".png"), p, width = 8, height = 6)
  
  return(p)
}

# Main analysis
main <- function() {
  # Part 1: Simulation Study
  cat("\nRunning simulation study...\n")
  sim_data <- simulate_simple_data(n_samples = 2000, true_effect = 0.1)
  
  X_sim <- as.matrix(sim_data[, c("Age", "Stage", "HPV")])
  W_sim <- sim_data$high_dose
  Y_sim <- sim_data$time
  D_sim <- sim_data$event
  
  horizons <- c(1, 2, 3)
  sim_results <- list()
  
  # Initial simulation plots
  plot_survival_curves(Y_sim, D_sim, W_sim, "Simulation Data - Overall Survival")
  plot_covariate_balance(X_sim, W_sim, colnames(X_sim))
  
  for (h in horizons) {
    cat(sprintf("\nAnalyzing horizon: %d years\n", h))
    
    # Run both methods
    prop_results <- run_analysis(X_sim, W_sim, Y_sim, D_sim, h, "propensity")
    orig_results <- run_analysis(X_sim, W_sim, Y_sim, D_sim, h, "original")
    
    sim_results[[as.character(h)]] <- list(
      propensity = prop_results,
      original = orig_results,
      true_effect = 0.1
    )
    
    # Print results
    cat("\nPropensity Score Method:\n")
    cat(sprintf("  Estimated effect: %.1f%%\n", prop_results$ate[1] * 100))
    cat(sprintf("  Standard error: %.1f%%\n", prop_results$ate[2] * 100))
    cat(sprintf("  Samples excluded: %d (%.1f%%)\n", 
                prop_results$n_excluded,
                prop_results$n_excluded/nrow(sim_data) * 100))
    
    cat("\nOriginal Method:\n")
    cat(sprintf("  Estimated effect: %.1f%%\n", orig_results$ate[1] * 100))
    cat(sprintf("  Standard error: %.1f%%\n", orig_results$ate[2] * 100))
    cat(sprintf("True effect: 10.0%%\n"))
    
    # Generate plots for simulation data
    plot_effects_distribution(prop_results, "propensity", h)
    plot_effects_distribution(orig_results, "original", h)
    
    plot_propensity_distributions(
      prop_results$forest$W.hat, 
      W_sim, 
      keep_indices = prop_results$keep_indices,
      paste("Simulation Data - Horizon", h)
    )
    
    # Plot treatment effects by covariates
    for(var in c("Age", "Stage", "HPV")) {
      plot_treatment_effects_by_covariate(
        prop_results$predictions, 
        X_sim, 
        var,
        keep_indices = prop_results$keep_indices,
        paste("Simulation - Horizon", h)
      )
    }
  }
  
  # Part 2: Real Data Analysis
  cat("\nAnalyzing real data...\n")
  real_data <- read.csv("selected_data_with_sites.csv")
  
  # Prepare real data
  real_data <- real_data %>%
    mutate(
      dose_intensity = d_Frac * Fx / Total_days_RT,
      high_dose = as.numeric(dose_intensity > median(dose_intensity, na.rm = TRUE)),
      event = as.numeric(Cause_of_Death_Status == 1),
      survival_time = Length_FU
    )
  
  X_real <- real_data %>%
    select(Age, Sex, Smoking_PY, Stage_numeric, 
           HPV_Positive, HPV_Unknown, Chemo, RT_year) %>%
    as.matrix()
  
  # Scale numeric columns
  numeric_cols <- c("Age", "Smoking_PY", "RT_year")
  X_real[, numeric_cols] <- scale(X_real[, numeric_cols])
  
  W_real <- real_data$high_dose
  Y_real <- real_data$survival_time
  D_real <- real_data$event
  
  # Initial real data plots
  plot_survival_curves(Y_real, D_real, W_real, "Real Data - Overall Survival")
  plot_covariate_balance(X_real, W_real, colnames(X_real))
  
  real_results <- list()
  
  for (h in c(1, 2, 3)) {
    cat(sprintf("\nAnalyzing horizon: %d years\n", h))
    
    # Run both methods
    prop_results <- run_analysis(X_real, W_real, Y_real, D_real, h, "propensity")
    orig_results <- run_analysis(X_real, W_real, Y_real, D_real, h, "original")
    
    real_results[[as.character(h)]] <- list(
      propensity = prop_results,
      original = orig_results
    )
    
    cat("\nPropensity Score Method:\n")
    cat(sprintf("  Estimated effect: %.1f%%\n", prop_results$ate[1] * 100))
    cat(sprintf("  Standard error: %.1f%%\n", prop_results$ate[2] * 100))
    cat(sprintf("  Samples excluded: %d (%.1f%%)\n", 
                prop_results$n_excluded,
                prop_results$n_excluded/nrow(real_data) * 100))
    
    cat("\nOriginal Method:\n")
    cat(sprintf("  Estimated effect: %.1f%%\n", orig_results$ate[1] * 100))
    cat(sprintf("  Standard error: %.1f%%\n", orig_results$ate[2] * 100))
    
    # Generate plots for real data
    plot_effects_distribution(prop_results, "propensity_real", h)
    plot_effects_distribution(orig_results, "original_real", h)
    plot_propensity_distributions(
      prop_results$forest$W.hat, 
      W_real, 
      keep_indices = prop_results$keep_indices,
      paste("Real Data - Horizon", h)
    )
    
    for(var in c("Age", "Stage_numeric", "HPV_Positive")) {
      plot_treatment_effects_by_covariate(
        prop_results$predictions, 
        X_real, 
        var,
        keep_indices = prop_results$keep_indices,
        paste("Real Data - Horizon", h)
      )
    }
  }
  
  return(list(
    simulation_results = sim_results,
    real_data_results = real_results
  ))
}

results <- main()
