if (!require("grf")) install.packages("grf")
if (!require("survival")) install.packages("survival")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("survminer")) install.packages("survminer")

library(grf)
library(survival)
library(dplyr)
library(ggplot2)
library(survminer)

prepare_data <- function(df) {
  df <- df %>%
    mutate(
      dose_intensity = d_Frac * Fx / Total_days_RT,
      high_dose = as.numeric(dose_intensity > median(dose_intensity)),
      event = as.numeric(Cause_of_Death_Status == 1),
      survival_time = Length_FU
    )
  return(df)
}

# Function to check censoring and determine valid horizon
check_censoring <- function(Y, D, horizon) {
  # Calculate censoring probability at horizon
  km_cens <- survfit(Surv(Y, 1-D) ~ 1)
  
  # Get probabilities at all times
  times <- sort(unique(Y[Y <= horizon]))
  probs <- summary(km_cens, times = times)$surv
  
  # Find valid horizons where:
  # 1. We have at least 5%/10% censoring probability (i tried both)
  # 2. We have at least 10% of original sample size
  n_total <- length(Y)
  valid_horizons <- times[which(
    probs > 0.05 & #i tried 0.05 and 0.1 but it was not working for either which means our data lacks for horizon = 3 months
    sapply(times, function(t) sum(Y >= t)) > 0.1 * n_total
  )]
  
  if (length(valid_horizons) == 0) {
    return(NULL)
  }
  
  # Return largest valid horizon up to requested horizon
  valid_horizon <- max(valid_horizons[valid_horizons <= horizon])
  if (is.null(valid_horizon) || is.infinite(valid_horizon)) {
    valid_horizon <- max(valid_horizons)
  }
  
  return(valid_horizon)
}

# Main analysis function
run_causal_survival_analysis <- function(df, horizon = 24) {
  # Prepare features
  X <- df %>%
    select(Age, Sex, Smoking_PY, Stage_numeric, 
           HPV_Positive, HPV_Unknown, Chemo, RT_year) %>%
    as.matrix()
  
  # Standardize numeric features
  numeric_cols <- c("Age", "Smoking_PY", "RT_year")
  X[, numeric_cols] <- scale(X[, numeric_cols])
  
  # Binary treatment variable
  W <- as.numeric(df$high_dose)
  
  # Survival outcome
  Y <- df$survival_time
  D <- df$event
  
  # Check censoring and get valid horizon
  valid_horizon <- check_censoring(Y, D, horizon)
  if (is.null(valid_horizon)) {
    stop(sprintf("Horizon %d months does not meet censoring requirements", horizon))
  }
  if (valid_horizon != horizon) {
    message(sprintf("Adjusting horizon from %d to %d months due to censoring constraints", 
                   horizon, valid_horizon))
    horizon <- valid_horizon
  }
  
  # Add stabilized weights for treatment
  # IPCW (Inverse Probability of Censoring Weighting) through stabilized weights:
  get_stabilized_weights <- function(W, X) {
    ps_model <- glm(W ~ X, family = binomial())
    ps <- predict(ps_model, type = "response")
    w_mean <- mean(W)
    sw <- ifelse(W == 1, w_mean/ps, (1-w_mean)/(1-ps))
    sw <- pmin(sw, quantile(sw, 0.99))
    sw <- pmax(sw, quantile(sw, 0.01))
    return(sw)
  }
  
  # Get stabilized weights
  sw <- get_stabilized_weights(W, X)
  
  # Fit the causal survival forest
  cs_forest <- causal_survival_forest(
    X = X,
    W = W,
    Y = Y,
    D = D,
    horizon = horizon,
    num.trees = 2000,
    honesty = TRUE,
    sample.weights = sw
  )
  
  # Get predictions and treatment effects
  predictions <- predict(cs_forest)
  ate <- average_treatment_effect(cs_forest)
  
  # Create visualization of treatment effects
  effects_df <- data.frame(effects = predictions$predictions)
  p1 <- ggplot(effects_df, aes(x = effects)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    geom_vline(xintercept = ate[1], color = "red", linetype = "dashed") +
    labs(title = sprintf("Distribution of Treatment Effects at %d Month Horizon", horizon),
         x = "Effect on Survival Probability",
         y = "Count") +
    theme_bw()
  ggsave(sprintf("treatment_effects_distribution_%dm.png", horizon), p1)
  
  return(list(
    forest = cs_forest,
    predictions = predictions,
    ate = ate,
    horizon = horizon,
    weights = sw
  ))
}

# Plot Kaplan-Meier curves
plot_km_curves <- function(df) {
  fit <- survfit(Surv(survival_time, event) ~ high_dose, data = df)
  km_plot <- ggsurvplot(fit, data = df, 
                        risk.table = TRUE,
                        pval = TRUE,
                        conf.int = TRUE,
                        title = "Kaplan-Meier Curves by Dose Intensity",
                        legend.labs = c("Low Dose", "High Dose"))
  ggsave("km_curves.png", km_plot$plot + theme_bw())
}

# Plot effects by stage
plot_effects_by_stage <- function(results, df) {
  stage_effects <- data.frame(
    Stage = df$Stage_numeric,
    Effect = results$predictions$predictions
  ) %>%
    group_by(Stage) %>%
    summarise(
      mean_effect = mean(Effect),
      se = sd(Effect) / sqrt(n())
    )
  
  p <- ggplot(stage_effects, aes(x = factor(Stage), y = mean_effect)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean_effect - 1.96*se, 
                      ymax = mean_effect + 1.96*se), width = 0.2) +
    labs(title = sprintf("Treatment Effects by Stage at %d Month Horizon", 
                        results$horizon),
         x = "Stage",
         y = "Effect on Survival Probability") +
    theme_bw()
  ggsave(sprintf("treatment_effects_by_stage_%dm.png", results$horizon), p)
}

# Plot variable importance
plot_variable_importance <- function(results) {
  var_imp <- variable_importance(results$forest)
  var_names <- c("Age", "Sex", "Smoking_PY", "Stage", 
                 "HPV+", "HPV Unknown", "Chemo", "RT Year")
  importance_df <- data.frame(
    Variable = var_names,
    Importance = var_imp
  )
  p <- ggplot(importance_df, aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = sprintf("Variable Importance in Causal Forest at %d Month Horizon",
                        results$horizon),
         x = "Variable", 
         y = "Importance") +
    theme_bw()
  ggsave(sprintf("variable_importance_%dm.png", results$horizon), p)
}

# Plot hazard ratios
plot_hazard_ratios <- function(df) {
  cox_model <- coxph(Surv(survival_time, event) ~ Age + HPV_Positive + 
                     Smoking_PY + Chemo + high_dose, data = df)
  cox_summary <- summary(cox_model)
  
  hr_df <- data.frame(
    Variable = rownames(cox_summary$coefficients),
    HR = exp(cox_summary$coefficients[, 1]),
    LowerCI = exp(cox_summary$conf.int[, 3]),
    UpperCI = exp(cox_summary$conf.int[, 4])
  )
  
  p <- ggplot(hr_df, aes(x = reorder(Variable, HR), y = HR)) +
    geom_point() +
    geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +
    coord_flip() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title = "Hazard Ratios from Cox Model",
         x = "Variable",
         y = "Hazard Ratio (95% CI)") +
    theme_bw()
  ggsave("hazard_ratios.png", p)
}

# Subgroup analysis function
run_subgroup_analysis <- function(results, df) {
  calc_subgroup_effect <- function(mask) {
    ate <- average_treatment_effect(results$forest, subset = mask)
    return(c(ate[1], ate[2]))
  }
  
  # HPV status analysis
  print(sprintf("\nHPV Status Analysis at %d month horizon:", results$horizon))
  for (hpv in c(0, 1)) {
    mask <- df$HPV_Positive == hpv
    effect <- calc_subgroup_effect(mask)
    cat(sprintf("HPV %d: Effect = %.3f (95%% CI: %.3f, %.3f), n=%d\n",
                hpv, effect[1], 
                effect[1] - 1.96*effect[2],
                effect[1] + 1.96*effect[2],
                sum(mask)))
  }
  
  # Stage analysis
  print(sprintf("\nStage Analysis at %d month horizon:", results$horizon))
  for (stage in sort(unique(df$Stage_numeric))) {
    mask <- df$Stage_numeric == stage
    if (sum(mask) >= 10) {
      effect <- calc_subgroup_effect(mask)
      cat(sprintf("Stage %d: Effect = %.3f (95%% CI: %.3f, %.3f), n=%d\n",
                  stage, effect[1],
                  effect[1] - 1.96*effect[2],
                  effect[1] + 1.96*effect[2],
                  sum(mask)))
    }
  }
}

# Main execution
main_analysis <- function() {
  # Read and prepare data
  df <- read.csv("selected_data_with_sites.csv")
  df <- prepare_data(df)
  
  # Basic data summary
  print("Data Summary:")
  print(summary(df[c("survival_time", "dose_intensity", "event", "high_dose")]))
  
  # Check valid horizons first
  suggested_horizons <- c(1, 2, 3)
  valid_horizons <- numeric(0)
  
  # Pre-check which horizons are valid
  Y <- df$survival_time
  D <- df$event
  
  for (h in suggested_horizons) {
    if (!is.null(check_censoring(Y, D, h))) {
      valid_horizons <- c(valid_horizons, h)
    } else {
      cat(sprintf("\nSkipping horizon %d months due to insufficient censoring\n", h))
    }
  }
  
  if (length(valid_horizons) == 0) {
    stop("No valid horizons found with sufficient censoring probability")
  }
  
  cat("\nAnalyzing at valid horizons:", valid_horizons, "months\n")
  
  # Run analysis for each valid horizon
  all_results <- list()
  for (horizon in valid_horizons) {
    tryCatch({
      cat(sprintf("\nAnalyzing horizon: %d months\n", horizon))
      results <- run_causal_survival_analysis(df, horizon)
      
      cat(sprintf("\nTreatment Effect at %d months: %.3f (95%% CI: %.3f, %.3f)\n",
                  results$horizon,
                  results$ate[1],
                  results$ate[1] - 1.96*results$ate[2],
                  results$ate[1] + 1.96*results$ate[2]))
      
      # Run subgroup analysis
      run_subgroup_analysis(results, df)
      
      # Generate visualizations
      plot_effects_by_stage(results, df)
      plot_variable_importance(results)
      
      all_results[[as.character(results$horizon)]] <- results
    }, error = function(e) {
      cat(sprintf("\nError analyzing horizon %d months: %s\n", horizon, e$message))
    })
  }
  
  # Generate overall visualizations
  plot_km_curves(df)
  plot_hazard_ratios(df)
  
  return(all_results)
}

# Run the analysis
results <- main_analysis()