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

simulate_simple_data <- function(n_samples = 2000, true_effect = 0.1, seed = 42) {
  set.seed(seed)
  
  # Create basic covariates
  data <- data.frame(
    Age = rnorm(n_samples, mean = 60, sd = 10),
    Stage = sample(1:4, n_samples, replace = TRUE),
    HPV = rbinom(n_samples, 1, 0.6)
  )
  
  # Assign treatment (high/low dose)
  data$high_dose <- rbinom(n_samples, 1, 0.5)
  
  baseline_risk <- exp(0.01 * (data$Age - 60) + 0.2 * data$Stage) * 0.2

  treatment_effect <- -log(1 - true_effect) * data$high_dose
  data$true_time <- rexp(n_samples, rate = baseline_risk * exp(treatment_effect))
  admin_censor <- rep(10, n_samples)
  
  random_censor <- rexp(n_samples, rate = 0.05) #i changed this down from 0.1
  
  data$time <- pmin(data$true_time, random_censor, admin_censor)
  data$event <- as.numeric(data$true_time <= pmin(random_censor, admin_censor))
  
  return(data)
}

validate_methodology <- function() {

  # Create simulated data with known 10% treatment effect
  sim_data <- simulate_simple_data(n_samples = 2000, true_effect = 0.10)
  
  time_points <- c(1, 2, 3, 4, 5, 6) 
  cat("\nCensoring pattern in simulated data:\n")
  for(t in time_points) {
    n_risk <- sum(sim_data$time >= t)
    pct <- (n_risk/nrow(sim_data)) * 100
    cat(sprintf("%d years: %.1f%% remaining (%d/%d)\n", 
                t, pct, n_risk, nrow(sim_data)))
  }
  
  horizons <- c(1, 2, 3, 4, 5, 6)
  results <- list()
  
  for(h in horizons) {
    cat(sprintf("\nTesting horizon: %d years\n", h))
    
    valid_horizon <- check_censoring(sim_data$time, sim_data$event, h)
    if(is.null(valid_horizon)) {
      cat("Invalid horizon due to censoring\n")
      next
    }
    
    X <- as.matrix(sim_data[, c("Age", "Stage", "HPV")])
    W <- sim_data$high_dose
    Y <- sim_data$time
    D <- sim_data$event
    
    weights <- calculate_ipcw(Y, D, X)
    
    cs_forest <- causal_survival_forest(
      X = X,
      W = W,
      Y = Y,
      D = D,
      horizon = valid_horizon,
      num.trees = 1000,
      sample.weights = weights
    )
    
    ate <- average_treatment_effect(cs_forest)
    
    results[[as.character(h)]] <- list(
      horizon = h,
      estimated_effect = ate[1],
      se = ate[2],
      n_at_risk = sum(Y >= h),
      percent_remaining = mean(Y >= h) * 100
    )
    
    cat("\nResults for horizon", h, "years:\n")
    cat(sprintf("  Estimated effect: %.1f%%\n", ate[1] * 100))
    cat(sprintf("  Standard error: %.1f%%\n", ate[2] * 100))
    cat(sprintf("  Patients at risk: %d (%.1f%%)\n", 
                sum(Y >= h), mean(Y >= h) * 100))
    cat(sprintf("  True effect: 10.0%%\n"))
    cat("\n")
  }
  
  return(results)
}


#real dataset/anaylsis

prepare_data <- function(df) {
  df %>%
    mutate(
      dose_intensity = d_Frac * Fx / Total_days_RT,
      high_dose = as.numeric(dose_intensity > median(dose_intensity, na.rm = TRUE)),
      event = as.numeric(Cause_of_Death_Status == 1),
      survival_time = Length_FU
    )
}

calculate_ipcw <- function(Y, D, X, batch_size = 100) {
  X_df <- as.data.frame(X)
  formula_str <- paste("Surv(Y, 1-D) ~", paste(colnames(X_df), collapse = " + "))
  simple_model <- coxph(as.formula(formula_str), data = data.frame(Y = Y, D = D, X_df))
  
  n <- length(Y)
  n_batches <- ceiling(n/batch_size)
  ipcw <- numeric(n)
  
  for(i in 1:n_batches) {
    start_idx <- (i-1)*batch_size + 1
    end_idx <- min(i*batch_size, n)
    batch_indices <- start_idx:end_idx
    
    surv_prob <- survfit(simple_model, 
                        newdata = X_df[batch_indices,],
                        y = TRUE)
    
    for(j in seq_along(batch_indices)) {
      idx <- batch_indices[j]
      if(D[idx] == 0) {
        ipcw[idx] <- 1
      } else {
        prob <- summary(surv_prob, times = Y[idx])$surv[1]
        ipcw[idx] <- 1/max(prob, 0.05)
      }
    }
  }
  
  ipcw <- ipcw/mean(ipcw)
  ipcw <- pmin(ipcw, quantile(ipcw, 0.99))
  ipcw <- pmax(ipcw, quantile(ipcw, 0.01))
  
  return(ipcw)
}

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
    labs(title = sprintf("Treatment Effects by Stage at %d Year Horizon", 
                        results$horizon),
         x = "Stage",
         y = "Effect on Survival Probability") +
    theme_bw()
  ggsave(sprintf("treatment_effects_by_stage_%dy.png", results$horizon), p)
}

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
    labs(title = sprintf("Variable Importance in Causal Forest at %d Year Horizon",
                        results$horizon),
         x = "Variable", 
         y = "Importance") +
    theme_bw()
  ggsave(sprintf("variable_importance_%dy.png", results$horizon), p)
}

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

run_subgroup_analysis <- function(results, df) {
  calc_subgroup_effect <- function(mask) {
    ate <- average_treatment_effect(results$forest, subset = mask)
    return(c(ate[1], ate[2]))
  }
  
  subgroup_results <- list()
  
  # HPV status analysis
  for (hpv in c(0, 1)) {
    mask <- df$HPV_Positive == hpv
    effect <- calc_subgroup_effect(mask)
    subgroup_results[[paste0("HPV_", hpv)]] <- list(
      effect = effect[1],
      ci_lower = effect[1] - 1.96*effect[2],
      ci_upper = effect[1] + 1.96*effect[2],
      n = sum(mask)
    )
  }
  
  # Stage analysis
  for (stage in sort(unique(df$Stage_numeric))) {
    mask <- df$Stage_numeric == stage
    if (sum(mask) >= 10) {
      effect <- calc_subgroup_effect(mask)
      subgroup_results[[paste0("Stage_", stage)]] <- list(
        effect = effect[1],
        ci_lower = effect[1] - 1.96*effect[2],
        ci_upper = effect[1] + 1.96*effect[2],
        n = sum(mask)
      )
    }
  }
  
  return(subgroup_results)
}

run_causal_survival_analysis <- function(df, horizon = 24) {
  X <- df %>%
    select(Age, Sex, Smoking_PY, Stage_numeric, 
           HPV_Positive, HPV_Unknown, Chemo, RT_year) %>%
    as.matrix()
  
  numeric_cols <- c("Age", "Smoking_PY", "RT_year")
  X[, numeric_cols] <- scale(X[, numeric_cols])
  
  W <- as.numeric(df$high_dose)
  Y <- df$survival_time
  D <- df$event
  
  valid_horizon <- check_censoring(Y, D, horizon)
  if (is.null(valid_horizon)) {
    stop(sprintf("Horizon %d years does not meet censoring requirements", horizon))
  }
  
  weights <- calculate_ipcw(Y, D, X)
  
  cs_forest <- causal_survival_forest(
    X = X,
    W = W,
    Y = Y,
    D = D,
    horizon = valid_horizon,
    num.trees = 2000,
    honesty = TRUE,
    sample.weights = weights,
    min.node.size = 10,
    seed = 42
  )
  
  predictions <- predict(cs_forest)
  ate <- average_treatment_effect(cs_forest)
  
  effects_df <- data.frame(effects = predictions$predictions)
  p1 <- ggplot(effects_df, aes(x = effects)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    geom_vline(xintercept = ate[1], color = "red", linetype = "dashed") +
    labs(title = sprintf("Distribution of Treatment Effects at %d Year Horizon", valid_horizon),
         x = "Effect on Survival Probability",
         y = "Count") +
    theme_bw()
  ggsave(sprintf("treatment_effects_distribution_%dy.png", valid_horizon), p1)
  
  return(list(
    forest = cs_forest,
    predictions = predictions,
    ate = ate,
    horizon = valid_horizon,
    weights = weights
  ))
}

analyze_followup <- function(df) {

  stats <- summary(df$Length_FU)
  print("Follow-up Length Summary (years):")
  print(stats)
  
  p <- ggplot(df, aes(x = Length_FU)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    geom_vline(xintercept = median(df$Length_FU, na.rm = TRUE), 
               color = "red", linetype = "dashed") +
    labs(title = "Distribution of Follow-up Length",
         x = "Follow-up Length (years)",
         y = "Count") +
    theme_bw()
  
  # censoring information at different time points
  time_points <- c(1, 2, 3, 6, 12, 24)
  cat("\nCensoring proportions at different time points:\n")
  for(t in time_points) {
    n_risk <- sum(df$Length_FU >= t)
    n_total <- nrow(df)
    cat(sprintf("%d years: %.1f%% remaining (%d/%d)\n", 
                t, 100*n_risk/n_total, n_risk, n_total))
  }
  n_events <- sum(df$event)
  cat(sprintf("\nTotal events: %d (%.1f%%)\n", 
              n_events, 100*n_events/nrow(df)))
  
  return(p)
}

main_analysis <- function() {
  df <- read.csv("selected_data_with_sites.csv")
  df <- prepare_data(df)

  cat("\n----------------------------------------\n")
  cat("Analyzing Follow-up Length Distribution\n")
  cat("----------------------------------------\n")
  p <- analyze_followup(df)
  ggsave("followup_distribution.png", p)
  
  suggested_horizons <- c(1, 2, 3)  # or any horizons you want to test (in years)
  valid_horizons <- numeric(0)
  
  Y <- df$survival_time
  D <- df$event
  
  for (h in suggested_horizons) {
    if (!is.null(check_censoring(Y, D, h))) {
      valid_horizons <- c(valid_horizons, h)
    }
  }
  
  if (length(valid_horizons) == 0) {
    stop("No valid horizons found with sufficient censoring probability")
  }
  
  all_results <- list()
  for (horizon in valid_horizons) {
    cat("\n----------------------------------------\n")
    cat(sprintf("Analysis for horizon: %d years\n", horizon))
    cat("----------------------------------------\n")
    
    results <- run_causal_survival_analysis(df, horizon)
    
    cat(sprintf("Number of observations: %d\n", nrow(df)))
    cat(sprintf("Number of events: %d\n", sum(df$event)))
    cat(sprintf("Censoring rate: %.2f%%\n", (1 - mean(df$event)) * 100))
    
    cat("\nAverage Treatment Effect:\n")
    cat(sprintf("Estimate: %.3f\n", results$ate[1]))
    cat(sprintf("Standard Error: %.3f\n", results$ate[2]))
    cat(sprintf("95%% CI: (%.3f, %.3f)\n", 
                results$ate[1] - 1.96*results$ate[2],
                results$ate[1] + 1.96*results$ate[2]))
    
    plot_effects_by_stage(results, df)
    plot_variable_importance(results)
    
    subgroup_results <- run_subgroup_analysis(results, df)
    results$subgroup_analysis <- subgroup_results
    
    cat("\nSubgroup Analysis Results:\n")
    for (group_name in names(subgroup_results)) {
      cat(sprintf("\n%s:\n", group_name))
      cat(sprintf("  Effect: %.3f\n", subgroup_results[[group_name]]$effect))
      cat(sprintf("  95%% CI: (%.3f, %.3f)\n", 
                  subgroup_results[[group_name]]$ci_lower,
                  subgroup_results[[group_name]]$ci_upper))
      cat(sprintf("  N: %d\n", subgroup_results[[group_name]]$n))
    }

    all_results[[as.character(results$horizon)]] <- results
  }

  plot_km_curves(df)
  plot_hazard_ratios(df)
  return(all_results)
}

cat("\nRunning methodology validation with simulated data...\n")
validation_results <- validate_methodology()

cat("\nRunning main analysis with real data...\n")
results <- main_analysis()
