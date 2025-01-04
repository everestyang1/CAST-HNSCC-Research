# Install and load required packages
if (!require("grf")) install.packages("grf")
if (!require("survival")) install.packages("survival")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")

library(grf)
library(survival)
library(dplyr)
library(ggplot2)

# Read data
df <- read.csv("selected_data_with_sites.csv")

# Data preparation
prepare_data <- function(df) {
  # Calculate dose intensity and create necessary variables
  df <- df %>%
    mutate(
      dose_intensity = d_Frac * Fx / Total_days_RT,
      dose_intensity_scaled = as.vector(scale(dose_intensity)),
      event = as.numeric(Cause_of_Death_Status == 1),
      survival_time = Length_FU  # Using Length of Follow-up as survival time
    )
  
  return(df)
}

# Main analysis function
run_causal_analysis <- function(df) {
  # Prepare features
  X <- df %>%
    select(Age, Sex, Smoking_PY, Stage_numeric, 
           HPV_Positive, HPV_Unknown, Chemo, RT_year) %>%
    as.matrix()
  
  # Standardize numeric features
  numeric_cols <- c("Age", "Smoking_PY", "RT_year")
  X[, numeric_cols] <- scale(X[, numeric_cols])
  
  # Treatment variable
  W <- as.vector(df$dose_intensity_scaled)
  
  # Outcome variable - use survival time weighted by event status
  Y <- as.vector(df$survival_time * (1 + df$event))  # Weight events more heavily
  
  # Fit the causal forest
  cf <- causal_forest(
    X = X,
    Y = Y,
    W = W,
    num.trees = 2000,
    honesty = TRUE
  )
  
  # Get treatment effects
  tau.hat <- predict(cf)
  
  # Calculate average treatment effect
  ate <- average_treatment_effect(cf)
  
  # Plot treatment effects
  p1 <- ggplot(data.frame(effects = tau.hat$predictions), aes(x = effects)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    geom_vline(xintercept = ate[1], color = "red", linetype = "dashed") +
    labs(title = "Distribution of Treatment Effects",
         x = "Effect on Survival Time",
         y = "Count") +
    theme_minimal()
  print(p1)
  
  # Plot by stage
  stage_effects <- data.frame(
    Stage = df$Stage_numeric,
    Effect = tau.hat$predictions
  ) %>%
    group_by(Stage) %>%
    summarise(
      mean_effect = mean(Effect),
      se = sd(Effect) / sqrt(n())
    )
  
  p2 <- ggplot(stage_effects, aes(x = factor(Stage), y = mean_effect)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean_effect - 1.96*se, 
                      ymax = mean_effect + 1.96*se), width = 0.2) +
    labs(title = "Treatment Effects by Stage",
         x = "Stage",
         y = "Effect on Survival Time") +
    theme_minimal()
  print(p2)
  
  # Return results
  return(list(
    forest = cf,
    tau.hat = tau.hat,
    ate = ate
  ))
}

# Subgroup analysis function
run_subgroup_analysis <- function(results, df) {
  # Function to calculate subgroup effects
  calc_subgroup_effect <- function(mask) {
    ate <- average_treatment_effect(results$forest, subset = mask)
    return(c(ate[1], ate[2]))  # Returns estimate and standard error
  }
  
  # HPV status
  print("\nHPV Status Analysis:")
  for (hpv in c(0, 1)) {
    mask <- df$HPV_Positive == hpv
    effect <- calc_subgroup_effect(mask)
    cat(sprintf("HPV %d: Effect = %.3f (95%% CI: %.3f, %.3f), n=%d\n",
                hpv, effect[1], 
                effect[1] - 1.96*effect[2],
                effect[1] + 1.96*effect[2],
                sum(mask)))
  }
  
  # Stage
  print("\nStage Analysis:")
  for (stage in sort(unique(df$Stage_numeric))) {
    mask <- df$Stage_numeric == stage
    if (sum(mask) >= 10) {  # Only analyze if enough samples
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
  # Prepare data
  df <- prepare_data(df)
  
  # Basic data summary
  print("Data Summary:")
  print(summary(df[c("survival_time", "dose_intensity", "event")]))
  
  # Run analysis
  results <- run_causal_analysis(df)
  
  # Print main results
  cat("\nAverage Treatment Effect:", results$ate[1], 
      "\n95% CI: (", results$ate[1] - 1.96*results$ate[2], ",", 
      results$ate[1] + 1.96*results$ate[2], ")\n")
  
  # Run subgroup analysis
  run_subgroup_analysis(results, df)
  
  # Variable importance
  var_imp <- variable_importance(results$forest)
  var_names <- c("Age", "Sex", "Smoking_PY", "Stage", 
                 "HPV+", "HPV Unknown", "Chemo", "RT Year")
  importance_df <- data.frame(
    Variable = var_names,
    Importance = var_imp
  )
  print("\nVariable Importance:")
  print(importance_df[order(-importance_df$Importance), ])
  
  return(results)
}

# Run everything
results <- main_analysis()