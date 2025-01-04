# if (!require("grf")) install.packages("grf")
# if (!require("survival")) install.packages("survival")
# if (!require("dplyr")) install.packages("dplyr")
# if (!require("ggplot2")) install.packages("ggplot2")

# library(grf)
# library(survival)
# library(dplyr)
# library(ggplot2)

# # Read data
# df <- read.csv("selected_data_with_sites.csv")

# # Data preparation
# prepare_data <- function(df) {
#   # Calculate dose intensity and create necessary variables
#   df <- df %>%
#     mutate(
#       dose_intensity = d_Frac * Fx / Total_days_RT,
#       dose_intensity_scaled = as.vector(scale(dose_intensity)),
#       event = as.numeric(Cause_of_Death_Status == 1),
#       survival_time = Length_FU  # Using Length of Follow-up as survival time
#     )
  
#   return(df)
# }

# # Main analysis function
# run_causal_analysis <- function(df) {
#   # Prepare features
#   X <- df %>%
#     select(Age, Sex, Smoking_PY, Stage_numeric, 
#            HPV_Positive, HPV_Unknown, Chemo, RT_year) %>%
#     as.matrix()
  
#   # Standardize numeric features
#   numeric_cols <- c("Age", "Smoking_PY", "RT_year")
#   X[, numeric_cols] <- scale(X[, numeric_cols])
  
#   # Treatment variable
#   W <- as.vector(df$dose_intensity_scaled)
  
#   # Outcome variable - use survival time weighted by event status
#   Y <- as.vector(df$survival_time * (1 + df$event))  # Weight events more heavily
  
#   # Fit the causal forest
#   cf <- causal_forest(
#     X = X,
#     Y = Y,
#     W = W,
#     num.trees = 2000,
#     honesty = TRUE
#   )
  
#   # Get treatment effects
#   tau.hat <- predict(cf)
  
#   # Calculate average treatment effect
#   ate <- average_treatment_effect(cf)
  
#   # Plot treatment effects
#   p1 <- ggplot(data.frame(effects = tau.hat$predictions), aes(x = effects)) +
#     geom_histogram(bins = 30, fill = "lightblue", color = "black") +
#     geom_vline(xintercept = ate[1], color = "red", linetype = "dashed") +
#     labs(title = "Distribution of Treatment Effects",
#          x = "Effect on Survival Time",
#          y = "Count") +
#     theme_minimal()
#   print(p1)
  
#   # Plot by stage
#   stage_effects <- data.frame(
#     Stage = df$Stage_numeric,
#     Effect = tau.hat$predictions
#   ) %>%
#     group_by(Stage) %>%
#     summarise(
#       mean_effect = mean(Effect),
#       se = sd(Effect) / sqrt(n())
#     )
  
#   p2 <- ggplot(stage_effects, aes(x = factor(Stage), y = mean_effect)) +
#     geom_point() +
#     geom_errorbar(aes(ymin = mean_effect - 1.96*se, 
#                       ymax = mean_effect + 1.96*se), width = 0.2) +
#     labs(title = "Treatment Effects by Stage",
#          x = "Stage",
#          y = "Effect on Survival Time") +
#     theme_minimal()
#   print(p2)
  
#   # Return results
#   return(list(
#     forest = cf,
#     tau.hat = tau.hat,
#     ate = ate
#   ))
# }

# # Subgroup analysis function
# run_subgroup_analysis <- function(results, df) {
#   # Function to calculate subgroup effects
#   calc_subgroup_effect <- function(mask) {
#     ate <- average_treatment_effect(results$forest, subset = mask)
#     return(c(ate[1], ate[2]))  # Returns estimate and standard error
#   }
  
#   # HPV status
#   print("\nHPV Status Analysis:")
#   for (hpv in c(0, 1)) {
#     mask <- df$HPV_Positive == hpv
#     effect <- calc_subgroup_effect(mask)
#     cat(sprintf("HPV %d: Effect = %.3f (95%% CI: %.3f, %.3f), n=%d\n",
#                 hpv, effect[1], 
#                 effect[1] - 1.96*effect[2],
#                 effect[1] + 1.96*effect[2],
#                 sum(mask)))
#   }
  
#   # Stage
#   print("\nStage Analysis:")
#   for (stage in sort(unique(df$Stage_numeric))) {
#     mask <- df$Stage_numeric == stage
#     if (sum(mask) >= 10) {  # Only analyze if enough samples
#       effect <- calc_subgroup_effect(mask)
#       cat(sprintf("Stage %d: Effect = %.3f (95%% CI: %.3f, %.3f), n=%d\n",
#                   stage, effect[1],
#                   effect[1] - 1.96*effect[2],
#                   effect[1] + 1.96*effect[2],
#                   sum(mask)))
#     }
#   }
# }

# # Main execution
# main_analysis <- function() {
#   # Prepare data
#   df <- prepare_data(df)
  
#   # Basic data summary
#   print("Data Summary:")
#   print(summary(df[c("survival_time", "dose_intensity", "event")]))
  
#   # Run analysis
#   results <- run_causal_analysis(df)
  
#   # Print main results
#   cat("\nAverage Treatment Effect:", results$ate[1], 
#       "\n95% CI: (", results$ate[1] - 1.96*results$ate[2], ",", 
#       results$ate[1] + 1.96*results$ate[2], ")\n")
  
#   # Run subgroup analysis
#   run_subgroup_analysis(results, df)
  
#   # Variable importance
#   var_imp <- variable_importance(results$forest)
#   var_names <- c("Age", "Sex", "Smoking_PY", "Stage", 
#                  "HPV+", "HPV Unknown", "Chemo", "RT Year")
#   importance_df <- data.frame(
#     Variable = var_names,
#     Importance = var_imp
#   )
#   print("\nVariable Importance:")
#   print(importance_df[order(-importance_df$Importance), ])
  
#   return(results)
# }

# results <- main_analysis()


if (!require("grf")) install.packages("grf")
if (!require("survival")) install.packages("survival")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("survminer")) install.packages("survminer")
if (!require("cmprsk")) install.packages("cmprsk")

library(grf)
library(survival)
library(dplyr)
library(ggplot2)
library(survminer)
library(cmprsk)

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
  
  # Plot treatment effects distribution with white background
  p1 <- ggplot(data.frame(effects = tau.hat$predictions), aes(x = effects)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    geom_vline(xintercept = ate[1], color = "red", linetype = "dashed") +
    labs(title = "Distribution of Treatment Effects",
         x = "Effect on Survival Time",
         y = "Count") +
    theme_bw()
  ggsave("treatment_effects_distribution.png", p1)
  
  # Plot by stage with white background
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
    theme_bw()
  ggsave("treatment_effects_by_stage.png", p2)
  
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

# Plot survival curves by treatment (dose intensity)
plot_survival_by_dose_intensity <- function(df) {
  fit <- survfit(Surv(survival_time, event) ~ dose_intensity, data = df)
  ggsurv <- ggsurvplot(fit, data = df, risk.table = TRUE, pval = TRUE, 
                       title = "Survival Curves by Dose Intensity")
  ggsave("survival_by_dose_intensity.png", ggsurv$plot + theme_bw())
}

# Plot cumulative incidence function for competing risks
plot_cif <- function(df) {
  cr_model <- cuminc(ftime = df$survival_time, fstatus = df$event)
  plot(cr_model, main = "Cumulative Incidence Function", xlab = "Time", ylab = "CIF")
  ggsave("cif_plot.png", plot = last_plot() + theme_bw())
}

# Hazard ratio plot using Cox model
plot_hazard_ratios <- function(df) {
  cox_model <- coxph(Surv(survival_time, event) ~ Age + HPV_Positive + Smoking_PY + Chemo, data = df)
  cox_summary <- summary(cox_model)
  hr_df <- data.frame(
    Variable = rownames(cox_summary$coefficients),
    HR = exp(cox_summary$coefficients[, 1]),
    LowerCI = exp(cox_summary$conf.int[, 1]),
    UpperCI = exp(cox_summary$conf.int[, 3])
  )
  p <- ggplot(hr_df, aes(x = Variable, y = HR)) +
    geom_point() +
    geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI)) +
    labs(title = "Hazard Ratios for Covariates", x = "Variable", y = "Hazard Ratio") +
    theme_bw()
  ggsave("hazard_ratios.png", p)
}

# Causal Forest variable importance plot
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
    labs(title = "Variable Importance in Causal Forest", x = "Variable", y = "Importance") +
    theme_bw()
  ggsave("variable_importance.png", p)
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
  
  # Plot additional visualizations
  plot_survival_by_dose_intensity(df)
  plot_cif(df)
  plot_hazard_ratios(df)
  plot_variable_importance(results)
  
  return(results)
}

results <- main_analysis()
