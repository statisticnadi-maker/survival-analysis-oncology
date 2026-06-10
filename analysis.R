# ==============================================================================
# Project: Advanced Survival Analysis Pipeline for Oncology Research
# Author: Senior Biostatistician & Epidemiologist
# Description: Standard Cox Proportional Hazards and Kaplan-Meier Estimation
#              for Colorectal Cancer (CRC) Epidemiological Data.
# ==============================================================================

# Load required libraries
library(survival)
library(survminer)
library(tidyverse)

# 1. Data Simulation (Emulating Clinical Trial / Cohort Study)
set.seed(42)
n <- 250
crc_data <- tibble(
  patient_id = 1:n,
  age        = round(rnorm(n, mean = 62, sd = 10)),
  gender     = factor(sample(c("Male", "Female"), n, replace = TRUE)),
  stage      = factor(sample(c("Stage III", "Stage IV"), n, replace = TRUE, prob = c(0.6, 0.4))),
  treatment  = factor(sample(c("Standard Chemo", "Targeted Therapy"), n, replace = TRUE)),
  time       = round(rexp(n, rate = 0.02) + runif(n, 1, 10)),
  status     = sample(c(0, 1), n, replace = TRUE, prob = c(0.3, 0.7)) # 1 = Event (Death), 0 = Censored
)

# 2. Kaplan-Meier Survival Analysis
km_fit <- survfit(Surv(time, status) ~ treatment, data = crc_data)

# Print Kaplan-Meier Summary Statistics
print("--- Kaplan-Meier Estimation Summary ---")
print(summary(km_fit)$table)

# 3. Cox Proportional Hazards Regression (Multivariable Model)
cox_model <- coxph(Surv(time, status) ~ treatment + age + gender + stage, data = crc_data)

print("--- Multivariable Cox Proportional Hazards Model ---")
print(summary(cox_model))

# 4. Methodological Diagnostics: Assessing Proportional Hazards (PH) Assumption
# Schoenfeld Residuals Test
ph_test <- cox.zph(cox_model)
print("--- Schoenfeld Residuals Test (PH Assumption) ---")
print(ph_test)

# 5. Output Generation (Simulated Plot Export)
# In production, this saves the publication-ready Kaplan-Meier curve
km_plot <- ggsurvplot(
  km_fit, 
  data = crc_data, 
  pval = TRUE, 
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#E7B800", "#2E9FDF"),
  title = "Overall Survival by Treatment Regimen (Colorectal Cancer Cohort)"
)
