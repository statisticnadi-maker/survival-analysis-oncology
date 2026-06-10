# Advanced Survival Analysis Pipeline for Oncology Research

This repository contains a production-ready statistical and epidemiological pipeline developed in R for conducting survival analysis and hazard modeling in clinical oncology datasets, with a focus on colorectal cancer (CRC) cohorts.

## Methodological Architecture
1. **Non-parametric Estimation:** Kaplan-Meier survival curves with log-rank testing for treatment efficacy.
2. **Semi-parametric Modeling:** Multivariable Cox Proportional Hazards regression estimating Adjusted Hazard Ratios (aHR) controlled for age, gender, and clinical stage.
3. **Statistical Diagnostics:** Evaluation of the Proportional Hazards assumption utilizing Schoenfeld residuals ($cox.zph$).

## Dependencies
Ensure you have R (version $\ge$ 4.0) and the following libraries installed:
- `survival`
- `survminer`
- `tidyverse`

## Usage
Execute the primary analytical pipeline via:
```R
source("analysis.R")
