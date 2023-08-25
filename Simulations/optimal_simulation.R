# Title: optimal_simulation.R
# Created by: Caleb Leedy
# Created on: August 19, 2023
# Purpose: This file contains the code to test the proposed estimator with a
# optimal estimator in the case where everything is normal.

# ***********
# * Outline *
# ***********
# 1. Load Libraries and code.
# 2. Run Monte Carlo simulation.

# *************
# * Libraries *
# *************

library(dplyr)
library(stringr)

library(doParallel)
library(doRNG)
library(parallelly)

source("R/opt_est.R")

# ***************************
# * Monte Carlo Simulations *
# ***************************

clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
registerDoParallel(clust)

B <- 1000
n_obs <- 1000
true_theta <- 0

mc_theta <-
  foreach(iter = 1:B, .options.RNG = 1, .packages = c("dplyr", "stringr")) %dorng% {

    # Generate Data
    df <- gen_optsim_data(n = n_obs, theta = true_theta)

    # Get Estimates
    oracle_est <- mean(df$Y2)
    ipw_est <- 
      filter(df, delta_1 == 1, delta_2 == 1) |>
      mutate(w_est = Y2 / prob_11) |>
      pull(w_est) |>
      sum() / nrow(df)
    cc_est <- 
      filter(df, delta_2 == 1) |>
      pull(Y2) |>
      mean()
    # em_est <- em_optsim(df)
    wls_est <- opt_lin_est(df)
    prop_est <- prop_nmono_est(df)
    prop2_est <- prop2_nmono_est(df)

    return(tibble(oracle = oracle_est,
                  cc = cc_est,
                  ipw = ipw_est,
                  # em = em_est,
                  wls = wls_est,
                  prop = prop_est,
                  prop2 = prop2_est))
  } |> bind_rows()

stopCluster(clust)

# *******************
# * Analyze Results *
# *******************

mc_theta |>
  summarize(
    bias_oracle = mean(oracle) - true_theta,
    bias_cc = mean(cc) - true_theta,
    bias_ipw = mean(ipw) - true_theta,
    bias_wls = mean(wls) - true_theta,
    bias_prop = mean(prop) - true_theta,
    bias_prop2 = mean(prop2) - true_theta,
    sd_oracle = sd(oracle),
    sd_cc = sd(cc),
    sd_ipw = sd(ipw),
    sd_wls = sd(wls),
    sd_prop = sd(prop),
    sd_prop2 = sd(prop2),
    tstat_oracle = (mean(oracle) - true_theta) / sqrt(var(oracle) / B),
    tstat_cc = (mean(cc) - true_theta) / sqrt(var(cc) / B),
    tstat_ipw = (mean(ipw) - true_theta) / sqrt(var(ipw) / B),
    tstat_wls = (mean(wls) - true_theta) / sqrt(var(wls) / B),
    tstat_prop = (mean(prop) - true_theta) / sqrt(var(prop) / B),
    tstat_prop2 = (mean(prop2) - true_theta) / sqrt(var(prop2) / B)
  ) |>
  tidyr::pivot_longer(cols = everything(),
               names_to = c(".value", "algorithm"),
               names_pattern = "(.*)_(.*)") |>
  mutate(pval = pt(-abs(tstat), df = B)) |>
  knitr::kable("latex", booktabs = TRUE,
               digits = 3, caption = paste0("True Value is ", true_theta))
