# Title: Optimal_Nonmonotone_Est/comp_ests_g.R
# Created by: Caleb Leedy
# Created on: September 06, 2023
# Purpose: This file contains the code to test the proposed estimator with a 
# g function of g = "Y1^2Y2".

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

B <- 3000
n_obs <- 1000
true_theta <- 5
cor_e1e2 <- 0.5

mc_theta <-
  foreach(iter = 1:B,
          .options.RNG = 1,
          .packages = c("dplyr", "stringr")) %dorng% {

    # g_fun <- "Y1^2Y2"
    
    # Generate Data
    df <- gen_optsim_data(n = n_obs, theta = true_theta, cor_e1e2 = cor_e1e2)

    # Get Estimates
    oracle_est <- mean(df$Y1^2 * df$Y2)
    oraclex_est <- mean((df$Y1 - df$X)^2 * (df$Y2 - df$X))
    ipw_est <- 
      filter(df, delta_1 == 1, delta_2 == 1) |>
      mutate(g_f = Y1^2 * Y2) |>
      mutate(g_est = g_f / prob_11) |>
      pull(g_est) |>
      sum() / nrow(df)
    cc_est <- 
      filter(df, delta_2 == 1) |>
      mutate(g_f = Y1^2 * Y2) |>
      pull(g_f) |>
      mean()
    # FIXME: Change g_fun
    wls_est <- opt_lin_est(df, cov_y1y2 = cor_e1e2)
    # FIXME: Change g_fun
    prop_est <- prop_nmono_est(df)

    return(tibble(oracle = oracle_est,
                  oraclex = oraclex_est,
                  cc = cc_est,
                  ipw = ipw_est,
                  # em = em_est,
                  wls = wls_est,
                  prop = prop_est))
  } |> 
  bind_rows()

stopCluster(clust)

# *******************
# * Analyze Results *
# *******************

mc_theta |>
  summarize(
    bias_oracle = mean(oracle) - true_theta,
    bias_oraclex = mean(oraclex) - true_theta,
    bias_cc = mean(cc) - true_theta,
    bias_ipw = mean(ipw) - true_theta,
    bias_wls = mean(wls) - true_theta,
    bias_prop = mean(prop) - true_theta,
    bias_prop2 = mean(prop2) - true_theta,
    sd_oracle = sd(oracle),
    sd_oraclex = sd(oraclex),
    sd_cc = sd(cc),
    sd_ipw = sd(ipw),
    sd_wls = sd(wls),
    sd_prop = sd(prop),
    sd_prop2 = sd(prop2),
    tstat_oracle = (mean(oracle) - true_theta) / sqrt(var(oracle) / B),
    tstat_oraclex = (mean(oraclex) - true_theta) / sqrt(var(oraclex) / B),
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
               digits = 3,
               caption = paste0("True Theta is ", true_theta,
                                ". Cov_e1e2 = ", cor_e1e2 ))

# * Testing the standard deviations
# When sigma_11 = 1 and sigma_22 = 1, we have:
# Expected sd: 
# Oracle:
sqrt(2 / n_obs)
# OracleX: 
sqrt(1 / n_obs)
# CC:
sqrt(2 / (n_obs * 0.6))
# IPW:
sqrt((2 + (0.6 * true_theta^2)) / (n_obs * 0.4))
# WLS:
m_mat <- matrix(c(1, 0, 1, 0, 0, 1, 0, 1), ncol = 2)
v_mat <- matrix(c(1 / (n_obs * 0.4), cor_e1e2 / (n_obs * 0.4), 0, 0,
                  cor_e1e2 / (n_obs * 0.4), 1 / (n_obs * 0.4), 0, 0,
                  0, 0, 1 / (n_obs * 0.2), 0,
                  0, 0, 0, 1 / (n_obs * 0.2)), nrow = 4)
exp_inv_cov <- solve(t(m_mat) %*% solve(v_mat) %*% m_mat)
sqrt(exp_inv_cov[2, 2])
# Prop
exp_var <- (1 + cor_e1e2^2 * (1 / 0.4 - 1 / 0.6) + 1 / 0.6 + 
            2 * cor_e1e2^2 * (1 / 0.6 - 1) * 0.4 / 0.6) / n_obs

sqrt(exp_var)
