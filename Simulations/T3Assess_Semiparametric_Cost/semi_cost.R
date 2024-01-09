# Title: semi_cost.R
# Created by: Caleb Leedy
# Created on: January 9, 2024
# Purpose: This file contains the code to compare the semiparametric 
# estimators with parametric estimators. It largely follows the same 
# design as ../T2Optimal_Nonmonotone_Est/optimal_simulation.R.

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

    # Generate Data
    df <- gen_optsim_data(n = n_obs, theta = true_theta, cor_e1e2 = cor_e1e2)

    # Get Estimates
    oracle_est <- mean(df$Y2)
    oraclex_est <- mean(df$Y2 - df$X)
    ipw_est <- 
      filter(df, delta_1 == 1, delta_2 == 1) |>
      mutate(g_est = Y2 / prob_11) |>
      pull(g_est) |>
      sum() / nrow(df)
    cc_est <- 
      filter(df, delta_2 == 1) |>
      pull(Y2) |>
      mean()
    # em_est <- em_optsim(df)
    wls_est <- opt_lin_est(df, cov_y1y2 = cor_e1e2)
    wlsalt_est <- comb_lin_est_lin(df, theta2 = true_theta, cov_e1e2 = cor_e1e2)
    prop_est <- prop_nmono_est(df)
    propind_est <- prop_nmono_est(df, prop_ind = TRUE)
    propopt_est <- opt_theta_c(df)
    propoptdef_est <- opt_theta_c(df, est = "default")
    semiopt_est <- opt_semi_est(df)
    semidef_est <- opt_semi_est(df, est = "default")
    semidel_est <- opt_delta_c(df)
    semideldef_est <- opt_delta_c(df, est = "default")

    return(tibble(oracle = oracle_est,
                  oraclex = oraclex_est,
                  cc = cc_est,
                  ipw = ipw_est,
                  # em = em_est,
                  wls = wls_est,
                  wlsalt = wlsalt_est,
                  prop = prop_est,
                  propind = propind_est,
                  propopt = propopt_est,
                  propoptdef = propoptdef_est,
                  semiopt = semiopt_est,
                  semidef = semidef_est,
                  semidel = semidel_est,
                  semideldef = semideldef_est
    ))
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
    bias_wlsalt = mean(wlsalt) - true_theta,
    bias_prop = mean(prop) - true_theta,
    bias_propind = mean(propind) - true_theta,
    bias_propopt = mean(propopt) - true_theta,
    bias_propoptdef = mean(propoptdef) - true_theta,
    bias_semiopt = mean(semiopt) - true_theta,
    bias_semidef = mean(semidef) - true_theta,
    bias_semidel = mean(semidel) - true_theta,
    bias_semideldef = mean(semideldef) - true_theta,
    sd_oracle = sd(oracle),
    sd_oraclex = sd(oraclex),
    sd_cc = sd(cc),
    sd_ipw = sd(ipw),
    sd_wls = sd(wls),
    sd_wlsalt = sd(wlsalt),
    sd_prop = sd(prop),
    sd_propind = sd(propind),
    sd_propopt = sd(propopt),
    sd_propoptdef = sd(propoptdef),
    sd_semiopt = sd(semiopt),
    sd_semidef = sd(semidef),
    sd_semidel = sd(semidel),
    sd_semideldef = sd(semideldef),
    tstat_oracle = (mean(oracle) - true_theta) / sqrt(var(oracle) / B),
    tstat_oraclex = (mean(oraclex) - true_theta) / sqrt(var(oraclex) / B),
    tstat_cc = (mean(cc) - true_theta) / sqrt(var(cc) / B),
    tstat_ipw = (mean(ipw) - true_theta) / sqrt(var(ipw) / B),
    tstat_wls = (mean(wls) - true_theta) / sqrt(var(wls) / B),
    tstat_wlsalt = (mean(wlsalt) - true_theta) / sqrt(var(wlsalt) / B),
    tstat_prop = (mean(prop) - true_theta) / sqrt(var(prop) / B),
    tstat_propind = (mean(propind) - true_theta) / sqrt(var(propind) / B),
    tstat_propopt = (mean(propopt) - true_theta) / sqrt(var(propopt) / B),
    tstat_propoptdef = (mean(propoptdef) - true_theta) / sqrt(var(propoptdef) / B),
    tstat_semiopt = (mean(semiopt) - true_theta) / sqrt(var(semiopt) / B),
    tstat_semidef = (mean(semidef) - true_theta) / sqrt(var(semidef) / B),
    tstat_semidel = (mean(semidel) - true_theta) / sqrt(var(semidel) / B),
    tstat_semideldef = (mean(semideldef) - true_theta) / sqrt(var(semideldef) / B)
  ) |>
  tidyr::pivot_longer(cols = everything(),
                      names_to = c(".value", "algorithm"),
                      names_pattern = "(.*)_(.*)") |>
  mutate(pval = pt(-abs(tstat), df = B)) |>
  knitr::kable(#"latex", booktabs = TRUE,
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
exp_var <- 
  (1 + cor_e1e2^2 * (1 / 0.4 - 1 / 0.6) + 1 / 0.6 + 
   2 * cor_e1e2^2 * (1 / 0.6 - 1) * 0.4 / 0.6) / n_obs

sqrt(exp_var)
# Semiopt
## Default
sqrt((6 + 7.5 * cor_e1e2^2) / n_obs)
