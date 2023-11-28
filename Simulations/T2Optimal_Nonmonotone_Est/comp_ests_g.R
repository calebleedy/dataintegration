# Title: Optimal_Nonmonotone_Est/comp_ests_g.R
# Created by: Caleb Leedy
# Created on: September 06, 2023
# Purpose: This file contains the code to test the proposed estimator with a 
# g function of g = "Y1^2Y2". This code is very similar to 
# ./optimal_simulation.R except that the g function considered is Y1^2 Y2 
# instead of Y2.

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
    oraclex_est <- mean((df$Y1)^2 * (df$Y2 - df$X))
    ipw_est <- 
      filter(df, delta_1 == 1, delta_2 == 1) |>
      mutate(g_f = Y1^2 * Y2) |>
      mutate(g_est = g_f / prob_11) |>
      pull(g_est) |>
      sum() / nrow(df)
    cc_est <- 
      filter(df, delta_1 == 1, delta_2 == 1) |>
      mutate(g_f = Y1^2 * Y2) |>
      pull(g_f) |>
      mean()
    wls_est <- comb_lin_est(df, gfun = "Y1^2 * Y2", cov_e1e2 = cor_e1e2)
    wlstt_est <- 
      comb_lin_est(df, gfun = "Y1^2 * Y2", cov_e1e2 = cor_e1e2, theta2 = true_theta)
    prop_est <- prop_nmono_est(df, gfun = "Y1^2 * Y2", pow = 3)
    propopt_est <- opt_theta_c(df, gfun = "Y1^2 * Y2", pow = 3)
    semiopt_est <- opt_semi_est(df, gfun = "Y1^2 * Y2", pow = 3)
    semidef_est <- opt_semi_est(df, gfun = "Y1^2 * Y2", est = "default", pow = 3)
    semidel_est <- opt_delta_c(df, gfun = "Y1^2 * Y2", pow = 3)
    semideldef_est <- opt_delta_c(df, gfun = "Y1^2 * Y2", est = "default", pow = 3)

    return(tibble(oracle = oracle_est,
                  oraclex = oraclex_est,
                  cc = cc_est,
                  ipw = ipw_est,
                  wls = wls_est,
                  wlstt = wlstt_est,
                  prop = prop_est,
                  propopt = propopt_est,
                  semiopt = semiopt_est,
                  semidef = semidef_est,
                  semidel = semidel_est,
                  semideldef = semideldef_est))
  } |> 
  bind_rows()

stopCluster(clust)

# *******************
# * Analyze Results *
# *******************
true_g <- 2 * true_theta

mc_theta |>
  summarize(
    bias_oracle = mean(oracle) - true_g,
    bias_oraclex = mean(oraclex) - true_g,
    bias_cc = mean(cc) - true_g,
    bias_ipw = mean(ipw) - true_g,
    bias_wls = mean(wls) - true_g,
    bias_wlstt = mean(wlstt) - true_g,
    bias_prop = mean(prop) - true_g,
    bias_propopt = mean(propopt) - true_g,
    bias_semiopt = mean(semiopt) - true_g,
    bias_semidef = mean(semidef) - true_g,
    bias_semidel = mean(semidel) - true_g,
    bias_semideldef = mean(semideldef) - true_g,
    sd_oracle = sd(oracle),
    sd_oraclex = sd(oraclex),
    sd_cc = sd(cc),
    sd_ipw = sd(ipw),
    sd_wls = sd(wls),
    sd_wlstt = sd(wlstt),
    sd_prop = sd(prop),
    sd_propopt = sd(propopt),
    sd_semiopt = sd(semiopt),
    sd_semidef = sd(semidef),
    sd_semidel = sd(semidel),
    sd_semideldef = sd(semideldef),
    tstat_oracle = (mean(oracle) - true_g) / sqrt(var(oracle) / B),
    tstat_oraclex = (mean(oraclex) - true_g) / sqrt(var(oraclex) / B),
    tstat_cc = (mean(cc) - true_g) / sqrt(var(cc) / B),
    tstat_ipw = (mean(ipw) - true_g) / sqrt(var(ipw) / B),
    tstat_wls = (mean(wls) - true_g) / sqrt(var(wls) / B),
    tstat_wlstt = (mean(wlstt) - true_g) / sqrt(var(wlstt) / B),
    tstat_prop = (mean(prop) - true_g) / sqrt(var(prop) / B),
    tstat_propopt = (mean(propopt) - true_g) / sqrt(var(propopt) / B),
    tstat_semiopt = (mean(semiopt) - true_g) / sqrt(var(semiopt) / B),
    tstat_semidef = (mean(semidef) - true_g) / sqrt(var(semidef) / B),
    tstat_semidel = (mean(semidel) - true_g) / sqrt(var(semidel) / B),
    tstat_semideldef = (mean(semideldef) - true_g) / sqrt(var(semideldef) / B)
  ) |>
  tidyr::pivot_longer(cols = everything(),
               names_to = c(".value", "algorithm"),
               names_pattern = "(.*)_(.*)") |>
  mutate(pval = pt(-abs(tstat), df = B)) |>
  knitr::kable(#"latex", booktabs = TRUE,
               digits = 3,
               caption = paste0("True g is ", true_g,
                                ". Cov_e1e2 = ", cor_e1e2 ))
