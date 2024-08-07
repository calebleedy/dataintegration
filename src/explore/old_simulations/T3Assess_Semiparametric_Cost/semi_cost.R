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
library(CVXR)

library(doParallel)
library(doRNG)
library(parallelly)

source("R/opt_est.R")
source("Simulations/T3Assess_Semiparametric_Cost/compare_algs.R")

# ***************************
# * Monte Carlo Simulations *
# ***************************

# ************************
# * theta = E[g] = E[Y2] *
# ************************

clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
registerDoParallel(clust)

B <- 1000
n_obs <- 1000
true_theta <- 5 
cov_e1e2 <- 0.5

mc_theta <-
  foreach(iter = 1:B,
          .options.RNG = 1,
          .packages = c("dplyr", "stringr", "CVXR")) %dorng% {

    # Generate Data
    # df <- gen_optsim_data(n = n_obs, theta = true_theta, cov_e1e2 = cov_e1e2)
    df <- gen_unbal_data(n = n_obs, theta = true_theta, cov_e1e2 = cov_e1e2)

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
    wls_est <- opt_lin_est(df, cov_y1y2 = cov_e1e2)
    wlsalt_est <- comb_lin_est_lin(df, theta2 = true_theta, cov_e1e2 = cov_e1e2)
    prop_est <- prop_nmono_est(df)
    propind_est <- prop_nmono_est(df, prop_ind = TRUE)
    propopt_est <- opt_theta_c(df)
    propoptdef_est <- opt_theta_c(df, est = "default")
    semiopt_est <- opt_semi_est(df)
    semidef_est <- opt_semi_est(df, est = "default")
    semidel_est <- opt_delta_c(df)
    semideldef_est <- opt_delta_c(df, est = "default")

    conalg_para <- min_var_alg(df, mod = "para", cov_e1e2 = cov_e1e2)
    conalg_sumzero <- min_var_alg(df, mod = "sumzero", cov_e1e2 = cov_e1e2)
    conalg_outrob <- min_var_alg(df, mod = "out_rob", cov_e1e2 = cov_e1e2)
    conalg_resprob <- min_var_alg(df, mod = "resp_rob", cov_e1e2 = cov_e1e2)
    conalg_doubrob <- min_var_alg(df, mod = "double_rob", cov_e1e2 = cov_e1e2)

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
                  semideldef = semideldef_est,
                  conpara = conalg_para$theta_est,
                  conzero = conalg_sumzero$theta_est,
                  conout = conalg_outrob$theta_est,
                  conresp = conalg_resprob$theta_est,
                  condoub = conalg_doubrob$theta_est
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
    bias_conpara = mean(conpara) - true_theta,
    bias_conzero = mean(conzero) - true_theta,
    bias_conout = mean(conout) - true_theta,
    bias_conresp = mean(conresp) - true_theta,
    bias_condoub = mean(condoub) - true_theta,
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
    sd_conpara = sd(conpara),
    sd_conzero = sd(conzero),
    sd_conout = sd(conout),
    sd_conresp = sd(conresp),
    sd_condoub = sd(condoub),
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
    tstat_semideldef = (mean(semideldef) - true_theta) / sqrt(var(semideldef) / B),
    tstat_conpara = (mean(conpara) - true_theta) / sqrt(var(conpara) / B),
    tstat_conzero = (mean(conzero) - true_theta) / sqrt(var(conzero) / B),
    tstat_conout = (mean(conout) - true_theta) / sqrt(var(conout) / B),
    tstat_conresp = (mean(conresp) - true_theta) / sqrt(var(conresp) / B),
    tstat_condoub = (mean(condoub) - true_theta) / sqrt(var(condoub) / B)
  ) |>
  tidyr::pivot_longer(cols = everything(),
                      names_to = c(".value", "algorithm"),
                      names_pattern = "(.*)_(.*)") |>
  mutate(pval = pt(-abs(tstat), df = B)) |>
  knitr::kable(#"latex", booktabs = TRUE,
               digits = 3,
               caption = paste0("True Theta is ", true_theta,
                                ". Cov_e1e2 = ", cov_e1e2 ))

# *******************************
# * theta = E[g] = E[Y1^2 * Y2] *
# *******************************

clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
registerDoParallel(clust)

mc_theta <-
  foreach(iter = 1:B,
          .options.RNG = 1,
          .packages = c("dplyr", "stringr", "CVXR")) %dorng% {

    # g_fun <- "Y1^2Y2"
    
    # Generate Data
    # df <- gen_optsim_data(n = n_obs, theta = true_theta, cov_e1e2 = cov_e1e2)
    df <- gen_unbal_data(n = n_obs, theta = true_theta, cov_e1e2 = cov_e1e2)

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
    wls_est <- comb_lin_est(df, gfun = "Y1^2 * Y2", cov_e1e2 = cov_e1e2)
    wlstt_est <- 
      comb_lin_est(df, gfun = "Y1^2 * Y2", cov_e1e2 = cov_e1e2, theta2 = true_theta)
    prop_est <- prop_nmono_est(df, gfun = "Y1^2 * Y2", pow = 3)
    propopt_est <- opt_theta_c(df, gfun = "Y1^2 * Y2", pow = 3)
    semiopt_est <- opt_semi_est(df, gfun = "Y1^2 * Y2", pow = 3)
    semidef_est <- opt_semi_est(df, gfun = "Y1^2 * Y2", est = "default", pow = 3)
    semidel_est <- opt_delta_c(df, gfun = "Y1^2 * Y2", pow = 3)
    semideldef_est <- opt_delta_c(df, gfun = "Y1^2 * Y2", est = "default", pow = 3)

    conalg_para <- 
      min_var_alg(df, gfun = "Y1^2 * Y2", mod = "para", cov_e1e2 = cov_e1e2)
    conalg_sumzero <- 
      min_var_alg(df, gfun = "Y1^2 * Y2", mod = "sumzero", cov_e1e2 = cov_e1e2)
    conalg_outrob <- 
      min_var_alg(df, gfun = "Y1^2 * Y2", mod = "out_rob", cov_e1e2 = cov_e1e2)
    conalg_resprob <- 
      min_var_alg(df, gfun = "Y1^2 * Y2", mod = "resp_rob", cov_e1e2 = cov_e1e2)
    conalg_doubrob <- 
      min_var_alg(df, gfun = "Y1^2 * Y2", mod = "double_rob", cov_e1e2 = cov_e1e2)

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
                  semideldef = semideldef_est,
                  conpara = conalg_para$theta_est,
                  conzero = conalg_sumzero$theta_est,
                  conout = conalg_outrob$theta_est,
                  conresp = conalg_resprob$theta_est,
                  condoub = conalg_doubrob$theta_est
    ))
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
    bias_conpara = mean(conpara) - true_g,
    bias_conzero = mean(conzero) - true_g,
    bias_conout = mean(conout) - true_g,
    bias_conresp = mean(conresp) - true_g,
    bias_condoub = mean(condoub) - true_g,
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
    sd_conpara = sd(conpara),
    sd_conzero = sd(conzero),
    sd_conout = sd(conout),
    sd_conresp = sd(conresp),
    sd_condoub = sd(condoub),
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
    tstat_semideldef = (mean(semideldef) - true_g) / sqrt(var(semideldef) / B),
    tstat_conpara = (mean(conpara) - true_g) / sqrt(var(conpara) / B),
    tstat_conzero = (mean(conzero) - true_g) / sqrt(var(conzero) / B),
    tstat_conout = (mean(conout) - true_g) / sqrt(var(conout) / B),
    tstat_conresp = (mean(conresp) - true_g) / sqrt(var(conresp) / B),
    tstat_condoub = (mean(condoub) - true_g) / sqrt(var(condoub) / B)
  ) |>
  tidyr::pivot_longer(cols = everything(),
               names_to = c(".value", "algorithm"),
               names_pattern = "(.*)_(.*)") |>
  mutate(pval = pt(-abs(tstat), df = B)) |>
  knitr::kable(#"latex", booktabs = TRUE,
               digits = 3,
               caption = paste0("True g is ", true_g,
                                ". Cov_e1e2 = ", cov_e1e2 ))

