# Title: analyze_semiconstraints.R
# Created by: Caleb Leedy
# Created on: January 17, 2024
# Purpose: This file aims to explore the results of `semi_cost.R` in more
# detail. I do not understand why adding a constraint does not change the 
# estimation in the first simulation.

# *************
# * Libraries *
# *************

library(dplyr)
library(tidyr)
library(stringr)
library(CVXR)
library(ggplot2)

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

B <- 100
n_obs <- 1000
true_theta <- 5 
cov_e1e2 <- 0.5

mc_theta <-
  foreach(iter = 1:B,
          .options.RNG = 1,
          .packages = c("dplyr", "stringr", "CVXR", "tidyr")) %dorng% {

    # Generate Data
    # df <- gen_optsim_data(n = n_obs, theta = true_theta, cov_e1e2 = cov_e1e2)
    df <- gen_unbal_data(n = n_obs, theta = true_theta, cov_e1e2 = cov_e1e2)

    # Get Estimates
    conalg_para <- min_var_alg(df, mod = "para", cov_e1e2 = cov_e1e2)
    conalg_sumzero <- min_var_alg(df, mod = "sumzero", cov_e1e2 = cov_e1e2)
    conalg_outrob <- min_var_alg(df, mod = "out_rob", cov_e1e2 = cov_e1e2)
    conalg_resprob <- min_var_alg(df, mod = "resp_rob", cov_e1e2 = cov_e1e2)
    conalg_doubrob <- min_var_alg(df, mod = "double_rob", cov_e1e2 = cov_e1e2)

    return(
      tibble(conpara = conalg_para$theta_est,
             conzero = conalg_sumzero$theta_est,
             conout = conalg_outrob$theta_est,
             conresp = conalg_resprob$theta_est,
             condoub = conalg_doubrob$theta_est,
             parac = conalg_para$c_hat,
             zeroc = conalg_sumzero$c_hat,
             outc = conalg_outrob$c_hat,
             respc = conalg_resprob$c_hat,
             doubc = conalg_doubrob$c_hat,
             iter = iter,
             c_i = 1:9
      ) |> 
      group_by(conpara, conzero, conout, conresp, condoub) |> 
      nest()
    )
  } |>
  bind_rows()

stopCluster(clust)

# ****************
# * Analyze Data *
# ****************

c_tab <-
  bind_rows(mc_theta$data) |>
  pivot_longer(cols = parac:doubc, names_to = "estimator", values_to = "c_val") |>
  mutate(c_i = paste0("c_", c_i)) |>
  pivot_wider(names_from = "c_i", values_from = "c_val")

# Variance of each coefficient per estimator
c_tab |>
  group_by(estimator) |>
  summarize(across(starts_with("c_"), var), .groups = "drop")


get_f_statistic <- function(c_tab, reduced_mod, full_mod) {

}
