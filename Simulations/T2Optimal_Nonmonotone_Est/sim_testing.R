# Title: sim_testing.R
# Created by: Caleb Leedy
# Created on: November 01, 2023
# Purpose: This script tests the different components of the variance estimate 
# of \hat \theta_c. We want to ensure that each component is estimated 
# correctly.

# ***********
# * Outline *
# ***********
# 1. Load Libraries and code.
# 2. Run Monte Carlo simulation for each variance component

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

B <- 5000
n_obs <- 1000
true_theta <- 5
cor_e1e2 <- 0.5

mc_ests <-
  foreach(iter = 1:B,
          .options.RNG = 1,
          .packages = c("dplyr", "stringr")) %dorng% {

    # Generate Data
    df <- gen_optsim_data(n = n_obs, theta = true_theta, cor_e1e2 = cor_e1e2)

    # Estimate different parts of theta_hat_c
    # semidef <- opt_semi_est(df, est = "opt", test = TRUE)
    semidef <- opt_theta_c(df, est = "opt", test = TRUE)

    return(tibble(est = semidef$est,
                  ipw = semidef$ipw,
                  A0 = semidef$A0,
                  A1 = semidef$A1,
                  A2 = semidef$A2,
                  ipwA0 = semidef$ipwA0,
                  ipwA1 = semidef$ipwA1,
                  ipwA2 = semidef$ipwA2,
                  A0A1 = semidef$A0A1,
                  A0A2 = semidef$A0A2,
                  A1A2 = semidef$A1A2
    ))

  } |>
  bind_rows()

stopCluster(clust)

# *******************
# * Analyze Results *
# *******************

# True Values
pi_11 <- 0.4
pi_10 <- 0.2
pi_01 <- 0.2
pi_00 <- 0.2

# For theta_c
c_mat <- 
solve(matrix(c(11.5 * (true_theta^2 + 1),
               -7.5 * (true_theta^2 + 1),
               -7.5 * (true_theta^2 + 1),
               -7.5 * (true_theta^2 + 1),
               7.5 * (true_theta^2 + 1 + cor_e1e2^2),
               2.5 * (true_theta^2 + 1 + cor_e1e2^2),
               -7.5 * (true_theta^2 + 1),
               2.5 * (true_theta^2 + 1 + cor_e1e2^2),
               7.5 * (true_theta^2 + 2)), nrow = 3),
      matrix(c(-3.5 * (true_theta^2 + 1),
               2.5 * (true_theta^2 + 1 + cor_e1e2^2),
               2.5 * (true_theta^2 + 2)), nrow = 3))

c0 <- c_mat[1]
c1 <- c_mat[2]
c2 <- c_mat[3]

var_ipw <- (2 + (1 - pi_11) / pi_11 * (2 + true_theta^2)) / n_obs
var_A0 <- (c0^2 * (true_theta^2 + 1) * 
              (1 / pi_10 + 1 / pi_01 + 1 / pi_11 - 1)) / n_obs
var_A1 <- (c1^2 * (true_theta^2 + 1 + cor_e1e2^2) * 
              (1 / pi_10 + 1 / pi_11)) / n_obs
var_A2 <- (c2^2 * (true_theta^2 + 2) * 
              (1 / pi_01 + 1 / pi_11)) / n_obs
cov_ipwA0 <- (c0 * (true_theta^2 + 1) * (1 + 1 / pi_11)) / n_obs
cov_ipwA1 <- (-c1 * (true_theta^2 + 1 + cor_e1e2^2) * (1 / pi_11)) / n_obs
cov_ipwA2 <- (-c2 * (true_theta^2 + 2) * (1 / pi_11)) / n_obs
cov_A0A1 <- (-(1 / pi_10 + 1 / pi_11) * c0 * c1 * (true_theta^2 + 1)) / n_obs
cov_A0A2 <- (-(1 / pi_01 + 1 / pi_11) * c0 * c2 * (true_theta^2 + 1)) / n_obs
cov_A1A2 <- (1 / pi_11) * c1 * c2 * (true_theta^2 + 1 + cor_e1e2^2) / n_obs

exp_var <- 
  c(var_ipw + var_A0 + var_A1 + var_A2 + 
      2 * (cov_ipwA0 + cov_ipwA1 + cov_ipwA2 + cov_A0A1 + cov_A0A2 + cov_A1A2),
    var_ipw,
    var_A0,
    var_A1,
    var_A2,
    var_ipw + var_A0 + 2 * cov_ipwA0,
    var_ipw + var_A1 + 2 * cov_ipwA1,
    var_ipw + var_A2 + 2 * cov_ipwA2,
    var_A0 + var_A1 + 2 * cov_A0A1,
    var_A0 + var_A2 + 2 * cov_A0A2,
    var_A1 + var_A2 + 2 * cov_A1A2
  )

# For theta_delta
# c_mat <- matrix(c(1, 1, 1), nrow = 3)
V_mat <- 
  matrix(c(pi_00^2 * (1 / pi_11 + 1 / pi_00) * (true_theta^2 + 1),
           pi_00 * pi_10 / pi_11 * (true_theta^2 + 1),
           pi_00 * pi_01 / pi_11 * (true_theta^2 + 1),
           pi_00 * pi_10 / pi_11 * (true_theta^2 + 1),
           pi_10^2 * (1 / pi_11 + 1 / pi_10) * (true_theta^2 + 1 + cor_e1e2^2),
           pi_10 * pi_01 / pi_11 * (true_theta^2 + 1 + cor_e1e2^2),
           pi_00 * pi_01 / pi_11 * (true_theta^2 + 1),
           pi_10 * pi_01 / pi_11 * (true_theta^2 + 1 + cor_e1e2^2),
           pi_01^2 * (1 / pi_11 + 1 / pi_01) * (true_theta^2 + 2)),
  nrow = 3)

C_mat <- matrix(c(-pi_00 / pi_11 * (true_theta^2 + 1),
                  -pi_10 / pi_11 * (true_theta^2 + 1 + cor_e1e2^2),
                  -pi_01 / pi_11 * (true_theta^2 + 2)),
         nrow = 3)

c_mat <- solve(V_mat, C_mat)
c0 <- c_mat[1]
c1 <- c_mat[2]
c2 <- c_mat[3]

var_ipw <- (2 + (1 - pi_11) / pi_11 * (2 + true_theta^2)) / n_obs
var_A0 <- (c0^2 * (true_theta^2 + 1) * pi_00^2 * (1 / pi_11 + 1 / pi_00)) / n_obs
var_A1 <- (c1^2 * (true_theta^2 + 1 + cor_e1e2^2) * 
            pi_10^2 * (1 / pi_10 + 1 / pi_11)) / n_obs
var_A2 <- (c2^2 * (true_theta^2 + 2) * pi_01^2 * 
              (1 / pi_01 + 1 / pi_11)) / n_obs
cov_ipwA0 <- (c0 * (true_theta^2 + 1) * (pi_00 / pi_11)) / n_obs
cov_ipwA1 <- (c1 * (true_theta^2 + 1 + cor_e1e2^2) * (pi_10 / pi_11)) / n_obs
cov_ipwA2 <- (c2 * (true_theta^2 + 2) * (pi_01 / pi_11)) / n_obs
cov_A0A1 <- ((pi_00 * pi_10 / pi_11) * c0 * c1 * (true_theta^2 + 1)) / n_obs
cov_A0A2 <- ((pi_00 * pi_01 / pi_11) * c0 * c2 * (true_theta^2 + 1)) / n_obs
cov_A1A2 <- (pi_10 * pi_01 / pi_11) * c1 * c2 * (true_theta^2 + 1 + cor_e1e2^2) / n_obs

exp_var <- 
  c(var_ipw + var_A0 + var_A1 + var_A2 + 
      2 * (cov_ipwA0 + cov_ipwA1 + cov_ipwA2 + cov_A0A1 + cov_A0A2 + cov_A1A2),
    var_ipw,
    var_A0,
    var_A1,
    var_A2,
    var_ipw + var_A0 + 2 * cov_ipwA0,
    var_ipw + var_A1 + 2 * cov_ipwA1,
    var_ipw + var_A2 + 2 * cov_ipwA2,
    var_A0 + var_A1 + 2 * cov_A0A1,
    var_A0 + var_A2 + 2 * cov_A0A2,
    var_A1 + var_A2 + 2 * cov_A1A2
  )

# Estimated Values
mc_ests |>
  summarize(
    var_est = var(est),
    var_ipw = var(ipw),
    var_A0 = var(A0),
    var_A1 = var(A1),
    var_A2 = var(A2),
    var_ipwA0 = var(ipwA0),
    var_ipwA1 = var(ipwA1),
    var_ipwA2 = var(ipwA2),
    var_A0A1 = var(A0A1),
    var_A0A2 = var(A0A2),
    var_A1A2 = var(A1A2)
  ) |>
  tidyr::pivot_longer(cols = everything(),
                      names_to = c(".value", "algorithm"),
                      names_pattern = "(.*)_(.*)") |>
  mutate(exp_var = exp_var) |>
  knitr::kable(#"latex", booktabs = TRUE,
               digits = 6,
               caption = paste0("True Theta is ", true_theta,
                                ". Cov_e1e2 = ", cor_e1e2 ))

# For theta_c_opt
cov_mat <- 
  matrix(c((true_theta^2 + 1) * ((pi_10 * pi_01 + pi_11^2) / pi_1112 - 1),
           -pi_10 * pi_01 / pi_1112 * (true_theta^2 + 1),
           -pi_10 * pi_01 / pi_1112 * (true_theta^2 + 1),
           -pi_10 * pi_01 / pi_1112 * (true_theta^2 + 1),
           (1 / pi_11 - 1 / (pi_10 + pi_11)) * (true_theta^2 + 1 + cor_e1e2^2),
           pi_10 * pi_01 / pi_1112 * (true_theta^2 + 1 + cor_e1e2^2),
           -pi_10 * pi_01 / pi_1112 * (true_theta^2 + 1),
           pi_10 * pi_01 / pi_1112 * (true_theta^2 + 1 + cor_e1e2^2),
           (1 / pi_11 - 1 / (pi_01 + pi_11)) * (true_theta^2 + 2)),
  nrow = 3)

resp_mat <- 
  matrix(c(-(1 - 1 / (pi_10 + pi_11) - 1 / (pi_01 + pi_11) + 1 / pi_11) * (true_theta^2 + 1),
           -(1 / (pi_10 + pi_11) - 1 / pi_11) * (true_theta^2 + 1 + cor_e1e2^2),
           -(1 / (pi_01 + pi_11) - 1 / pi_11) * (true_theta^2 + 2)),
  nrow = 3)
c_mat <- solve(cov_mat, resp_mat)
#c_mat <- matrix(c(1, 1, 1), nrow = 3)

c0 <- c_mat[1]
c1 <- c_mat[2]
c2 <- c_mat[3]

pi_1112 <- pi_11 * (pi_10 + pi_11) * (pi_01 + pi_11)
var_ipw <- (2 + (1 - pi_11) / pi_11 * (2 + true_theta^2)) / n_obs
var_A0 <- (c0^2 * (true_theta^2 + 1) * 
              ((pi_10 * pi_01 + pi_11^2) / pi_1112 - 1)) / n_obs
var_A1 <- (c1^2 * (true_theta^2 + 1 + cor_e1e2^2) * 
              (1 / pi_11 - 1 / (pi_10 + pi_11))) / n_obs
var_A2 <- (c2^2 * (true_theta^2 + 2) * 
              (1 / pi_11 - 1 / (pi_01 + pi_11))) / n_obs
cov_ipwA0 <- (c0 * (true_theta^2 + 1) * (1 - 1 / (pi_10 + pi_11) - 1 / (pi_01 + pi_11) + 1 / pi_11)) / n_obs
cov_ipwA1 <- (c1 * (true_theta^2 + 1 + cor_e1e2^2) * (1 / (pi_10 + pi_11) - 1 / pi_11)) / n_obs
cov_ipwA2 <- (c2 * (true_theta^2 + 2) * (1 / (pi_01 + pi_11) - 1 / pi_11)) / n_obs
cov_A0A1 <- ((-pi_10 * pi_01 / pi_1112) * c0 * c1 * (true_theta^2 + 1)) / n_obs
cov_A0A2 <- ((-pi_10 * pi_01 / pi_1112) * c0 * c2 * (true_theta^2 + 1)) / n_obs
cov_A1A2 <- (pi_10 * pi_01 / pi_1112) * c1 * c2 * (true_theta^2 + 1 + cor_e1e2^2) / n_obs

exp_var <- 
  c(var_ipw + var_A0 + var_A1 + var_A2 + 
      2 * (cov_ipwA0 + cov_ipwA1 + cov_ipwA2 + cov_A0A1 + cov_A0A2 + cov_A1A2),
    var_ipw,
    var_A0,
    var_A1,
    var_A2,
    var_ipw + var_A0 + 2 * cov_ipwA0,
    var_ipw + var_A1 + 2 * cov_ipwA1,
    var_ipw + var_A2 + 2 * cov_ipwA2,
    var_A0 + var_A1 + 2 * cov_A0A1,
    var_A0 + var_A2 + 2 * cov_A0A2,
    var_A1 + var_A2 + 2 * cov_A1A2
  )

mc_ests |>
  summarize(
    var_est = var(est),
    var_ipw = var(ipw),
    var_A0 = var(A0),
    var_A1 = var(A1),
    var_A2 = var(A2),
    var_ipwA0 = var(ipwA0),
    var_ipwA1 = var(ipwA1),
    var_ipwA2 = var(ipwA2),
    var_A0A1 = var(A0A1),
    var_A0A2 = var(A0A2),
    var_A1A2 = var(A1A2)
  ) |>
  tidyr::pivot_longer(cols = everything(),
                      names_to = c(".value", "algorithm"),
                      names_pattern = "(.*)_(.*)") |>
  mutate(exp_var = exp_var) |>
  knitr::kable(#"latex", booktabs = TRUE,
               digits = 6,
               caption = paste0("True Theta is ", true_theta,
                                ". Cov_e1e2 = ", cor_e1e2 ))
