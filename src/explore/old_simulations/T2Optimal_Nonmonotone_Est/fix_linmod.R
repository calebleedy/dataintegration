# Title: fix_linmod.R
# Created by: Caleb Leedy
# Created on: December 8, 2023
# Purpose: This file contains the code to test if the optimal linear model 
# contains bias for the non-linear case of g = Y1^2 * Y2.

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
cov_e1e2 <- 0.5

mc_theta <-
  foreach(iter = 1:B,
          .options.RNG = 1,
          .packages = c("dplyr", "stringr")) %dorng% {

    # Generate Data
    df <- gen_optsim_data(n = n_obs, theta = true_theta, cov_e1e2 = cov_e1e2)

    # Get Estimates
    wls_est <- comb_lin_est(df, gfun = "Y1^2 * Y2", cov_e1e2 = cov_e1e2)
    wls_est_test <- 
      comb_lin_est(df, gfun = "Y1^2 * Y2", cov_e1e2 = cov_e1e2, test = TRUE)
    wlstt_est <- 
      comb_lin_est(df, gfun = "Y1^2 * Y2", cov_e1e2 = cov_e1e2, theta2 = true_theta)
    wlstt_est_test <- 
      comb_lin_est(df, gfun = "Y1^2 * Y2", cov_e1e2 = cov_e1e2,
                   theta2 = true_theta, test = TRUE)

    return(tibble(wls = wls_est,
                  wls11 = wls_est_test$est_vec[1],
                  wls10 = wls_est_test$est_vec[2],
                  wls01 = wls_est_test$est_vec[3],
                  wls00 = wls_est_test$est_vec[4],
                  w_11 = wls_est_test$w_vec[1],
                  w_10 = wls_est_test$w_vec[2],
                  w_01 = wls_est_test$w_vec[3],
                  w_00 = wls_est_test$w_vec[4],
                  v_11 = wls_est_test$v_vec[1],
                  v_10 = wls_est_test$v_vec[2],
                  v_01 = wls_est_test$v_vec[3],
                  v_00 = wls_est_test$v_vec[4],
                  wlstt = wlstt_est,
                  wlstt11 = wlstt_est_test$est_vec[1],
                  wlstt10 = wlstt_est_test$est_vec[2],
                  wlstt01 = wlstt_est_test$est_vec[3],
                  wlstt00 = wlstt_est_test$est_vec[4],
                  wtt_11 = wlstt_est_test$w_vec[1],
                  wtt_10 = wlstt_est_test$w_vec[2],
                  wtt_01 = wlstt_est_test$w_vec[3],
                  wtt_00 = wlstt_est_test$w_vec[4]))
  } |> 
  bind_rows()

stopCluster(clust)

# ******************************
# * Check Variance and Weights *
# ******************************

library(ggplot2)

p1 <- 
  mc_theta |>
    mutate(bias_w11 = wls11 - 2 * true_theta) |>
  ggplot(aes(x = bias_w11, y = v_11)) +
    geom_point()

ggsave("Images/biasvar_11.png", p1)

# This shows why we have bias. While we do have 4 unbiased estimators,
# when we weight them together, we use the variance of each estimator 
# for their weight function. Unfortunately, the variance is proportional 
# to the estimated value. This causes us to weight smaller estimates more 
# than larger estimates which lead to a negative bias.

# *******************
# * Analyze Results *
# *******************
true_g <- 2 * true_theta

mc_theta |>
  summarize(
    bias_wls = mean(wls) - true_g,
    bias_wls11 = mean(wls11) - true_g,
    bias_wls10 = mean(wls10) - true_g,
    bias_wls01 = mean(wls01) - true_g,
    bias_wls00 = mean(wls00) - true_g,
    bias_wlstt = mean(wlstt) - true_g,
    bias_wlstt11 = mean(wlstt11) - true_g,
    bias_wlstt10 = mean(wlstt10) - true_g,
    bias_wlstt01 = mean(wlstt01) - true_g,
    bias_wlstt00 = mean(wlstt00) - true_g,
    sd_wls = sd(wls),
    sd_wls11 = sd(wls11),
    sd_wls10 = sd(wls10),
    sd_wls01 = sd(wls01),
    sd_wls00 = sd(wls00),
    sd_wlstt = sd(wlstt),
    sd_wlstt11 = sd(wlstt11),
    sd_wlstt10 = sd(wlstt10),
    sd_wlstt01 = sd(wlstt01),
    sd_wlstt00 = sd(wlstt00),
    tstat_wls = (mean(wls) - true_g) / sqrt(var(wls) / B),
    tstat_wls11 = (mean(wls11) - true_g) / sqrt(var(wls11) / B),
    tstat_wls10 = (mean(wls10) - true_g) / sqrt(var(wls10) / B),
    tstat_wls01 = (mean(wls01) - true_g) / sqrt(var(wls01) / B),
    tstat_wls00 = (mean(wls00) - true_g) / sqrt(var(wls00) / B),
    tstat_wlstt = (mean(wlstt) - true_g) / sqrt(var(wlstt) / B),
    tstat_wlstt11 = (mean(wlstt11) - true_g) / sqrt(var(wlstt11) / B),
    tstat_wlstt10 = (mean(wlstt10) - true_g) / sqrt(var(wlstt10) / B),
    tstat_wlstt01 = (mean(wlstt01) - true_g) / sqrt(var(wlstt01) / B),
    tstat_wlstt00 = (mean(wlstt00) - true_g) / sqrt(var(wlstt00) / B),
  ) |>
  tidyr::pivot_longer(cols = everything(),
               names_to = c(".value", "algorithm"),
               names_pattern = "(.*)_(.*)") |>
  mutate(pval = pt(-abs(tstat), df = B)) |>
  knitr::kable(#"latex", booktabs = TRUE,
               digits = 3,
               caption = paste0("True g is ", true_g,
                                ". Cov_e1e2 = ", cov_e1e2 ))



mean((mc_theta$wls11 + mc_theta$wls10 + mc_theta$wls01 + mc_theta$wls00) / 4)
mean(mc_theta$wls)

