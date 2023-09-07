# Title: Comparing_Estimator_Types/multiphase_est.R
# Author: Caleb Leedy
# Date Created: May 3, 2023
# Purpose: This document first compares the two three-phase sampling estimators:
# F3P and K3P (see Notes/multiphase_estimators.tex or page 33 in the research
# journal).

# ******************
# * Script Outline *
# ******************
# 1. Setup
# 2. Compare three-phase estimators

# *************
# * Libraries *
# *************

library(dplyr)

library(doParallel)
library(doRNG)
library(parallelly)

# Run from Simulations/
source("generate_data.R")
source("nmono_est_funs.R")

# ********************
# * Estimator 1: F3P *
# ********************
f3p <- threephase_reg

# ********************
# * Estimator 2: K3P *
# ********************
k3p <- function(df) {

  gfun <- "y2"

  # Steps:
  # 1. Get weights
  # 2. Regression estimator Phase 2
  # 3. Regression estimator Phase 1

  # 1. Get weights
  pi_11 <- df$p12 * df$p1 + df$p21 * df$p2
  df$w <- 1 / pi_11

  df$w_y1 <- 1 / df$p1

  # 2. Regression estimator Phase 2
  df_11 <- filter(df, r1 == 1, r2 == 1)
  mean_y2 <- sum(df_11$w * df_11$y2) / sum(df_11$w)
  mean_x2 <- sum(df_11$w * df_11$x) / sum(df_11$w)
  mean_y1 <- sum(df_11$w * df_11$y1) / sum(df_11$w)

  df_11[[gfun]] <- df_11[[gfun]] - mean_y2
  df_11$x <- df_11$x - mean_x2
  df_11$y1 <- df_11$y1 - mean_y1

  mod <- lm(as.formula(paste0(gfun, " ~ x + y1")), weights = w, data = df_11)

  # 3. Regression estimator Phase 1
  mod_1 <- lm(as.formula(paste0(gfun, " ~ x")), weights = w, data = df_11)

  df_1 <- 
    filter(df, r1 == 1) %>%
    mutate(w_y1 = 1 / p1) 

  mean_x <- mean(df$x)
  mean_x1 <- sum(df_1$x * df_1$w_y1) / sum(df_1$w_y1)
  mean_y11 <- sum(df_1$y1 * df_1$w_y1) / sum(df_1$w_y1)

  mean_y2 + mod$coefficients[2] * (mean_x1 - mean_x2) + 
    mod$coefficients[3] * (mean_y11 - mean_y1) + 
    mod_1$coefficients[2] * (mean_x - mean_x1)

}

# ****************************
# * Comparing the Estimators *
# ****************************

df <- mono_mar(n = 1000)
f3p(df)
k3p(df)

df <- nonmono_mar(n = 1000)
f3p(df)
k3p(df)

# The estimators are different, but how?

# ****************************
# * Default Simulation Study *
# ****************************

# *****************************
# * Monotone Simulation Study *
# *****************************

clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
registerDoParallel(clust)

B <- 1000
true_mean <- 5

mc_theta <- foreach(b = 1:B, 
                  .options.RNG = 1,
                  .packages = c("dplyr", "CVXR")) %dorng% {

  # Parameters
  n_sim <- 1000

  df <- mono_mar(n_sim, mean_y2 = true_mean)
  tibble(oracle = mean(df$y2),
         ipworacle = ipw_mono_est(df, est = FALSE),
         ipwest = ipw_mono_est(df, est = TRUE),
         semi = est_mono(df),
         reg2p = twophase_reg(df),
         reg3p = threephase_reg(df),
         calib = mono_est_weights(df),
         k3p = k3p(df),
         calibw = f3p_cali(df))
         
} %>% bind_rows()

stopCluster(clust)

# ****************
# * Analyze Data *
# ****************

mc_theta %>%
 summarize(
   bias_oracle = mean(oracle) - true_mean,
   bias_ipworacle = mean(ipworacle) - true_mean,
   bias_ipwest = mean(ipwest) - true_mean,
   bias_semi = mean(semi) - true_mean,
   bias_reg2p = mean(reg2p) - true_mean,
   bias_reg3p = mean(reg3p) - true_mean,
   bias_calib = mean(calib) - true_mean,
   bias_k3p = mean(k3p) - true_mean,
   bias_calibw = mean(calibw) - true_mean,
   sd_oracle = sd(oracle),
   sd_ipworacle = sd(ipworacle),
   sd_ipwest = sd(ipwest),
   sd_semi = sd(semi),
   sd_reg2p = sd(reg2p),
   sd_reg3p = sd(reg3p),
   sd_calib = sd(calib),
   sd_k3p = sd(k3p),
   sd_calibw = sd(calibw),
   tstat_oracle = (mean(oracle) - true_mean) / sqrt(var(oracle) / B),
   tstat_ipworacle = (mean(ipworacle) - true_mean) / sqrt(var(ipworacle) / B),
   tstat_ipwest = (mean(ipwest) - true_mean) / sqrt(var(ipwest) / B),
   tstat_semi = (mean(semi) - true_mean) / sqrt(var(semi) / B), 
   tstat_reg2p = (mean(reg2p) - true_mean) / sqrt(var(reg2p) / B), 
   tstat_reg3p = (mean(reg3p) - true_mean) / sqrt(var(reg3p) / B), 
   tstat_calib = (mean(calib) - true_mean) / sqrt(var(calib) / B), 
   tstat_k3p = (mean(k3p) - true_mean) / sqrt(var(k3p) / B),
   tstat_calibw = (mean(calibw) - true_mean) / sqrt(var(calibw) / B)) %>%
  tidyr::pivot_longer(cols = everything(),
               names_to = c(".value", "algorithm"),
               names_pattern = "(.*)_(.*)") %>%
  mutate(pval = pt(-abs(tstat), df = B)) 

#%>%
#  knitr::kable("latex", booktabs = TRUE,
#               digits = 3, caption = paste0("True Value is ", true_mean)) %>%
#  cat(., file = paste0("../Tables/calimono_t" , true_mean, ".tex"))
  
# Overall, ipworacle looks odd, but the others all look good. They work and
# exhibit no bias.

# ********************************
# * Nonmonotone Simulation Study *
# ********************************

run_cali_sims <- function(n_sim, B, true_mean, cor_y1y2,
                          miss_out, est_weights, seed, r_ind_y = FALSE) {


  clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
  registerDoParallel(clust)

  mc_theta <- 
    foreach(b = 1:B, 
           .options.RNG = seed,
           .export = c("nonmono_mar", "ipw_cc_oracle", "est_nonmono", "expit",
                       "nonmono_est_weights", "twophase_reg", "threephase_reg",
                       "k3p", "f3p_cali"),
           .packages = c("dplyr", "CVXR")) %dorng% {

    # Parameters
    df <- nonmono_mar(n_sim, mean_y2 = true_mean, cor_y1y2 = cor_y1y2,
                      r_ind_y = r_ind_y, miss_out = miss_out)

    # Estimation
    # We want to compute:
    # 1. Oracle estimate
    # 2. pi_11 estimate (ipw using fully observed cases)
    # 3. New estimate

    tibble(oracle = mean(df$y2),
           ipworacle = ipw_cc_oracle(df),
           reg2p = twophase_reg(df, reg_on_y1 = FALSE),
           reg3p = threephase_reg(df),
           proposed = est_nonmono(df, oracle = !est_weights),
           calib = nonmono_est_weights(df),
           k3p = k3p(df),
           calibw = f3p_cali(df))

  } %>% bind_rows()

  stopCluster(clust)

  # ****************
  # * Analyze Data *
  # ****************

  mc_theta %>%
   summarize(
     bias_oracle = mean(oracle) - true_mean,
     bias_ipworacle = mean(ipworacle) - true_mean,
     bias_proposed = mean(proposed) - true_mean,
     bias_reg2p = mean(reg2p) - true_mean,
     bias_reg3p = mean(reg3p) - true_mean,
     bias_calib = mean(calib) - true_mean,
     bias_k3p = mean(k3p) - true_mean,
     bias_calibw = mean(calibw) - true_mean,
     sd_oracle = sd(oracle),
     sd_ipworacle = sd(ipworacle),
     sd_proposed = sd(proposed),
     sd_reg2p = sd(reg2p),
     sd_reg3p = sd(reg3p),
     sd_calib = sd(calib),
     sd_k3p = sd(k3p),
     sd_calibw = sd(calibw),
     tstat_oracle = (mean(oracle) - true_mean) / sqrt(var(oracle) / B),
     tstat_ipworacle = (mean(ipworacle) - true_mean) / sqrt(var(ipworacle) / B),
     tstat_proposed = (mean(proposed) - true_mean) / sqrt(var(proposed) / B),
     tstat_reg2p = (mean(reg2p) - true_mean) / sqrt(var(reg2p) / B),
     tstat_reg3p = (mean(reg3p) - true_mean) / sqrt(var(reg3p) / B),
     tstat_calib = (mean(calib) - true_mean) / sqrt(var(calib) / B),
     tstat_k3p = (mean(k3p) - true_mean) / sqrt(var(k3p) / B), 
     tstat_calibw = (mean(calibw) - true_mean) / sqrt(var(calibw) / B)) %>%
    tidyr::pivot_longer(cols = everything(),
                 names_to = c(".value", "algorithm"),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(pval = pt(-abs(tstat), df = B)) 

}

# *******************************
# * Run Nonmonotone Simulations *
# *******************************

tr_m <- 0
cy1y2 <- 0
run_cali_sims(n_sim = 1000, B = 1000,
              true_mean = tr_m, cor_y1y2 = cy1y2,
              miss_out = FALSE, est_weights = FALSE, seed = 2) 

# reg3p and k3p show similar behavior with a significant bias in the estimate.
# calib and proposed work really well. calibw has larger bias and variance but
# is unbiased.


