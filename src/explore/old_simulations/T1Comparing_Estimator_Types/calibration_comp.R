# Title: Comparing_Estimator_Types/calibration_comp.R
# Author: Caleb Leedy
# Date Created: April 25, 2023
# Purpose: This script contains the code to compare the calibration estimators
# with the regression estimators. The regression estimators can be originally
# found in Simulations/nonmonotone.R but they have been imported here.

# ******************
# * Script Outline *
# ******************
# A. Monotone Case
#
# B. Nonmonotone Case

# *************
# * Libraries *
# *************

library(dplyr)
library(CVXR)

library(doParallel)
library(doRNG)
library(parallelly)

# *****************************
# * Data Generation Functions *
# *****************************
source("generate_data.R")

# *************************************
# * Non/monotone Estimation Functions *
# *************************************
source("nmono_est_funs.R")

# *******************************
# * Monotone: Current Functions *
# *******************************

# *****************************
# * Monotone Simulation Study *
# *****************************

clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
registerDoParallel(clust)

B <- 1000
true_mean <- 5

mc_theta <- foreach(b = 1:B, 
                  #.options.RNG = 1,
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
         calib = mono_est_weights(df))
         
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
   sd_oracle = sd(oracle),
   sd_ipworacle = sd(ipworacle),
   sd_ipwest = sd(ipwest),
   sd_semi = sd(semi),
   sd_reg2p = sd(reg2p),
   sd_reg3p = sd(reg3p),
   sd_calib = sd(calib),
   tstat_oracle = (mean(oracle) - true_mean) / sqrt(var(oracle) / B),
   tstat_ipworacle = (mean(ipworacle) - true_mean) / sqrt(var(ipworacle) / B),
   tstat_ipwest = (mean(ipwest) - true_mean) / sqrt(var(ipwest) / B),
   tstat_semi = (mean(semi) - true_mean) / sqrt(var(semi) / B), 
   tstat_reg2p = (mean(reg2p) - true_mean) / sqrt(var(reg2p) / B), 
   tstat_reg3p = (mean(reg3p) - true_mean) / sqrt(var(reg3p) / B), 
   tstat_calib = (mean(calib) - true_mean) / sqrt(var(calib) / B)) %>%
  tidyr::pivot_longer(cols = everything(),
               names_to = c(".value", "algorithm"),
               names_pattern = "(.*)_(.*)") %>%
  mutate(pval = pt(-abs(tstat), df = B)) %>%
  knitr::kable("latex", booktabs = TRUE,
               digits = 3, caption = paste0("True Value is ", true_mean)) %>%
  cat(., file = paste0("../Tables/calimono_t" , true_mean, ".tex"))
  
# Overall, pretty good. It looks like this is working as expected.

# ***************************
# * Nonmonotone Calibration *
# ***************************

# **************************
# * Nonmonotone Simulation *
# **************************

run_cali_sims <- function(n_sim, B, true_mean, cor_y1y2,
                          miss_out, est_weights, seed, r_ind_y = FALSE) {


  clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
  registerDoParallel(clust)

  mc_theta <- 
    foreach(b = 1:B, 
           .options.RNG = seed,
           .export = c("nonmono_mar", "ipw_cc_oracle", "est_nonmono", "expit",
                       "nonmono_est_weights", "twophase_reg", "threephase_reg"),
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
           calib = nonmono_est_weights(df))

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
     sd_oracle = sd(oracle),
     sd_ipworacle = sd(ipworacle),
     sd_proposed = sd(proposed),
     sd_reg2p = sd(reg2p),
     sd_reg3p = sd(reg3p),
     sd_calib = sd(calib),
     tstat_oracle = (mean(oracle) - true_mean) / sqrt(var(oracle) / B),
     tstat_ipworacle = (mean(ipworacle) - true_mean) / sqrt(var(ipworacle) / B),
     tstat_proposed = (mean(proposed) - true_mean) / sqrt(var(proposed) / B),
     tstat_reg2p = (mean(reg2p) - true_mean) / sqrt(var(reg2p) / B),
     tstat_reg3p = (mean(reg3p) - true_mean) / sqrt(var(reg3p) / B),
     tstat_calib = (mean(calib) - true_mean) / sqrt(var(calib) / B)) %>%
    tidyr::pivot_longer(cols = everything(),
                 names_to = c(".value", "algorithm"),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(pval = pt(-abs(tstat), df = B)) 

}

# **********************************
# * Simulations: Sim 1 Nonmonotone *
# **********************************

tr_m <- -5
cy1y2 <- 0
run_cali_sims(n_sim = 1000, B = 1000,
              true_mean = tr_m, cor_y1y2 = cy1y2,
              miss_out = FALSE, est_weights = FALSE, seed = 1) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = paste0("True Value is ", tr_m, ". Cor(Y1, Y2) = ", cy1y2)) %>%
  cat(., file = paste0("../Tables/calinonmono1" , "m", tr_m, ".tex"))

tr_m <- 0
cy1y2 <- 0
run_cali_sims(n_sim = 1000, B = 1000,
              true_mean = tr_m, cor_y1y2 = cy1y2,
              miss_out = FALSE, est_weights = FALSE, seed = 2) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = paste0("True Value is ", tr_m, ". Cor(Y1, Y2) = ", cy1y2)) %>%
  cat(., file = paste0("../Tables/calinonmono1" , "m", tr_m, ".tex"))

tr_m <- 5
cy1y2 <- 0
run_cali_sims(n_sim = 1000, B = 1000,
              true_mean = tr_m, cor_y1y2 = cy1y2,
              miss_out = FALSE, est_weights = FALSE, seed = 4) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = paste0("True Value is ", tr_m, ". Cor(Y1, Y2) = ", cy1y2)) %>%
  cat(., file = paste0("../Tables/calinonmono1" , "m", tr_m, ".tex"))




# ****************
# * Clean Tables *
# ****************
# The tables look better if they have the [h!] attribute, so we correct them.
# To be run from Simulations/
source("clean_tables.R")

