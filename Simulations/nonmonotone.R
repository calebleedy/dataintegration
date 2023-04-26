# Title: nonmonotone.R
# Author: Caleb Leedy
# Date Created: March 29, 2023
# Purpose: This document contains the simulations for computing estimating
# equations with both monotone and nonmonotone missingness. See
# Notes/Non-monotone.pdf for more information.
#
# Questions:
# 1. What estimating functions should we use? 
# - U(\theta, Y) = (\bar y_1 + \bar y_2) / 2 - \theta
# - U(\theta, Y) = \bar y_1 + \bar y_2 - \theta
# - U(\theta, Y) = \bar y_1 - \bar y_2 - \theta
# Other ideas include moments of a distribution of quantiles.
#
# 2. How should nonmonotone missingness data be generated?
# In the monotone setting we think about data being generated conditionally (or
# at least the missingness is conditionally generated). This is because we want
# to have data MAR instead of MCAR. The challenge in this is that we also have
# to think about how to ensure that 
#
# \[p(R_1, R_2 | X) = p(R_1 | R_2, X) p(R_2 | X) = p(R_2 | R_1, X) p(R_1 | X).\]
# 
# A further challenge is that these conditional probabilites can depend on Y_1
# and Y_2 if they have been observed. However, we cannot predetermine if an
# observation should have wlog R_2 depend on Y_1 because we don't know if Y_1 is
# observed for that element.
#
# Answers:
# 1. Currently, we will use U(\theta, Y) = \bar y_2 - \theta.
#
# 2. While one cannot have a global model this is MAR, using the ideas of Robins
# and Gill (1998), we can construct an MAR model where each probability is
# conditional.

# *************
# * Libraries *
# *************

library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

library(doParallel)
library(doRNG)
library(parallelly)

# *****************
# * Generate Data *
# *****************

# Run from Simulations/
source("generate_data.R")

# ************************
# * Estimating Functions *
# ************************

source("nmono_est_funs.R")

# ****************************************
# * Run Simulation: Monotone Missingness *
# ****************************************

clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
registerDoParallel(clust)

B <- 1000
true_mean <- -5

mc_theta <- foreach(b = 1:B, 
                  #.options.RNG = 1,
                  .packages = c("dplyr")) %dorng% {

  # Parameters
  n_sim <- 1000

  df <- mono_mar(n_sim, mean_y2 = true_mean)
  tibble(oracle = mean(df$y2),
         ipworacle = ipw_mono_est(df, est = FALSE),
         ipwest = ipw_mono_est(df, est = TRUE),
         semi = est_mono(df))
         
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
   sd_oracle = sd(oracle),
   sd_ipworacle = sd(ipworacle),
   sd_ipwest = sd(ipwest),
   sd_semi = sd(semi),
   tstat_oracle = (mean(oracle) - true_mean) / sqrt(var(oracle) / B),
   tstat_ipworacle = (mean(ipworacle) - true_mean) / sqrt(var(ipworacle) / B),
   tstat_ipwest = (mean(ipwest) - true_mean) / sqrt(var(ipwest) / B),
   tstat_semi = (mean(semi) - true_mean) / sqrt(var(semi) / B)) %>%
  tidyr::pivot_longer(cols = everything(),
               names_to = c(".value", "algorithm"),
               names_pattern = "(.*)_(.*)") %>%
  mutate(pval = pt(-abs(tstat), df = B)) %>%
  knitr::kable("latex", booktabs = TRUE,
               digits = 3, caption = paste0("True Value is ", true_mean)) %>%
  cat(., file = paste0("../Tables/monomar_t" , true_mean, ".tex"))
  

# *******************************************
# * Run Simulation: Nonmonotone Missingness *
# *******************************************

# FIXME: miss_resp is misleading. This should be a parameter of nonmono_mar
# instead it is a parameter for est_nonmono. This should be changed to something
# like est_weights.
run_init_sims <- function(n_sim, B, true_mean, cor_y1y2,
                          miss_out, miss_resp, seed, r_ind_y = FALSE) {


  clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
  registerDoParallel(clust)

  mc_theta <- 
    foreach(b = 1:B, 
            .options.RNG = seed,
            .export = c("nonmono_mar", "ipw_cc_oracle", "est_nonmono", "expit"),
            .packages = c("dplyr")) %dorng% {

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
           proposed = est_nonmono(df, oracle = !miss_resp))

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
     sd_oracle = sd(oracle),
     sd_ipworacle = sd(ipworacle),
     sd_proposed = sd(proposed),
     tstat_oracle = (mean(oracle) - true_mean) / sqrt(var(oracle) / B),
     tstat_ipworacle = (mean(ipworacle) - true_mean) / sqrt(var(ipworacle) / B),
     tstat_proposed = (mean(proposed) - true_mean) / sqrt(var(proposed) / B))%>%
    tidyr::pivot_longer(cols = everything(),
                 names_to = c(".value", "algorithm"),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(pval = pt(-abs(tstat), df = B)) 

}

# **********************************
# * Simulations: Sim 1 Nonmonotone *
# **********************************

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = -5, cor_y1y2 = 0,
              miss_out = FALSE, miss_resp = FALSE, seed = 1) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is -5. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/nonmonosim1" , "m-5", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0,
              miss_out = FALSE, miss_resp = FALSE, seed = 1) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/nonmonosim1" , "m0", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 5, cor_y1y2 = 0,
              miss_out = FALSE, miss_resp = FALSE, seed = 1) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 5. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/nonmonosim1" , "m5", ".tex"))

# **********************************
# * Simulations: Sim 2 Nonmonotone *
# **********************************

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.1,
              miss_out = FALSE, miss_resp = FALSE, seed = 2) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.1") %>%
  cat(., file = paste0("../Tables/nonmonosim2" , "c0.1", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.5,
              miss_out = FALSE, miss_resp = FALSE, seed = 2) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.5") %>%
  cat(., file = paste0("../Tables/nonmonosim2" , "c0.5", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.9,
              miss_out = FALSE, miss_resp = FALSE, seed = 2) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.9") %>%
  cat(., file = paste0("../Tables/nonmonosim2" , "c0.9", ".tex"))

# **********************************
# * Simulations: Sim 3 Nonmonotone *
# **********************************

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0,
              miss_out = TRUE, miss_resp = FALSE, seed = 3) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/nonmonosim3" , "c0", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.1,
              miss_out = TRUE, miss_resp = FALSE, seed = 3) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.1") %>%
  cat(., file = paste0("../Tables/nonmonosim3" , "c0.1", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.5,
              miss_out = TRUE, miss_resp = FALSE, seed = 3) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.5") %>%
  cat(., file = paste0("../Tables/nonmonosim3" , "c0.5", ".tex"))

# **********************************
# * Simulations: Sim 4 Nonmonotone *
# **********************************

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0,
              miss_out = FALSE, miss_resp = TRUE, seed = 4, r_ind_y = TRUE) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/nonmonosim4" , "c0", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.1,
              miss_out = FALSE, miss_resp = TRUE, seed = 4, r_ind_y = TRUE) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.1") %>%
  cat(., file = paste0("../Tables/nonmonosim4" , "c0.1", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.5,
              miss_out = FALSE, miss_resp = TRUE, seed = 4, r_ind_y = TRUE) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.5") %>%
  cat(., file = paste0("../Tables/nonmonosim4" , "c0.5", ".tex"))

# **********************************
# * Simulations: Sim 5 Nonmonotone *
# **********************************

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0,
              miss_out = FALSE, miss_resp = TRUE, seed = 4) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/nonmonosim5" , "c0", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.1,
              miss_out = FALSE, miss_resp = TRUE, seed = 4) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.1") %>%
  cat(., file = paste0("../Tables/nonmonosim5" , "c0.1", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.5,
              miss_out = FALSE, miss_resp = TRUE, seed = 4) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.5") %>%
  cat(., file = paste0("../Tables/nonmonosim5" , "c0.5", ".tex"))

# %>%
#   knitr::kable("latex", booktabs = TRUE, digits = 3,
#   caption = paste0("True Value is ", true_mean, ". Cor(Y1, Y2) = ", cor_y1y2)) %>%
#   cat(., file = paste0("../Tables/nonmonocmar_t" , true_mean,
#                        "cy1y2", cor_y1y2, ".tex"))

# *********************************
# * Run Simulation: Estimate pi_i *
# *********************************

run_resp_sims <- function(n_sim, B, true_mean, cor_y1y2,
                          miss_out, r_ind_y, seed) {


  clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
  registerDoParallel(clust)

  mc_theta <- 
    foreach(b = 1:B, 
            .options.RNG = seed,
            .export = c("nonmono_mar", "ipw_cc_oracle",
                        "ipw_cc_est", "est_nonmono", "expit"),
            .packages = c("dplyr")) %dorng% {

    # Parameters
    df <- nonmono_mar(n_sim, mean_y2 = true_mean, cor_y1y2 = cor_y1y2,
                      r_ind_y = r_ind_y, miss_out = miss_out,
                      miss_resp = FALSE, mcar = TRUE)

    # Estimation
    # We want to compute:
    # 1. Oracle estimate
    # 2. pi_11 estimate (ipw using fully observed cases)
    # 3. New estimate

    tibble(oracle = mean(df$y2),
           ipw.or = ipw_cc_oracle(df),
           prop.0 = est_nonmono(df, oracle = FALSE, alpha = 0),
           ipw.0 = ipw_cc_est(df, oracle = FALSE, alpha = 0),
           ipw.half = ipw_cc_est(df, oracle = FALSE, alpha = 0.5),
           ipw.1 = ipw_cc_est(df, oracle = FALSE, alpha = 1))

  } %>% bind_rows()

  stopCluster(clust)

  # ****************
  # * Analyze Data *
  # ****************

  mc_theta %>%
   summarize(
     bias_oracle = mean(oracle) - true_mean,
     bias_ipw.or= mean(ipw.or) - true_mean,
     bias_prop.0 = mean(prop.0) - true_mean,
     bias_ipw.0 = mean(ipw.0) - true_mean,
     bias_ipw.half = mean(ipw.half) - true_mean,
     bias_ipw.1 = mean(ipw.1) - true_mean,
     sd_oracle = sd(oracle),
     sd_ipw.or = sd(ipw.or),
     sd_prop.0 = sd(prop.0),
     sd_ipw.0 = sd(ipw.0),
     sd_ipw.half = sd(ipw.half),
     sd_ipw.1 = sd(ipw.1),
     tstat_oracle = (mean(oracle) - true_mean) / sqrt(var(oracle) / B),
     tstat_ipw.or= (mean(ipw.or) - true_mean) / sqrt(var(ipw.or) / B),
     tstat_prop.0 = (mean(prop.0) - true_mean) / sqrt(var(prop.0) / B),
     tstat_ipw.0 = (mean(ipw.0) - true_mean) / sqrt(var(ipw.0) / B),
     tstat_ipw.half = (mean(ipw.half) - true_mean) / sqrt(var(ipw.half) / B),
     tstat_ipw.1 = (mean(ipw.1) - true_mean) / sqrt(var(ipw.1) / B)) %>%
    tidyr::pivot_longer(cols = everything(),
                 names_to = c(".value", "algorithm"),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(pval = pt(-abs(tstat), df = B)) 

}


run_resp_sims(n_sim = 2000, B = 2000,
             true_mean = 0, cor_y1y2 = 0,
             miss_out = FALSE, r_ind_y = TRUE, seed = 4) %>%
  filter(algorithm %in% c("oracle", "ipw.or", "ipw.0", "prop.0")) 


%>%
  knitr::kable("latex", booktabs = TRUE,
               digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/ipwsim" , ".tex"))

# ******************************************************
# * Run Simulation: Nonmonotone Missingness with alpha *
# ******************************************************

run_alpha_sims <- function(n_sim, B, true_mean, cor_y1y2,
                          miss_out, seed) {


  clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
  registerDoParallel(clust)

  mc_theta <- 
    foreach(b = 1:B, 
            .options.RNG = seed,
            .export = c("nonmono_mar", "ipw_cc_oracle", "est_nonmono", "expit"),
            .packages = c("dplyr")) %dorng% {

    # Parameters
    df <- nonmono_mar(n_sim, mean_y2 = true_mean, cor_y1y2 = cor_y1y2,
                      r_ind_y = FALSE, miss_out = miss_out)

    # Estimation
    # We want to compute:
    # 1. Oracle estimate
    # 2. pi_11 estimate (ipw using fully observed cases)
    # 3. New estimate

    tibble(oracle = mean(df$y2),
           ipworacle = ipw_cc_oracle(df),
           prop.or = est_nonmono(df, oracle = TRUE),
           prop.0 = est_nonmono(df, oracle = FALSE, alpha = 0),
           prop.half = est_nonmono(df, oracle = FALSE, alpha = 0.5),
           prop.1 = est_nonmono(df, oracle = FALSE, alpha = 1))

  } %>% bind_rows()

  stopCluster(clust)

  # ****************
  # * Analyze Data *
  # ****************

  mc_theta %>%
   summarize(
     bias_oracle = mean(oracle) - true_mean,
     bias_ipworacle = mean(ipworacle) - true_mean,
     bias_prop.or = mean(prop.or) - true_mean,
     bias_prop.0 = mean(prop.0) - true_mean,
     bias_prop.half = mean(prop.half) - true_mean,
     bias_prop.1 = mean(prop.1) - true_mean,
     sd_oracle = sd(oracle),
     sd_ipworacle = sd(ipworacle),
     sd_prop.or = sd(prop.or),
     sd_prop.0 = sd(prop.0),
     sd_prop.half = sd(prop.half),
     sd_prop.1 = sd(prop.1),
     tstat_oracle = (mean(oracle) - true_mean) / sqrt(var(oracle) / B),
     tstat_ipworacle = (mean(ipworacle) - true_mean) / sqrt(var(ipworacle) / B),
     tstat_prop.or = (mean(prop.or) - true_mean) / sqrt(var(prop.or) / B),
     tstat_prop.0 = (mean(prop.0) - true_mean) / sqrt(var(prop.0) / B),
     tstat_prop.half = (mean(prop.half) - true_mean) / sqrt(var(prop.half) / B),
     tstat_prop.1 = (mean(prop.1) - true_mean) / sqrt(var(prop.1) / B)) %>%
    tidyr::pivot_longer(cols = everything(),
                 names_to = c(".value", "algorithm"),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(pval = pt(-abs(tstat), df = B)) 

}

run_alpha_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0,
              miss_out = FALSE, seed = 5) %>%
  knitr::kable("latex", booktabs = TRUE,
               digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/alphasim" , ".tex"))


# ****************
# * Clean Tables *
# ****************
# The tables look better if they have the [h!] attribute, so we correct them.
# To be run from Simulations/
source("clean_tables.R")

