# File: final_sims.R
# Created by: Caleb Leedy
# Created on: November 4, 2024
# Purpose:
# This file contains the code to run all of the simulations for the Data
# Integration paper. It also saves the output of all of the tables that are used
# in the final paper.
#
# The functions in this script have been developed and ran in the previous
# exploration documents:
# * ./explore/20240411-tpdcsim.qmd
# * ./explore/20240425-nndcsim.qmd
# * ./explore/20240506-msdcsim.qmd
# * ./explore/20240511-msdcsim2.qmd
# * ./explore/20241025-mcbssim.qmd
# * ./explore/20241101-msdcsim_estalpha.qmd
# * ./explore/proto-dctp.R

# **************
# * Parameters *
# **************

c_opts <- commandArgs(trailingOnly = TRUE)
if (length(c_opts) > 0) {
  sim_setup <- c_opts[1]
  pvars <- list(delta = as.numeric(c_opts[2]))
} else {
  # Options: "twophasesamp", "nonnestedtp", "multisource", "multisourcealpha"
  sim_setup <- "multisource"
  pvars <- list(delta = 1)
}


# Number of Monte Carlo replicates
B_sims <- 1000

# Finite population size
N_obs <- 10000

# ********************************
# * Custom Simulation Parameters *
# ********************************

if (sim_setup == "twophasesamp") {
  # delta controls model misspecification.
  # If delta = 1, the model is misspecified.
  name_add <- ifelse(pvars$delta == 1, "d1", "d0")
} else if (sim_setup == "nonnestedtp") {
  name_add <- ifelse(pvars$delta == 1, "d1", "d0")
} else if (sim_setup == "multisource") {
  pvars$n1_obs <- 2000
  name_add <- ifelse(pvars$delta == 1, "d1", "d0")
} else if (sim_setup == "multisourcealpha") {
  pvars$n1_obs <- 2000
  name_add <- ""
} else {

}

# *************
# * Libraries *
# *************

library(dplyr)
library(doParallel)
library(doRNG)
library(nleqslv)

# **************
# * Algorithms *
# **************

source("R/debiased_calibration.R")
source("R/integration_algorithms.R")
source("R/alternative_algorithms.R")
source("R/simulation_algorithms.R")

# ********************
# * Simulation Setup *
# ********************

if (sim_setup == "twophasesamp") {
  gen_pop <- tp_pop
  gen_samps <- tp_samps
} else if (sim_setup == "nonnestedtp") {
  gen_pop <- nn_pop
  gen_samps <- nn_samps
} else if (sim_setup == "multisource") {
  gen_pop <- ms_pop
  gen_samps <- ms_samps
} else if (sim_setup == "multisourcealpha") {
  gen_pop <- msa_pop
  gen_samps <- msa_samps
} else {

}


# ********************
# * Estimation Setup *
# ********************

set.seed(1)
pop_df <- gen_pop(N_obs, pvars)
true_theta <- mean(pop_df$Y)

clust <- makeCluster(min(detectCores() - 2, 100), outfile = "")
registerDoParallel(clust)

mc_res <- 
  foreach(iter=1:B_sims, .packages = c("dplyr", "nleqslv"), .options.RNG = 1) %dorng% {

    if (iter %% 100 == 0) {
      print(paste0("Iter: ", iter))
    }

    samps <- gen_samps(pop_df)

    if (sim_setup == "twophasesamp") {
      # Data
      p1_df <- samps[[1]]
      p2_df <- samps[[2]]
      alg_v <- list(p1_df, p2_df, pop_df, estT1 = TRUE)
      alg_v1 <- list(p1_df, p2_df, pop_df, estT1 = FALSE)

      # Algorithms
      pistar <- tp_pistar(alg_v)
      reg <- tp_reg(alg_v)
      est <- tp_est(alg_v)
      estp <- tp_est(alg_v1)

      # Returned tibble
      ret <- tibble(Est = c("PiStar", "Reg", "EstPop", "Est"),
                    Theta = c(pistar$theta, reg$theta, estp$theta, est$theta), 
                    Var = c(pistar$vhat, reg$vhat, estp$vhat, est$vhat), 
                    Iter = iter)

    } else if (sim_setup == "nonnestedtp") {

      # Data
      p1_df <- samps[[1]]
      p2_df <- samps[[2]]
      alg_v <- list(p1_df, p2_df, pop_df, estT1 = TRUE)
      alg_v1 <- list(p1_df, p2_df, pop_df, estT1 = FALSE)

      # Algorithms
      ht <- nn_ht(alg_v)
      reg <- nn_reg(alg_v)
      est <- nn_est(alg_v)
      estp <- nn_est(alg_v1)

      # Returned tibble
      ret <- tibble(Est = c("HT", "Reg", "EstPop", "Est"),
                    Theta = c(ht$theta, reg$theta, estp$theta, est$theta), 
                    Var = c(ht$vhat, reg$vhat, estp$vhat, est$vhat), 
                    Iter = iter)

    } else if (sim_setup == "multisource") {
     
      # Data
      p0_df <- samps[[1]]
      p1_df <- samps[[2]]
      p2_df <- samps[[3]]
      alg_v <- list(p0_df, p1_df, p2_df, pop_df, estT1 = TRUE)
      alg_v1 <- list(p0_df, p1_df, p2_df, pop_df, estT1 = FALSE)

      # Algorithms
      ht <- ms_ht(alg_v)
      reg <- ms_reg(alg_v)
      est <- ms_est(alg_v)
      estp <- ms_est(alg_v1)

      # Returned tibble
      ret <- tibble(Est = c("HT", "NNReg", "EstPop", "Est"),
                    Theta = c(ht$theta, reg$theta, estp$theta, est$theta), 
                    Var = c(ht$vhat, reg$vhat, estp$vhat, est$vhat), 
                    Iter = iter)

    } else if (sim_setup == "multisourcealpha") {

      # Data
      p0_df <- samps[[1]]
      p1_df <- samps[[2]]
      p2_df <- samps[[3]]
      alg_v <- list(p0_df, p1_df, p2_df, pop_df)

      # Algorithms
      ht <- msa_ht(alg_v)
      reg <- msa_reg(alg_v)
      est <- suppressWarnings(msa_est(alg_v))

      # Returned tibble
      ret <- tibble(Est = c("HT", "Reg", "Est"),
                    Theta = c(ht$theta, reg$theta, est$theta), 
                    Var = c(ht$vhat, reg$vhat, est$vhat), 
                    Iter = iter)
    } else {

    }

    return(ret)

  } %>% bind_rows()

stopCluster(clust)

# *******************
# * Analyze Results *
# *******************

mc_res %>%
  mutate(err = Theta - true_theta) %>%
  mutate(CI = abs(err) < qnorm(0.975) * sqrt(Var)) %>%
  group_by(Est) %>%
  summarize(Bias = mean(err),
            SE = sd(err),
            RMSE = sqrt(mean(err^2)),
            EmpCI = mean(CI),
            MCVar = var(Theta),
            EstVar = mean(Var),
            sdest = sd(Theta),
            RelBias = (EstVar - MCVar) / MCVar) %>%
  ungroup() %>%
  mutate(Ttest = abs(Bias) / sqrt(sdest^2 / B_sims)) %>%
  select(Est, Bias, SE, RMSE, EmpCI, Ttest, MCVar, EstVar, RelBias) %>%
  knitr::kable("latex", booktabs = TRUE,
              digits = c(0, 3, 3, 3, 3, 2, 3, 3, 3)) %>%
  cat(file = paste0("tables/", sim_setup, name_add, ".tex"))

