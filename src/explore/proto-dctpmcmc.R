# Title: proto-dctpmcmc.R
# Created by: Caleb Leedy
# Created on: April 24, 2024
# Purpose: This script runs the code for MCMC analysis of two-phase samples.

# *************
# * Libraries *
# *************

library(dplyr)
library(doParallel)
library(doRNG)

source("proto-dctp.R")

# **************
# * Parameters *
# **************

B_sims <- 1000
N_obs <- 10000

p1_type <- "srs"
# p1_pi <- "ifelse(pnorm(-pop_df$X2) < 0.10, 0.10, pnorm(-pop_df$X2))"
p1_pi <- "1200 / nrow(pop_df)"
p2_type <- "poisson"
# p2_pi <- "ifelse(pt(-pop_df$X1 + 1, 3) > 0.9, 0.9, pt(-pop_df$X1 + 1, 3))"
p2_pi <- "pmin(pmax(pt(pop_df$Z - 1, 3), 0.02), 0.7)"
comb <- paste0(p1_type, "-", p2_type)

# ************
# * Run MCMC *
# ************

set.seed(1)
pop_df <- gen_pop(N_obs)
true_theta <- mean(pop_df$Y)

upop_df <- update_pop(pop_df, p1_pi, p2_pi)

clust <- makeCluster(min(detectCores() - 2, 100), outfile = "")
registerDoParallel(clust)

mc_res <- 
  foreach(iter=1:B_sims,
          .packages = c("dplyr", "nleqslv", "rlang")) %dopar% {

  if (iter %% 100 == 0) {
    print(paste0("Iter: ", iter))
  }

  set.seed(iter)

  samps <- gen_samps(upop_df, p1_type, p2_type)
  p1_df <- samps[[1]]
  p2_df <- samps[[2]]

  pistar <- pi_star(p1_df, p2_df, N_obs, sampling = comb)
  tpreg <- tp_reg(p1_df, p2_df, upop_df, sampling = comb)
  dcpop <- dc_ybar(p1_df, p2_df, upop_df,
                   qi = 1, entropy = "EL", estT1 = FALSE, sampling = comb)
  dcest <- dc_ybar(p1_df, p2_df, upop_df,
                   qi = 1, entropy = "EL", estT1 = TRUE, sampling = comb)
  etapop <- 
    dc_ybar(p1_df, p2_df, upop_df,
            qi = 1, entropy = "EL", estT1 = FALSE, sampling = comb, linearized = TRUE)
  etaest <- 
    dc_ybar(p1_df, p2_df, upop_df,
            qi = 1, entropy = "EL", estT1 = TRUE, sampling = comb, linearized = TRUE)

  return(tibble(
  Est = c("pistar", "tpreg", "dcpop", "dcest", "etapop", "etaest"),
  Theta = c(pistar[[1]], tpreg[[1]], dcpop[[1]], dcest[[1]], etapop[[1]], etaest[[1]]),
  Var = c(pistar[[2]], tpreg[[2]], dcpop[[2]], dcest[[2]], etapop[[2]], etaest[[2]]),
  Iter = iter)
  )

} %>% bind_rows()

stopCluster(clust)

# ****************
# * Analyze Data *
# ****************

# Find any failed estimates
which(is.na(mc_res$Theta))
which(is.na(mc_res$Var))
mc_res %>%
  mutate(iter = rep(1:B_sims, each = 6)) %>%
  filter(is.na(Theta))

# Analyze Mean
mc_res %>%
  mutate(err = Theta - true_theta) %>%
  mutate(CI = abs(err) < qnorm(0.975) * sqrt(Var)) %>%
  group_by(Est) %>%
  summarize(Bias = mean(err, na.rm = TRUE),
            RMSE = sqrt(mean(err^2, na.rm = TRUE)),
            EmpCI = mean(CI, na.rm = TRUE),
            sdest = sd(Theta, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Ttest = abs(Bias) / sqrt(sdest^2 / B_sims)) %>%
  select(-sdest) %>%
  mutate(Ind = case_when(
    Est == "etaest" ~ 6,
    Est == "etapop" ~ 5,
    Est == "dcest" ~ 4,
    Est == "dcpop" ~ 3,
    Est == "tpreg" ~ 2,
    Est == "pistar" ~ 1)
  ) %>%
  mutate(Est = case_when(
    Est == "etaest" ~ "Eta-Est",
    Est == "etapop" ~ "Eta-Pop",
    Est == "dcest" ~ "DC-Est",
    Est == "dcpop" ~ "DC-Pop",
    Est == "tpreg" ~ "TP-Reg",
    Est == "pistar" ~ "$\\pi^*$")
  ) %>%
  arrange(Ind) %>%
  select(Est, Bias, RMSE, EmpCI, Ttest) %>%
  knitr::kable("latex", booktabs = TRUE, escape = FALSE, digits = 3) %>%
  cat(file = paste0("tables/tpdcsim",  "_mean.tex"))

# Analyze Variance
mc_res %>%
  group_by(Est) %>%
  summarize(MCVar = var(Theta, na.rm = TRUE),
            EstVar = mean(Var, na.rm = TRUE),
            VarVar = var(Var, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Ttest = abs(MCVar - EstVar) / sqrt(VarVar / B_sims)) %>%
  mutate(Ind = case_when(
    Est == "etaest" ~ 6,
    Est == "etapop" ~ 5,
    Est == "dcest" ~ 4,
    Est == "dcpop" ~ 3,
    Est == "tpreg" ~ 2,
    Est == "pistar" ~ 1)
  ) %>%
  mutate(Est = case_when(
    Est == "etaest" ~ "Eta-Est",
    Est == "etapop" ~ "Eta-Pop",
    Est == "dcest" ~ "DC-Est",
    Est == "dcpop" ~ "DC-Pop",
    Est == "tpreg" ~ "TP-Reg",
    Est == "pistar" ~ "$\\pi^*$")
  ) %>%
  arrange(Ind) %>%
  select(Est, MCVar, EstVar, VarVar, Ttest) 

%>%
  knitr::kable("latex", booktabs = TRUE, escape = FALSE, digits = 4) %>%
  cat(file = paste0("tables/tpdcsim", "_var.tex"))


mc_res %>%
  filter(Est == "dcpop") %>%
  pull(Var) %>% mean()
%>%
  arrange(desc(Var))

mc_res %>%
  filter(Est == "dcpop") %>%
  arrange(desc(Theta))

