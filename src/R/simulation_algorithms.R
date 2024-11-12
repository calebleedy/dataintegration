# File: simulation_algorithms.R
# Created by: Caleb Leedy
# Created on: November 4, 2024
# Purpose:
# This file contains all of the simulation algorithms that are used in the
# paper simulations procedure found in ../final_sims.R.

# *********************************
# * Two-Phase Sampling Algorithms *
# *********************************
# Exploration file: ../explore/proto-dctp.R

#' This function generates a population data frame with the requisite variables.
#'
#' @param obs - This is the number of observations in the population.
#'
#' @details - This function does not include any sampling probabilities because
#' it does not take into account the sampling type or formula.
tp_pop <- function(obs, vars) {
  delta <- vars$delta

  x1 <- rnorm(obs, 2, 1)
  x2 <- runif(obs, 0, 4)
  x3 <- rnorm(obs, 0, 1)
  x4 <- runif(obs, 0.1, 0.9)
  z <- rnorm(obs, 0, 1)
  eps <- rnorm(obs)

  y <- 3 * x1 + 2 * x2 + delta * 0.5 * z + eps
  n1 <- 1200
  p1_pi <- n1 / obs
  p2_pi <- pmin(pmax(pt(z - 1, 3), 0.02), 0.7)

  return(tibble(X1 = x1, X2 = x2, X3 = x3, X4 = x4, Z = z, Y = y,
                pi1 = p1_pi, pi2 = p2_pi))
}

#' This function generates Phase 1 and Phase 2 samples from an updated
#' population with selection probabilites.
tp_samps <- function(uppop_df, vars = list("srs", "poisson")) {
  p1_type <- vars[[1]]
  p2_type <- vars[[2]]

  if (p1_type == "poisson") {
    del1 <- rbinom(nrow(uppop_df), 1, uppop_df$pi1)
  } else if (p1_type == "srs") {
    ind <- sample(1:nrow(uppop_df), size = uppop_df$pi1[1] * nrow(uppop_df))
    del1 <- as.numeric(1:nrow(uppop_df) %in% ind)
  } else {
    stop("We have only implemented poisson and srs for Phase 1.")
  }

  p1_df <- mutate(uppop_df, del1 = del1) %>% filter(del1 == 1)

  if (p2_type == "poisson") {
    del2 <- rbinom(nrow(p1_df), 1, p1_df$pi2)
  } else if (p2_type == "srs") {
    ind <- sample(1:nrow(p1_df), size = round(p1_df$pi2[1] * nrow(p1_df)))
    del2 <- as.numeric(1:nrow(p1_df) %in% ind)
  } else {
    stop("We have only implemented poisson and srs for Phase 2.")
  }

  p2_df <- mutate(p1_df, del2 = del2) %>% filter(del2 == 1)

  return(list(p1_df, p2_df))
}


# *******************************************
# * Nonnested Two-Phase Sampling Algorithms *
# *******************************************
# Exploration file: ../explore/20240425-nndcsim.qmd

nn_pop <- function(obs, vars) {
  delta <- vars$delta

  x1 <- rnorm(obs, 2, 1)
  x2 <- runif(obs, 0, 4)
  x3 <- rnorm(obs, 0, 1)
  z <- rnorm(obs, 0, 1)
  x4 <- runif(obs, 0.1, 0.9)
  eps <- rnorm(obs)

  y <- 3 * x1 + 2 * x2 + delta * z + eps

  n1 <- 1000
  pi1 <- n1 / obs
  pi2 <- pmin(pmax(pt(z - 2.5, 3), 0.01), 0.9)

  return(tibble(X1 = x1, X2 = x2, X3 = x3, X4 = x4, Z = z, Y = y,
                pi1 = pi1, pi2 = pi2))
}

#' This function generates Phase 1 and Phase 2 samples from an updated
#' population with selection probabilites.
nn_samps <- function(upop_df, vars = list("srs", "poisson")) {
  p1_type <- vars[[1]]
  p2_type <- vars[[2]]

  if (p1_type == "poisson") {
    del1 <- rbinom(nrow(upop_df), 1, upop_df$pi1)
  } else if (p1_type == "srs") {
    ind <- sample(1:nrow(upop_df), size = upop_df$pi1[1] * nrow(upop_df))
    del1 <- as.numeric(1:nrow(upop_df) %in% ind)
  } else {
    stop("We have only implemented poisson and srs for Phase 1.")
  }

  p1_df <- mutate(upop_df, del1 = del1) %>% filter(del1 == 1)

  if (p2_type == "poisson") {
    del2 <- rbinom(nrow(upop_df), 1, upop_df$pi2)
  } else if (p2_type == "srs") {
    ind <- sample(1:nrow(upop_df), size = round(upop_df$pi2[1] * nrow(upop_df)))
    del2 <- as.numeric(1:nrow(upop_df) %in% ind)
  } else {
    stop("We have only implemented poisson and srs for Phase 2.")
  }

  p2_df <- mutate(upop_df, del2 = del2) %>% filter(del2 == 1)

  return(list(p1_df, p2_df))
}


# ***********************************
# * Multisource Sampling Algorithms *
# ***********************************
# Exploration file: ../explore/20240506-msdcsim.qmd, ../explore/20240511-msdcsim2.qmd

ms_pop <- function(N_obs, vars) {
  delta <- vars$delta
  n1_obs <- vars$n1_obs

  x1 <- rnorm(N_obs, 2, 1)
  x2 <- runif(N_obs, 0, 4)
  x3 <- rnorm(N_obs, 5, 1)
  z <- rnorm(N_obs, 0, 1)

  eps <- rnorm(N_obs)

  y <- 3 * x1 + 2 * x2 + delta * z + eps

  pi0 <- pmin(pmax(pnorm(z - 2), 0.02), 0.9)
  # pi1 <- pnorm(-z - 1) # Informative sampling. Do we want this?
  pi1 <- n1_obs / N_obs
  # pi2 <- pmin(pmax(pnorm(x2 - 4), 0.01), 0.9)
  pi2 <- pnorm(x2 - 2)

  return(tibble(X1 = x1, X2 = x2, X3 = x3, Y = y,
                pi0 = pi0, pi1 = pi1, pi2 = pi2))
}

#' This function generates Phase 1 and Phase 2 samples from an updated
#' population with selection probabilites.
ms_samps <- function(pop_df, vars = list("poisson", "srs", "poisson")) {

  p0_type <- vars[[1]]
  p1_type <- vars[[2]]
  p2_type <- vars[[3]]

  if (p0_type == "poisson") {
    del0 <- rbinom(nrow(pop_df), 1, pop_df$pi0)
  } else if (p0_type == "srs") {
    ind <- sample(1:nrow(pop_df), size = pop_df$pi0[1] * nrow(pop_df))
    del0 <- as.numeric(1:nrow(pop_df) %in% ind)
  } else {
    stop("We have only implemented poisson and srs for Sample 0.")
  }

  p0_df <- mutate(pop_df, del0 = del0) %>% filter(del0 == 1)

  if (p1_type == "poisson") {
    del1 <- rbinom(nrow(pop_df), 1, pop_df$pi1)
  } else if (p1_type == "srs") {
    ind <- sample(1:nrow(pop_df), size = pop_df$pi1[1] * nrow(pop_df))
    del1 <- as.numeric(1:nrow(pop_df) %in% ind)
  } else {
    stop("We have only implemented poisson and srs for Sample 1.")
  }

  p1_df <- mutate(pop_df, del1 = del1) %>% filter(del1 == 1)

  if (p2_type == "poisson") {
    del2 <- rbinom(nrow(pop_df), 1, pop_df$pi2)
  } else if (p2_type == "srs") {
    ind <- sample(1:nrow(pop_df), size = round(pop_df$pi2[1] * nrow(pop_df)))
    del2 <- as.numeric(1:nrow(pop_df) %in% ind)
  } else {
    stop("We have only implemented poisson and srs for Sample 2.")
  }

  p2_df <- mutate(pop_df, del2 = del2) %>% filter(del2 == 1)

  return(list(p0_df, p1_df, p2_df))
}


# ********************************************************
# * Multisource with Estimated Alpha Sampling Algorithms *
# ********************************************************
# Exploration file: ../explore/20241101-msdcsim_estalpha.qmd

msa_pop <- function(N_obs, vars) {
  delta <- vars$delta
  n1_obs <- vars$n1_obs

  x1 <- rnorm(N_obs, 2, 1)
  x2 <- runif(N_obs, 0, 4)
  x3 <- rnorm(N_obs, 5, 1)
  z <- rnorm(N_obs, 0, 1)

  eps <- rnorm(N_obs)

  y <- 3 * x1 + 2 * x2 + z + eps

  pi0 <- pmin(pmax(pnorm(z - 2), 0.02), 0.9)
  # pi1 <- pnorm(-z - 1) # Informative sampling. Do we want this?
  pi1 <- n1_obs / N_obs
  # pi2 <- pmin(pmax(pnorm(x2 - 4), 0.01), 0.9)
  pi2 <- pnorm(x2 - 2)

  return(tibble(X1 = x1, X2 = x2, X3 = x3, Y = y,
                pi0 = pi0, pi1 = pi1, pi2 = pi2))
}

#' This function generates Phase 1 and Phase 2 samples from an updated
#' population with selection probabilites.
msa_samps <- function(pop_df, vars = list("poisson", "srs", "poisson")) {
  p0_type <- vars[[1]]
  p1_type <- vars[[2]]
  p2_type <- vars[[3]]

  if (p0_type == "poisson") {
    del0 <- rbinom(nrow(pop_df), 1, pop_df$pi0)
  } else if (p0_type == "srs") {
    ind <- sample(1:nrow(pop_df), size = pop_df$pi0[1] * nrow(pop_df))
    del0 <- as.numeric(1:nrow(pop_df) %in% ind)
  } else {
    stop("We have only implemented poisson and srs for Sample 0.")
  }

  p0_df <- mutate(pop_df, del0 = del0) %>% filter(del0 == 1)

  if (p1_type == "poisson") {
    del1 <- rbinom(nrow(pop_df), 1, pop_df$pi1)
  } else if (p1_type == "srs") {
    ind <- sample(1:nrow(pop_df), size = pop_df$pi1[1] * nrow(pop_df))
    del1 <- as.numeric(1:nrow(pop_df) %in% ind)
  } else {
    stop("We have only implemented poisson and srs for Sample 1.")
  }

  p1_df <- mutate(pop_df, del1 = del1) %>% filter(del1 == 1)

  if (p2_type == "poisson") {
    del2 <- rbinom(nrow(pop_df), 1, pop_df$pi2)
  } else if (p2_type == "srs") {
    ind <- sample(1:nrow(pop_df), size = round(pop_df$pi2[1] * nrow(pop_df)))
    del2 <- as.numeric(1:nrow(pop_df) %in% ind)
  } else {
    stop("We have only implemented poisson and srs for Sample 2.")
  }

  p2_df <- mutate(pop_df, del2 = del2) %>% filter(del2 == 1)

  return(list(p0_df, p1_df, p2_df))
}

