# Title: ../Simulations/nmono_est_funs.R
# Author: Caleb Leedy
# Date Created: April 25, 2023
# Purpose: This script contains the functions used to estimate monotone and
# nonmonotone missingness. Like Simulations/generate_data.R, I wanted to
# separate the function estimation from the simulation study itself in order to
# gain more modular code.

# ------------------------------------------------------------------------------

# ********************
# * Script Functions *
# ********************
# A. Monotone estimators
#   1. est_mono -- Classical monotone estimator
#   2. ipw_mono_est -- IPW estimator with estimated weights
#   3. mono_est_weights -- Calibration estimator
# B. Nonmonotone estimators
#   1. est_nonmono -- The proposed estimator in the nonmonotone setting
#   2. ipw_cc_oracle -- IPW estimator with oracle weights
#   3. ipw_cc_est -- IPW estimator with weights estimated
#   4. twophase_reg -- Unfinished
#   5. nonmono_est_weights -- Calibration estimator

# ------------------------------------------------------------------------------

# **********************
# * Required Libraries *
# **********************

library(dplyr)
library(CVXR)

# ************************
# * Additional Functions *
# ************************

#' Thanks Hengfang Wang for this matrix inversion function!
hw_inv = function(X, eps = 1e-12) {
  eig.X = eigen(X, symmetric = TRUE)
  P = eig.X[[2]] 
  lambda = eig.X[[1]] 
  ind = lambda > eps
  lambda[ind] = 1 / lambda[ind] 
  lambda[!ind] = 0
  ans = P %*% diag(lambda) %*% t(P)
  return(ans)
}


# ------------------------------------------------------------------------------

# ***********************************
# * Estimating Monotone Missingness *
# ***********************************

est_mono <- function(df) {

  # g(x, y1, y2) = y2
  gfun <- "y2"

  # Algorithm Steps:
  # 1. Estimate pi_11, pi_1+
  # 2. Estimate b1 and b2
  # 3. Estimate \hat \theta

  df_1 <- filter(df, r1 == 1)  # Contains data from first stage
  df_11 <- filter(df, r2 == 1) # Contains data from first and second stage

  # 1. Estimate pi_11, pi_1+
  pi_11_mod <- glm(r2 ~ x + y1, data = df_1, family = binomial(link = logit))
  pi_11 <- pi_11_mod$fitted.values # NOTE: This is NOT of length n.

  df$pi_11 <- 0
  df$pi_11[df$r1 == 1] <- pi_11

  pi_1_mod <- glm(r1 ~ x, data = df, family = binomial(link = logit))
  pi_1 <- pi_1_mod$fitted.values

  # 2. Estimate b1 and b2
  b1_mod <- lm(as.formula(paste0(gfun, " ~ x")), data = df_11)
  b1 <- predict(b1_mod, newdata = df)
  b2_mod <- lm(as.formula(paste0(gfun, " ~ x + y1")), data = df_11)
  b2 <- predict(b2_mod, newdata = df_1) # NOTE: This is NOT of length n.
  df$b2 <- 0 
  df$b2[df$r1 == 1] <- b2


  # 3. Estimate \hat \theta
  df$r1r2pi_11 <- df$r1 * df$r2 / df$pi_11
  df$r1r2pi_11 <- ifelse(is.na(df$r1r2pi_11), 0, df$r1r2pi_11)
  theta <- 
    mean(df$r1r2pi_11 * df[[gfun]]) -
    mean((df$r1 / pi_1 - 1) * b1) - 
    mean((df$r1r2pi_11 - df$r1 / pi_1) * df$b2)

}

ipw_mono_est <- function(df, est = FALSE) {

  if (est) {
    pi_1_mod <- glm(r1 ~ x, data = df, family = binomial(link = logit))
    pi_1 <- pi_1_mod$fitted.values

    pi_11_mod <- 
      glm(r2~x+y1, data = filter(df, r1 == 1), family = binomial(link = logit))
    pi_11 <- numeric(nrow(df))
    pi_11[df$r1 == 1] <- pi_11_mod$fitted.values
    pi_11[df$r1 != 1] <- Inf

  } else {

    pi_1 <- expit(df$x)
    pi_11 <- expit(df$y1)

  }

  mean(df$y2 * df$r2 / (pi_1 * pi_11))

}

mono_est_weights <- function(df) {

  # Steps:
  # 1. Calibrate with w1 weights on Y_1
  # 2. Calibrate with w2 weights on Y_2 using w1 weights on Y_1.

  # 1. Calibrate with w1 weights on Y_1
  n1 <- sum(df$r1)
  w1 <- Variable(n1)
  objective <- Minimize(sum_squares(w1))
  Tx <- mean(df$x)
  xr1_vec <-
    filter(df, r1 == 1) %>%
    pull(x)
  constraints <- 
    list(t(w1) %*% matrix(xr1_vec, ncol = 1) == Tx,
         w1 >= 0,
         t(w1) %*% matrix(rep(1, n1), ncol = 1) == 1)
  p1 <- Problem(objective, constraints)
  res1 <- solve(p1)

  df$w1 <- 0
  df$w1[df$r1 == 1] <- res1$getValue(w1)

  # 2. Calibrate with w2 weights on Y_2 using w1 weights on Y_1.
  n2 <- sum(df$r2) # Only works because of monotone.
  w2 <- Variable(n2)
  objective <- Minimize(sum_squares(w2))
  Exw <- sum(df$w1 * df$x)
  Ey1w <- sum(df$w1 * df$y1)
  xr2_vec <- filter(df, r2 == 1) %>% pull(x)
  y1r2_vec <- filter(df, r2 == 1) %>% pull(y1)
  constraints <-
    list(w2 >= 0,
         t(w2) %*% matrix(rep(1, n2), ncol = 1) == 1,
         t(w2) %*% matrix(xr2_vec, ncol = 1) == Exw,
         t(w2) %*% matrix(y1r2_vec, ncol = 1) == Ey1w)
  p2 <- Problem(objective, constraints)
  res2 <- solve(p2)

  df$w2 <- 0
  df$w2[df$r2 == 1] <- res2$getValue(w2)

  # Objective is E[g] = E[y_2]
  sum(df$w2 * df$y2)

}

# **************************************
# * Estimating Nonmonotone Missingness *
# **************************************

est_nonmono <- function(df, oracle = FALSE, alpha = NA) {

  # g(x, y1, y2) = y2
  gfun <- "y2"

  # Algorithm Steps:
  # 1. Estimate pi_11, pi_1+, pi_2+
  # 2. Estimate b1 and b2, a2
  # 3. Estimate \hat \theta

  df_10 <- filter(df, r1 == 1)          # Contains data from first stage
  df_11 <- filter(df, r1 == 1, r2 == 1) # Contains data from both stages
  df_01 <- filter(df, r2 == 1)          # Contains data from first stage

  # 1. Estimate pi_11, pi_1+, pi_2+
  if (oracle) {
    pi_1 <- df$p1 + df$p2 * df$p21
    pi_2 <- df$p2 + df$p1 * df$p12
    pi_11 <- df$p12 * df$p1 + df$p21 * df$p2
    df$pi_11 <- pi_11

  } else {
    pi_11_mod <- glm(r2 ~ x + y1, data = df_10, family = binomial(link = logit))
    pi_12 <- pi_11_mod$fitted.values # NOTE: This is NOT of length n.

    df$pi_11 <- 0
    df$pi_11[df$r1 == 1] <- pi_12

    if (!is.na(alpha)) {
      pi_21_m <- glm(r1 ~ x + y2, data = df_01, family = binomial(link = logit))
      pi_21 <- pi_21_m$fitted.values
      tmp12_vec <- rep(0, nrow(df))
      tmp21_vec <- rep(0, nrow(df))
      tmp12_vec[df$r1 == 1] <- pi_12
      tmp21_vec[df$r2 == 1] <- pi_21
      df$pi_11 <- (1 - alpha) * tmp12_vec + (alpha) * tmp21_vec
    }

    pi_1_mod <- glm(r1 ~ x, data = df, family = binomial(link = logit))
    pi_1 <- pi_1_mod$fitted.values

    pi_2_mod <- glm(r2 ~ x, data = df, family = binomial(link = logit))
    pi_2 <- pi_2_mod$fitted.values
  }

  # 2. Estimate b1 and b2
  b1_mod <- lm(as.formula(paste0(gfun, " ~ x")), data = df_11)
  b1 <- predict(b1_mod, newdata = df)
  b2_mod <- lm(as.formula(paste0(gfun, " ~ x + y1")), data = df_11)
  b2 <- predict(b2_mod, newdata = df_10) # NOTE: This is NOT of length n.
  df$b2 <- 0 
  df$b2[df$r1 == 1] <- b2

  # a2_mod <- lm(as.formula(paste0(gfun, " ~ x + y2")), data = df_11)
  # a2 <- predict(a2_mod, newdata = df_01) # NOTE: This is NOT of length n.
  df$a2 <- df$y2
  df$a2[df$r2 == 0] <- 0

  # 3. Estimate \hat \theta
  df$r1r2pi_11 <- df$r1 * df$r2 / df$pi_11
  df$r1r2pi_11 <- ifelse(is.na(df$r1r2pi_11), 0, df$r1r2pi_11)

  theta <-
    mean(b1) +
    mean(df$r1 / pi_1 * (df$b2 - b1)) +
    mean(df$r2 / pi_2 * (df$a2 - b1)) +
    mean(df$r1r2pi_11 * (df[[gfun]] - df$b2 - df$a2 + b1))

}

# **************************
# * Alternative Algorithms *
# **************************

ipw_cc_oracle <- function(df) {

  gfun <- "y2"

  tmp <- df %>%
    filter(r1 == 1, r2 == 1) %>%
    mutate(pi_11 = p12 * p1 + p21 * p2) 

  sum(tmp[[gfun]] / tmp$pi_11) / nrow(df)

}

ipw_cc_est <- function(df, oracle = FALSE, alpha = NA) {

  gfun <- "y2"


  df_10 <- filter(df, r1 == 1)          # Contains data from first stage
  df_11 <- filter(df, r1 == 1, r2 == 1) # Contains data from both stages
  df_01 <- filter(df, r2 == 1)          # Contains data from first stage

  # 1. Estimate pi_11, pi_1+, pi_2+
  if (oracle) {
    pi_1 <- df$p1 + df$p2 * df$p21
    pi_2 <- df$p2 + df$p1 * df$p12
    pi_11 <- df$p12 * df$p1 + df$p21 * df$p2
    df$pi_11 <- pi_11

  } else {
    pi_11_mod <- glm(r2 ~ x + y1, data = df_10, family = binomial(link = logit))
    pi_12 <- pi_11_mod$fitted.values # NOTE: This is NOT of length n.

    df$pi_11 <- 0
    df$pi_11[df$r1 == 1] <- pi_12

    if (!is.na(alpha)) {
      pi_21_m <- glm(r1 ~ x + y2, data = df_01, family = binomial(link = logit))
      pi_21 <- pi_21_m$fitted.values
      tmp12_vec <- rep(0, nrow(df))
      tmp21_vec <- rep(0, nrow(df))
      tmp12_vec[df$r1 == 1] <- pi_12
      tmp21_vec[df$r2 == 1] <- pi_21
      df$pi_11 <- (1 - alpha) * tmp12_vec + (alpha) * tmp21_vec
    }
  }

  df$r1r2pi_11 <- df$r1 * df$r2 / df$pi_11
  df$r1r2pi_11 <- ifelse(is.na(df$r1r2pi_11), 0, df$r1r2pi_11)

  mean(df[[gfun]] * df$r1r2pi_11)

}

# png("diffhist.png")
# filter(df, r1 == 1, r2 == 1) %>%
#   dplyr::select(-p0, -p1, -p2, -p12, -p21) %>%
#   mutate(diff = pi_11 - true_pi11) %>%
#   pull(diff) %>% 
#   hist(main = "Histogram of Estimated - True") 
# dev.off()

twophase_reg <- function(df, reg_on_y1 = TRUE) {

  gfun <- "y2"

  # Steps:
  # 1. Get weights
  # 2. Estimate beta in complete case
  # 3. Impute missing values

  # Get weights
  pi_11 <- df$p12 * df$p1 + df$p21 * df$p2
  df$w <- 1 / pi_11

  if (reg_on_y1) {
    df_11 <- filter(df, r1 == 1, r2 == 1)
    mod <- 
      lm(as.formula(paste0(gfun, " ~ x + y1")), weights = w, data = df_11)

    df_1 <- filter(df, r1 == 1)

    return(mean(predict(mod, newdata = df_1)))

  } else {
    df_11 <- filter(df, r1 == 1, r2 == 1)
    mod <- lm(as.formula(paste0(gfun, " ~ x")), weights = w, data = df_11)

    return(mean(predict(mod, newdata = df)))

  }
}

threephase_reg <- function(df) {

  gfun <- "y2"

  # Steps:
  # 1. Get weights
  # 2. Estimate beta in complete case
  # 3. Estimate missing y1 values
  # 4. Impute missing values

  # Get weights
  pi_11 <- df$p12 * df$p1 + df$p21 * df$p2
  df$w <- 1 / pi_11

  # Estimate beta
  df_11 <- filter(df, r1 == 1, r2 == 1)
  mean_y2 <- sum(df_11$w * df_11$y2) / sum(df_11$w)
  mean_x2 <- sum(df_11$w * df_11$x) / sum(df_11$w)
  mean_y1 <- sum(df_11$w * df_11$y1) / sum(df_11$w)

  df_11[[gfun]] <- df_11[[gfun]] - mean_y2
  df_11$x <- df_11$x - mean_x2
  df_11$y1 <- df_11$y1 - mean_y1

  mod <- lm(as.formula(paste0(gfun, " ~ x + y1")), weights = w, data = df_11)

  df_1 <- 
    filter(df, r1 == 1) %>%
    mutate(w_y1 = 1 / p1) 

  mean_x1 <- sum(df_1$x * df_1$w_y1) / sum(df_1$w_y1)
  mean_y11 <- sum(df_1$y1 * df_1$w_y1) / sum(df_1$w_y1)

  df_1$x <- df_1$x - mean_x1
  df_1$y1 <- df_1$y1 - mean_y11

  mod_y1 <- lm(y1 ~ x, weights = w_y1, data = df_1)

  # 1. [X] Change estimate of \hat \beta_3 to be centered at 0.
  # 2. [X] Construct three phase estimator according to (3.3.38).
  #   * Need \bar y_2, \bar x_1, \bar x_2, \bar \hat y_1, \bar y_1
  bar_x1 <- sum(1 / df$p1 * df$x) / sum(1 / df$p1)
  bar_est_y1 <- mean(predict(mod_y1, newdata = df))

  mean_y2 + (bar_x1 - mean_x2) * mod$coefficients[2] + 
    (bar_est_y1 - mean_y1) * mod$coefficients[3]

}

nonmono_est_weights <- function(df) {

  gfun <- "y2"

  # Steps:
  # 1. Compute w1--the weights between U and A1.
  # 2. Compute w2--the weights between U and A2.
  # 3. Compute wc--the weights between A1, A2, U, and the core points.

  # 1. Compute w1--the weights between U and A1.
  #   a. Estimate E[g_i \mid X_i].
  #   b. Compute weights on A1.
  Eg_mod <- lm(as.formula(paste0(gfun, " ~ x")), data = df)
  df$Egx <- predict(Eg_mod)

  Mgx <- mean(df$Egx)
  n1 <- sum(df$r1)
  w1 <- Variable(n1)
  Egx_r1 <- filter(df, r1 == 1) %>% pull(Egx)
  const <- list(w1 >= 0,
                t(w1) %*% matrix(rep(1, n1), ncol = 1) == 1,
                t(w1) %*% matrix(Egx_r1, ncol = 1) == Mgx)

  p1 <- Problem(Minimize(sum_squares(w1)), const)
  res1 <- solve(p1)

  df$w1 <- 0
  df$w1[df$r1 == 1] <- res1$getValue(w1)

  # 2. Compute w2--the weights between U and A2.
  #   a. Compute weights on A2.
  n2 <- sum(df$r2)
  w2 <- Variable(n2)
  Egx_r2 <- filter(df, r2 == 1) %>% pull(Egx)
  const <- list(w2 >= 0,
                t(w2) %*% matrix(rep(1, n2), ncol = 1) == 1,
                t(w2) %*% matrix(Egx_r2, ncol = 1) == Mgx)

  p2 <- Problem(Minimize(sum_squares(w2)), const)
  res2 <- solve(p2)

  df$w2 <- 0
  df$w2[df$r2 == 1] <- res2$getValue(w2)

  # 3. Compute wc--the weights between A1, A2, U, and the core points.
  #   a. Estimate E[g_i \mid X_i, Y_{1i}].
  #   b. Estimate E[g_i \mid X_i, Y_{2i}].
  #   c. Compute core weights on each part simultaneously.
  df_1 <- filter(df, r1 == 1)
  Eg1_mod <- lm(as.formula(paste0(gfun, " ~ x + y1")), data = df_1)
  df$Eg1 <- 0
  df$Eg1[df$r1 == 1] <- predict(Eg1_mod)

  # E[y_2 | y_2 ] = y_2
  # df_2 <- filter(df, r2 == 1)
  # Eg2_mod <- lm(as.formula(paste0(gfun, " ~ x + y2")), data = df_2)

  df_11 <- filter(df, r1 == 1, r2 == 1)
  n_11 <- nrow(df_11)
  wc <- Variable(n_11)

  G1 <- sum(df$w1 * df$Eg1)
  G2 <- sum(df$w2 * df$y2)
  G3 <- mean(df$Egx)
  const <- list(wc >= 0,
                t(wc) %*% matrix(rep(1, n_11), ncol = 1) == 1,
                t(wc) %*% matrix(df_11$Eg1, ncol = 1) == G1,
                t(wc) %*% matrix(df_11$y2, ncol = 1) == G2,
                t(wc) %*% matrix(df_11$Egx, ncol = 1) == G3)

  p3 <- Problem(Minimize(sum_squares(wc)), const)
  res3 <- solve(p3)

  sum(df[[gfun]][df$r1 == 1 & df$r2 == 1] * res3$getValue(wc))

}

f3p_cali <- function(df) {

  # Steps
  # 1. Construct vectors / matrices.
  # 2. Run calibration

  # Do I need different steps if the data is monotone vs nonmonotone? No.
  pi_11 <- df$p12 * df$p1 + df$p21 * df$p2
  df$w <- 1 / pi_11
  df_11 <- filter(df, r1 == 1, r2 == 1)
  df_1 <- filter(df, r1 == 1)
  df_1$w_1 <- 1 / df_1$p1

  # tot_y2 <- sum(df_11$w * df_11$y2)
  tot_x <- sum(df$x)
  tot_z <- sum(df_1$w_1 * df_1$y1) 

  n2 <- nrow(df_11)

  # 1. Construct vectors / matrices.
  w_vec <- matrix(df_11$w, ncol = 1)
  L_mat <- base::diag(df_11$w)
  L_inv <- base::diag(1 / df_11$w)
  x_vec <- df_11$x
  z_vec <- df_11$y1
  c_vec <- Variable(n2)

  # 2. Run calibration
  const <- 
    list(c_vec >= 0,
         t(c_vec) %*% matrix(rep(1, n2), ncol = 1) == n2, # We estimate totals
         t(c_vec) %*% matrix(df_11$x, ncol = 1) == tot_x,
         t(c_vec) %*% matrix(df_11$y1, ncol = 1) == tot_z)

  prob <- Problem(Minimize(quad_form(c_vec - w_vec, L_inv)), const)
  res <- solve(prob)

  c_opt <- res$getValue(c_vec)

  return(mean(c_opt * df_11$y2))

}
