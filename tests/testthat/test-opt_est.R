# Title: test-opt_est.R
# Created by: Caleb Leedy
# Created on: August 21, 2023
# Purpose: This file contains the code to check to see if all of the functions
# in R/optimal_simulation.R work.

# *********
# * Tests *
# *********

# Test expit
# expit has been moved to R/utilities.R
test_that("expit works", {
  expect_equal(expit(0), 0.5)
  expect_equal(expit(Inf), 1)
  expect_equal(expit(-Inf), 0)
})

# Test expect_g
df <- 
  tibble(X = rep(c(1, 2), 5)) |>
  mutate(Y1 = 1 + 2 * X) |>
  mutate(Y2 = 2 + 3 * X) |>
  mutate(delta_1 = c(1, 1, 1, 0, 0, 0, 1, 1, 1, 0)) |>
  mutate(delta_2 = c(1, 1, 1, 0, 1, 1, 0, 0, 1, 0)) |>
  mutate(prob_11 = 0.4, prob_1 = 0.6, prob_2 = 0.6)

df_1 <- filter(df, delta_1 == 1)
df_2 <- filter(df, delta_2 == 1)
tmp1 <- rep(0, nrow(df))
tmp1[df$delta_1 == 1] <- df$Y2[df$delta_1 == 1]
tmp2 <- rep(0, nrow(df))
tmp2[df$delta_2 == 1] <- df$Y2[df$delta_2 == 1]

test_that("expect_g works", {
  expect_equal(expect_g(df, "Y2", c("X"), df, rep(1, nrow(df))), df$Y2)
  expect_equal(expect_g(df, "Y2", c("X", "Y1"), df_1, df$delta_1), tmp1)
  expect_equal(expect_g(df, "Y2", c("X", "Y2"), df_2, df$delta_2), tmp2)
})

# Test get_v_powers
test_that("get_v_powers works", {
  expect_equal(get_v_powers(c("X"), 1), c("I(X)"))
  expect_equal(get_v_powers(c("X", "Y1"), 1), c("I(X)", "I(Y1)"))
})

# Test get_cond_exp
test_that("get_cond_exp works with cols = 1", {
  df <- gen_optsim_data(n = 1000, theta = 0)
  obs_y2 <- matrix(df$Y2, ncol = 1)
  mean_y1 <- filter(df, delta_1 == 1) |> pull(Y1) |> mean()
  mean_y2 <- filter(df, delta_2 == 1) |> pull(Y2) |> mean()
  sigma_mat <- list(list(2, 1), list(1, 2))
  cond_exp <- get_cond_exp(obs_y2, mean_y1, mean_y2, sigma_mat)
  ans <- mean_y1 + 0.5 * (obs_y2 - mean_y2)

  expect_equal(class(cond_exp), class(ans))
  expect_equal(dim(cond_exp), dim(ans))
  expect_lt(sum(abs(as.numeric(cond_exp) - as.numeric(ans))), 1e-4)

  df <- gen_optsim_data(n = 1000, theta = 5)
  obs_y2 <- matrix(df$Y2, ncol = 1)
  mean_y1 <- filter(df, delta_1 == 1) |> pull(Y1) |> mean()
  mean_y2 <- filter(df, delta_2 == 1) |> pull(Y2) |> mean()
  cond_exp <- get_cond_exp(obs_y2, mean_y1, mean_y2, sigma_mat)
  ans <- mean_y1 + 0.5 * (obs_y2 - mean_y2)

  expect_equal(class(cond_exp), class(ans))
  expect_equal(dim(cond_exp), dim(ans))
  expect_lt(sum(abs(as.numeric(cond_exp) - as.numeric(ans))), 1e-4)

  df <- gen_optsim_data(n = 1000, theta = -5)
  obs_y2 <- matrix(df$Y2, ncol = 1)
  mean_y1 <- filter(df, delta_1 == 1) |> pull(Y1) |> mean()
  mean_y2 <- filter(df, delta_2 == 1) |> pull(Y2) |> mean()
  cond_exp <- get_cond_exp(obs_y2, mean_y1, mean_y2, sigma_mat)
  ans <- mean_y1 + 0.5 * (obs_y2 - mean_y2)

  expect_equal(class(cond_exp), class(ans))
  expect_equal(dim(cond_exp), dim(ans))
  expect_lt(sum(abs(as.numeric(cond_exp) - as.numeric(ans))), 1e-4)
})

test_that("get_cond_exp works with cols = 2", {

  df <- gen_optsim_data(n = 1000, theta = 0)
  obs_y2 <- filter(df, delta_2 == 1) %>% dplyr::select(X, Y2) %>% as.matrix()
  mean_y1 <- filter(df, delta_1 == 1) |> pull(Y1) |> mean()
  mean_y2 <- filter(df, delta_2 == 1) |> pull(Y2) |> mean()
  mean_x <- mean(df$X)
  y2_mat <- matrix(c(mean_x, mean_y2), nrow = nrow(obs_y2), ncol = 2, byrow = TRUE)

  sigma_mat <- list(list(2, 
                         matrix(c(1, 1), nrow = 1)),
                    list(matrix(c(1, 1), nrow = 2),
                         matrix(c(1, 1, 1, 2), nrow = 2)))
  cond_exp <- get_cond_exp(obs_y2, mean_y1, y2_mat, sigma_mat)
  ans <- mean_y1 + (matrix(obs_y2[, 1], ncol = 1) - mean_x)

  expect_equal(class(cond_exp), class(ans))
  expect_equal(dim(cond_exp), dim(ans))
  expect_lt(sum(abs(as.numeric(cond_exp) - as.numeric(ans))), 1e-4)

})
