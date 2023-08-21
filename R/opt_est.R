# Title: opt_est.R
# Created by: Caleb Leedy
# Created on: August 19, 2023
# Purpose: This file contains the functions needed for computation of the 
# optimal estimator in the non-monotone case.

# ************
# * Overview *
# ************
# The goal of this file is to create code to estimate the following:
# \hat \theta = n^{-1} \sum_{i = 1}^n E[g_i | X_i] +
#   n^{-1} \sum_{i = 1}^n (R_{1i} / \pi_{1+}) (E[g_i | X_i, Y_{1i}] - E[g_i | X_i]) +
#   n^{-1} \sum_{i = 1}^n (R_{2i} / \pi_{2+}) (E[g_i | X_i, Y_{2i}] - E[g_i | X_i]) +
#   n^{-1} \sum_{i = 1}^n (R_{1i}R_{2i} / \pi_{11}) (E[g_i | X_i, Y_{1i}, Y_{2i}] -
#     E[g_i | X_i, Y_{1i}] - E[g_i | X_i, Y_{2i}] + E[g_i | X_i])
