#!/bin/bash

# File: sim_script.sh
# Created by: Caleb Leedy
# Created on: November 12, 2024
# Purpose: This script runs all of the code for the simulations for the project:
# Debiased Calibration for Generalized Two-Phase Sampling.

Rscript final_sims.R "twophasesamp" 0
Rscript final_sims.R "twophasesamp" 1
Rscript final_sims.R "nonnestedtp" 0
Rscript final_sims.R "nonnestedtp" 1
Rscript final_sims.R "multisource" 0
Rscript final_sims.R "multisource" 1
Rscript final_sims.R "multisourcealpha" 0
Rscript final_sims.R "multisourcealpha" 1
