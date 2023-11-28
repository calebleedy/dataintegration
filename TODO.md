
# TODO

* [ ] Check to see if F3P is the calibration estimator as conjectured on pg. 35
      of notes.
* [ ] Understand the differences between estimators. Why do some work when
      others do not?
* [ ] Go through proof in Chen, Yang, Kim (2022)
* [ ] Look at differences between different three-phase estimators
* [ ] Recap Kim and Rao (2012, Biometrika). We could extend it to three phases
      or two step calibration.
* [ ] Recap Chen and Kim (2014)
* [ ] Show that proposed estimator is optimal

# Done

## November 9, 2023

* [X] Compare results to normal model
* [X] Unpool estimates
* [X] Run simulation study
* [X] Literature review of combining samples
* [X] Write implicit assumptions in simulation
## September 11, 2023
* [X] Non-linear simulation
* [X] Use independent samples
* [X] Email Dr. Kim results

## September 6, 2023
* [X] Optimality simulation. The idea is that we want to see if the proposed
estimator can recover the optimal estimator of the simulation in the case where 
everything is normal. The compute the optimal estimator we can use the EM 
algorithm.
## August 29, 2023

* [X] Fix notation
* [X] Compare coefficients

## August 19, 2023
* [X] Fix notation ($\sigma_{yy} \to \rho$) and ($E[\varepsilon] = 0$)
## May 4, 2023 

* [X] Proof in Kim et al. (2020)

## May 2, 2023

* [X] Recap Lohr and Raghunathan (2017)

## April 25, 2023

* [X] Fix simulation explanation
* [X] Efficiency comparison with existing method (compare with two-phase
      sampling) of sampling problem
* [X] Compare with calibration
* [X] Write up calibration idea

## April 20, 2023

* [X] Add papers from `references.bib` to the `Reference/` directory
* [X] Simulation study of Kim et al. (2020) 

## April 18, 2023

* [X] Compute bias of $\hat \theta_{eff}$. (We want to understand double
      robustness of estimator.) It is unbiased if either model is correct.
* [X] Fix $a_2$ (It needs to be $E[Y_2 \mid X]$)
* [X] What is the best choice of $a_2$ and $b_2$? (Find optimal.)

## April 12, 2023

* [X] Simulation study with modified estimator.

## April 11, 2023

* [X] Write up doubly robust estimator
* [X] Write up estimator with $\Pr(R_1 \mid X, Y_2)$
* [X] Read data integration papers:
  * [X] Yang and Kim (2020)
  * [X] Kim et al. (2020) (JRSSA)
  * [X] Chen, Yang, Kim (2022)

## April 4, 2023

* [X] Nonmonotone: make $Y_1, Y_2$ conditional on $X$.
* [X] Nonmonotone: use $R_1, R_2 | X$. The goal is to decide $a_2$ and
      $b_2$.
~~* [ ] Try using both $\theta_k = E[Y_k]$.~~
~~* [ ] Try adding two-phase regression estimator.~~
* [X] Test double robustness of proposed estimator.

## April 3, 2023

* [X] Comparisons for nonmonotone simulation:
  * [X] IPW with complete cases (oracle probabilites)
  * [X] IPW with complete cases (estimated probabilites)
* [X] Write up progress
* [X] Create simulation where I expect proposal to fail because of how it
      estimates $\pi_{11}$.

## March 29, 2023

* [X] Simulation study for nonmonotone missingness:
  * [X] Generate monotone simulations
  * [X] Generate nonmonotone simulations
  * [X] Estimate monotone simulation
  * [X] Estimate nonmonotone simulation
