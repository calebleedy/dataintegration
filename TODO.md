
# TODO

## April 20, 2023

* [ ] Add papers from `references.bib` to the `Reference/` directory
* [ ] Go through proof in Chen, Yang, Kim (2022)
* [ ] Simulation study of Kim et al. (2020) 
* [ ] Proof of Kim et al. (2020)

## April 18, 2023

* [X] Compute bias of $\hat \theta_{eff}$. (We want to understand double
      robustness of estimator.) It is unbiased if either model is correct.
* [X] Fix $a_2$ (It needs to be $E[Y_2 \mid X]$)
* [ ] What is the best choice of $a_2$ and $b_2$? (Find optimal.)

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
