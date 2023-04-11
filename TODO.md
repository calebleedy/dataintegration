
# TODO

## April 11, 2023

* [ ] Write up doubly robust estimator
* [ ] Write up estimator with $\Pr(R_1 | X, Y_2)$
* [ ] Read data integration papers:
  * [ ] Yang and Kim (2020)
  * [ ] Kim et al. (2020) (JRSSA)
  * [ ] Chen, Yang, Kim (2022)

## April 4, 2023

* [X] Nonmonotone: make $Y_1, Y_2$ conditional on $X$.
* [X] Nonmonotone: use $R_1, R_2 | X$. The goal is to decide $a_2$ and
      $b_2$.
* [ ] Try using both $\theta_k = E[Y_k]$.
* [ ] Try adding two-phase regression estimator.
* [X] Test double robustness of proposed estimator.

## April 3, 2023

* [o] Comparisons for nonmonotone simulation:
  * [X] IPW with complete cases (oracle probabilites)
  * [ ] IPW with complete cases (estimated probabilites)
* [X] Write up progress
* [X] Create simulation where I expect proposal to fail because of how it
      estimates $\pi_{11}$.

## March 29, 2023

* [X] Simulation study for nonmonotone missingness:
  * [X] Generate monotone simulations
  * [X] Generate nonmonotone simulations
  * [X] Estimate monotone simulation
  * [X] Estimate nonmonotone simulation
