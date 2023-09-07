
# Data Integration

The goal of this project is to be able to merge multiple (non-probability) data 
sets with a probability sample and still have valid estimates for the joint
distributions of the variables. The project focuses on when the missingness is 
non-monotone. Ultimately, this research should help NRI estimation. 
There are several steps to understanding an optimal estimator:

1. Merging two samples with non-monotone missingness where the inclusion 
probability ($\pi$) is known and neither a function of any of the missing 
variables $Y_k$.
2. Extending Step 1 to multiple samples with non-monotone missingness.
3. Merging one probability sample with a non-probability sample in which the 
missingness is non-monotone.
4. Extending Step 3 to mulitple non-probability samples.

# Related Literature

* Robins and Gill (1997) [[1]](#1) uses a non-monotone framework to merge data that 
is conditionally missing at random.
* Merkouris (2004) [[2]](#2) combines multiple data sets with non-monotone 
missingness under a normal model.
* Kim et al. (2021) [[3]](#3) and Chen, Yang, and Kim (2022) [[4]](#4) focus on
how to combine a probability sample with a non-probability sample.


## Tasks

For specific items see [TODO.md](TODO.md). Overall, we explore(d) the 
following:

1. Implementing and comparing different estimators under the Robins and Gill 
(1997) [[1]](#1) framework.
2. Finding the optimal estimator for non-monotone missing data.

## References

<a id="1">[1]</a>
Robins, J. M., & Gill, R. D. (1997). Non‐response models for the analysis of 
non‐monotone ignorable missing data. Statistics in medicine, 16(1), 39-56.

<a id="2">[2]</a>
Merkouris, Takis (2004). Combining independent regression estimators from
multiple surveys. Journal of the American Statistical Association, 99(468),
1131-1139.

<a id="3">[3]</a>
Kim, J. K., Park, S., Chen, Y., & Wu, C. (2021). Combining non-probability and
probability survey samples through mass imputation. Journal of the Royal
Statistical Society Series A: Statistics in Society, 184(3), 941-963.

<a id="4">[4]</a>
Chen, S., Yang, S., & Kim, J. K. (2022). Nonparametric mass imputation for data
integration. Journal of survey statistics and methodology, 10(1), 1-24.
