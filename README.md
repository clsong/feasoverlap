
<!-- README.md is generated from README.Rmd. Please edit that file -->

# feasoverlap

<!-- badges: start -->
<!-- badges: end -->

The goal of feasoverlap is to compute the overlap between two
feasibility domains.

## Installation

You can install the development version of feasoverlap from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("clsong/feasoverlap")
```

## Example of the overlap of two random interaction matrix

``` r
library(feasoverlap)

# generate two random interaction matrices
set.seed(1)
A <- interaction_matrix_random(num = 4, stren = 0.4, conne = 1)
set.seed(2)
B <- interaction_matrix_random(num = 4, stren = 0.4, conne = 1)

calculate_omega(A) # relative size of interaction matrix
#> [1] 0.6696702
calculate_omega_overlap(A, A) # overlap between the same matrix
#> [1] 0.6696702

calculate_omega(B) # relative size of interaction matrix
#> [1] 0.4679525
calculate_omega_overlap(B, B) # overlap between the same matrix
#> [1] 0.4679525

calculate_omega_overlap(A, B) # overlap of two interaction matrices
#> [1] 0.6130666
calculate_omega_overlap(B, A) # overlap of two interaction matrices
#> [1] 0.6120678
```
