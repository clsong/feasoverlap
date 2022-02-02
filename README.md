
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
#> [1] 0.4710153
calculate_omega_overlap(A, A) # overlap between the same matrix
#> [1] 0.4710153

calculate_omega(B) # relative size of interaction matrix
#> [1] 0.461008
calculate_omega_overlap(B, B) # overlap between the same matrix
#> [1] 0.461008

calculate_omega_overlap(A, B) # overlap of two interaction matrices
#> [1] 0.3115771
calculate_omega_overlap(B, A) # overlap of two interaction matrices
#> [1] 0.3134456
```

## Example of the normalized size of the feasibility domain of a random interaction matrix under linear biological inequalities

``` r
library(feasoverlap)

set.seed(4)
A <- interaction_matrix_random(3, 0.4, 1) #generate a random interaction matrix
C1 <- diag(c(-1,-1,-1), 3) #imposing a biological constraint. Here it refers to that the growth rates of all species have to be positive
C2 <- diag(c(1,-1,1), 3)  #imposing a biological constraint. Here it refers to that the growth rates of species 1 and 3 have to be negative, and the growth rates of species 2 has to be positive

calculate_omega(A) #relative size of the original interaction matrix
#> [1] 0.4954943
calculate_omega_overlap(A, C1) #the normalized size of the feasibility domain of a random interaction matrix under linear biological constriants C1
#> [1] 0.4222341
calculate_omega_overlap(A, C2) #the normalized size of the feasibility domain of a random interaction matrix under linear biological constriants C2
#> [1] 0.1584054
```
