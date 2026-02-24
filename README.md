
<!-- README.md is generated from README.Rmd. Please edit that file -->

# feasoverlap <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/clsong/feasoverlap/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/clsong/feasoverlap/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**feasoverlap** computes the overlap between high-dimensional feasibility domains — a key quantity for understanding species coexistence in ecological communities.

## 📦 Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("clsong/feasoverlap")
```

## 🚀 Quick Start

### Overlap between two interaction matrices

``` r
library(feasoverlap)

# Generate two random interaction matrices
set.seed(1)
A <- interaction_matrix_random(num = 4, stren = 0.4, conne = 1)
set.seed(2)
B <- interaction_matrix_random(num = 4, stren = 0.4, conne = 1)

calculate_omega(A) # feasibility domain size of A
#> [1] 0.4895205
calculate_omega(B) # feasibility domain size of B
#> [1] 0.4743946

calculate_omega_overlap(A, B) # overlap of A and B
#> [1] 0.316648
calculate_omega_overlap(B, A) # overlap of B and A
#> [1] 0.326564
```

### Feasibility under biological constraints

``` r
library(feasoverlap)

set.seed(4)
A <- interaction_matrix_random(3, 0.4, 1)

# Constraint: all growth rates must be positive
C1 <- diag(c(-1, -1, -1), 3)

# Constraint: species 1 & 3 negative, species 2 positive growth rate
C2 <- diag(c(1, -1, 1), 3)

calculate_omega(A)              # unconstrained feasibility
#> [1] 0.5073363
calculate_omega_overlap(A, C1)  # feasibility under C1
#> [1] 0.4306551
calculate_omega_overlap(A, C2)  # feasibility under C2
#> [1] 0.16287
```

## 📖 Key Concepts

| Function | Description |
|---|---|
| `interaction_matrix_random()` | Generate random interaction matrices |
| `calculate_omega()` | Compute the normalized size of a feasibility domain |
| `calculate_omega_overlap()` | Compute the overlap between two feasibility domains |

## 📄 Citation

If you use **feasoverlap** in your research, please cite:

> Song, C. & Long, C. (2025). feasoverlap: Measure the Overlap Between Two High-Dimensional Feasibility Domains. R package.

## 📝 License

MIT © [Chuliang Song](https://chuliangsong.com), [Chengyi Long](https://github.com/chengyilong)
