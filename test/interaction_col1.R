# nolint start
#This file shows that calculate_omega_overlap may have some problems
#when one of the interacrion matrix has an colunm proportional to rep(1,num).
library(geometry)
library(uniformly)
library(dplyr)
library(tidyverse)
library(mvtnorm)

library(feasoverlap)
# library(pracma)
# source("R/overlap.R")

# 2d case
set.seed(22)
A <- matrix(c(-1,0.375,-1,-1),2,2)
B <- matrix(c(-1,0,0,-1),2,2)
calculate_omega_overlap(A,B)
#> Error in UseMethod("determinant"

# 3d case
set.seed(32)
AA <- interaction_matrix_random(3,1,0.8)
BB <- matrix(c(0.375,-0.18,0.45,0.42,0.6,-2,-1,-1,-1),3,3)
calculate_omega_overlap(AA,BB)
calculate_omega_overlap(BB,AA)
#> Error in solve.default(coeff_matrix, coeff_vector)
# nolint end

