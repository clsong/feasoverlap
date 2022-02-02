library(tidyverse)
library(mvtnorm)
library(mgcv)
library(geometry)


num = 5
stren = 0.4
conne = 1

#generate two random interaction matrices
# set.seed(10)
A <- interaction_matrix_random(num, stren, conne)
# set.seed(2)
B <- interaction_matrix_random(num, stren, conne)

calculate_omega_overlap(A, B)
calculate_omega(A)
calculate_omega(B)
calculate_omega_overlap(B,B)


set.seed(10)
A <- interaction_matrix_random(3, 0.4, 1) #generate a random interaction matrix
C1 <- diag(c(-1,-1,-1), 3) #imposing a biological constraint. Here it refers to that the growth rates of all species have to be positive
C2 <- diag(c(1,-1,1), 3)  #imposing a biological constraint. Here it refers to that the growth rates of species 1 and 3 have to be negative, and the growth rates of species 2 has to be positive

calculate_omega(A) #relative size of the original interaction matrix
calculate_omega_overlap(A, C1) #the size of the feasibility domain of a random interaction matrix under linear biological constriants I1
calculate_omega_overlap(A, C2) #the size of the feasibility domain of a random interaction matrix under linear biological constriants I2
