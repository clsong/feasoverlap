# Producing nearly non-singular matrix
rm(list = ls())

library(feasoverlap)

num <- 3; stren <- 2; conne <- 1

interaction_matrix_ill <- function(num, stren, conne, epsilon = 10^(-5), threshold = 0.001) {
  inte <- interaction_matrix_random(num, stren, conne)
  new_col <- floor(num/2) + 1
  inte[new_col,1] <- (max(abs(inte[new_col,1]), threshold)) #set a threshold to avoid zero entry
  #inte[new_col,1] <- (max(inte[new_col,1], threshold))
  factor <- (-1)/inte[new_col,1]
  inte[,new_col] <- factor*inte[,1] + epsilon * rnorm(num)
  inte[new_col,new_col] <- -1
  return(inte)
}

set.seed(50)
A <- interaction_matrix_ill(num, stren, conne)

B <- A[1:(num-1),1:(num-1)] %>%
  cbind(c(rep(0,num-1))) %>%
  rbind(c(rep(0,num-1),-1))

calculate_omega_overlap(A,B)
min(calculate_omega(A), calculate_omega(B))

(ratio <- calculate_omega_overlap(A,B)/min(calculate_omega(A), calculate_omega(B)))
#ratio of raw-omega-value: seems to be scaled with epsilon
ratio^num
