library(feasoverlap)
library(tidyverse)

load("To_Chuliang.RData")

# original matrix
A1 <- A_comp[[1]]
# matrix without basal competition
A2 <- A1
A2[1:4, 1:4] <- 0
diag(A2)[1:4] <- -1


# Similar omega
calculate_omega(A1)
calculate_omega(A2)

# Very different omega overlap (a bug in the package)
calculate_omega_overlap(A1, B, nsamples = 100)
calculate_omega_overlap(A2, B, nsamples = 100)

# new function
calculate_omega_constraint(A1, B)
calculate_omega_constraint(A2, B)
