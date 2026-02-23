test_that("calculate_omega returns expected range", {
  set.seed(42)
  A <- interaction_matrix_random(num = 3, stren = 0.4, conne = 1)
  omega <- calculate_omega(A)
  expect_true(omega >= 0 && omega <= 1)
})

test_that("calculate_omega of identity (diagonal) matrix is positive", {
  omega_diag <- calculate_omega(diag(3))
  expect_true(omega_diag > 0) # Should be positive
  # For identity matrix in n=3: Sigma=I, P(X>0) = 1/2^3 = 1/8
  # Normalized: (1/8)^(1/3) = 0.5
  omega_sphere <- calculate_omega(diag(3), method = "sphere")
  expect_equal(omega_sphere, 0.5, tolerance = 0.05)
})

test_that("calculate_omega_overlap(A, A) ≈ calculate_omega(A) [sphere]", {
  set.seed(1)
  A <- interaction_matrix_random(num = 3, stren = 0.4, conne = 1)
  omega_single <- calculate_omega(A, method = "sphere")
  omega_self <- calculate_omega_overlap(A, A, method = "sphere")
  expect_equal(omega_single, omega_self, tolerance = 0.02)
})

test_that("calculate_omega_overlap is EXACTLY symmetric [sphere]", {
  set.seed(1)
  A <- interaction_matrix_random(num = 3, stren = 0.4, conne = 1)
  set.seed(2)
  B <- interaction_matrix_random(num = 3, stren = 0.4, conne = 1)
  omega_AB <- calculate_omega_overlap(A, B, method = "sphere")
  omega_BA <- calculate_omega_overlap(B, A, method = "sphere")
  # Sphere method is symmetric by construction (pmvnorm has slight internal MC)
  expect_equal(omega_AB, omega_BA, tolerance = 1e-3)
})

test_that("calculate_omega_overlap returns zero for non-overlapping cones", {
  A <- matrix(c(-1, 0, 0, -1), 2, 2)
  B <- matrix(c(1, 0, 0, 1), 2, 2)
  omega <- calculate_omega_overlap(A, B, method = "sphere")
  expect_equal(omega, 0, tolerance = 1e-6)
})

test_that("three methods agree on overlap [cross-validation]", {
  set.seed(1)
  A <- interaction_matrix_random(num = 3, stren = 0.4, conne = 1)
  set.seed(2)
  B <- interaction_matrix_random(num = 3, stren = 0.4, conne = 1)

  omega_sphere <- calculate_omega_overlap(A, B, method = "sphere")
  omega_mc <- calculate_omega_overlap(A, B, method = "montecarlo", nsamples = 10000)
  omega_ch <- calculate_omega_overlap(A, B, method = "convex_hull", nsamples = 100)

  # All three should agree within reasonable tolerance
  expect_equal(omega_sphere, omega_mc, tolerance = 0.1)
  expect_equal(omega_sphere, omega_ch, tolerance = 0.15)
})

test_that("sphere method agrees with single-cone sphere for self-overlap [4D]", {
  set.seed(10)
  A <- interaction_matrix_random(num = 4, stren = 0.4, conne = 1)
  omega_single <- calculate_omega(A, method = "sphere")
  omega_self <- calculate_omega_overlap(A, A, method = "sphere")
  expect_equal(omega_single, omega_self, tolerance = 0.02)
})

test_that("montecarlo method is approximately symmetric", {
  set.seed(1)
  A <- interaction_matrix_random(num = 3, stren = 0.4, conne = 1)
  set.seed(2)
  B <- interaction_matrix_random(num = 3, stren = 0.4, conne = 1)
  omega_AB <- calculate_omega_overlap(A, B, method = "montecarlo", nsamples = 5000)
  omega_BA <- calculate_omega_overlap(B, A, method = "montecarlo", nsamples = 5000)
  expect_equal(omega_AB, omega_BA, tolerance = 0.1)
})

test_that("calculate_omega_overlap falls back to MC for singular matrices", {
  # Singular matrix
  A <- matrix(c(-1, 0.5, -2, 1), 2, 2) # det = 0
  B <- diag(-1, 2)
  # Should not crash — falls back to montecarlo
  omega <- calculate_omega_overlap(A, B, method = "sphere")
  expect_true(is.numeric(omega))
})

test_that("inside_vertex_detection handles identical matrices", {
  set.seed(10)
  A <- interaction_matrix_random(num = 3, stren = 0.4, conne = 1)
  inside <- inside_vertex_detection(A, A)
  expect_true(length(inside) >= ncol(A))
})

test_that("normalize handles zero vector gracefully", {
  expect_warning(normalize(c(0, 0, 0)))
})

test_that("normalize returns unit vector", {
  v <- normalize(c(3, 4, 0))
  expect_equal(sqrt(sum(v^2)), 1, tolerance = 1e-10)
})

test_that("input validation catches dimension mismatches", {
  A <- matrix(1, 3, 3)
  B <- matrix(1, 4, 4)
  expect_error(inside_vertex_detection(A, B))
  expect_error(calculate_omega_overlap(A, B))
})

test_that("input validation catches non-square matrices in overlap", {
  A <- matrix(1, 3, 4)
  B <- matrix(1, 3, 4)
  expect_error(calculate_omega_overlap(A, B))
})

test_that("interaction_matrix_random validates inputs", {
  expect_error(interaction_matrix_random(num = 1, stren = 0.4, conne = 1))
  expect_error(interaction_matrix_random(num = 3, stren = -1, conne = 1))
  expect_error(interaction_matrix_random(num = 3, stren = 0.4, conne = 2))
})

test_that("interaction_matrix_random creates expected structure", {
  set.seed(1)
  A <- interaction_matrix_random(num = 4, stren = 0.4, conne = 1)
  expect_equal(nrow(A), 4)
  expect_equal(ncol(A), 4)
  expect_equal(diag(A), rep(-1, 4))
})

test_that("calculate_omega_overlap handles near-singular matrices", {
  A <- matrix(c(-1, 0.5, 0.5, -1), 2, 2)
  B <- diag(-1, 2)
  omega <- calculate_omega_overlap(A, B, method = "sphere")
  expect_true(is.numeric(omega))
  expect_true(omega >= 0)
})

test_that("calculate_omega does not corrupt global RNG state", {
  set.seed(123)
  r1 <- runif(1)
  set.seed(123)
  dummy <- calculate_omega(diag(3))
  r2 <- runif(1)
  expect_equal(r1, r2)
})

test_that("sphere method gives consistent result for single cone", {
  set.seed(1)
  A <- interaction_matrix_random(num = 3, stren = 0.4, conne = 1)
  omega1 <- calculate_omega(A, method = "sphere")
  omega2 <- calculate_omega(A, method = "sphere")
  expect_equal(omega1, omega2, tolerance = 0.01)
})

test_that("calculate_omega_constraint sphere method returns valid output", {
  set.seed(4)
  A <- interaction_matrix_random(3, 0.4, 1)
  C1 <- diag(c(-1, -1, -1), 3)
  omega <- calculate_omega_constraint(A, C1, method = "sphere")
  expect_true(is.numeric(omega))
  expect_true(omega >= 0)
})

test_that("calculate_omega_constraint methods agree", {
  set.seed(4)
  A <- interaction_matrix_random(3, 0.4, 1)
  C1 <- diag(c(-1, -1, -1), 3)
  omega_sphere <- calculate_omega_constraint(A, C1, method = "sphere")
  omega_mc <- calculate_omega_constraint(A, C1, method = "montecarlo", nsamples = 10000)
  expect_equal(omega_sphere, omega_mc, tolerance = 0.1)
})

test_that("2D case works correctly [sphere]", {
  set.seed(1)
  A <- interaction_matrix_random(num = 2, stren = 0.3, conne = 1)
  B <- diag(-1, 2)
  omega_A <- calculate_omega(A, method = "sphere")
  omega_overlap <- calculate_omega_overlap(A, B, method = "sphere")
  expect_true(omega_A >= 0)
  expect_true(omega_overlap >= 0)
  expect_true(omega_overlap <= max(omega_A, calculate_omega(B, method = "sphere")) + 0.05)
})

test_that("README example reproduces [sphere]", {
  set.seed(1)
  A <- interaction_matrix_random(num = 4, stren = 0.4, conne = 1)
  set.seed(2)
  B <- interaction_matrix_random(num = 4, stren = 0.4, conne = 1)

  # Self-overlap should equal omega
  omega_A <- calculate_omega(A, method = "sphere")
  omega_A_self <- calculate_omega_overlap(A, A, method = "sphere")
  expect_equal(omega_A, omega_A_self, tolerance = 0.02)

  # Overlap should be <= min(omega_A, omega_B)
  omega_B <- calculate_omega(B, method = "sphere")
  omega_AB <- calculate_omega_overlap(A, B, method = "sphere")
  expect_true(omega_AB <= min(omega_A, omega_B) + 0.05)
  expect_true(omega_AB >= 0)
})

test_that("method validation works", {
  A <- diag(-1, 3)
  B <- diag(-1, 3)
  expect_error(calculate_omega_overlap(A, B, method = "invalid"))
})

test_that("overlap is bounded by individual omegas [sphere, 5D]", {
  set.seed(100)
  A <- interaction_matrix_random(num = 5, stren = 0.3, conne = 1)
  set.seed(200)
  B <- interaction_matrix_random(num = 5, stren = 0.3, conne = 1)

  omega_A <- calculate_omega(A, method = "sphere")
  omega_B <- calculate_omega(B, method = "sphere")
  omega_AB <- calculate_omega_overlap(A, B, method = "sphere")

  expect_true(omega_AB <= min(omega_A, omega_B) + 0.05)
  expect_true(omega_AB >= 0)
})
