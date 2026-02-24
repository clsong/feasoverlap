# nolint start
#' function that generates a random matrix (as in May, Nature (1972))
#' @param num number of species
#' @param stren standard deviation of interaction strength
#' @param conne connectance of the interaction matrix
#' @return the generated random matrix
#' @importFrom stats rnorm
#' @export
interaction_matrix_random <- function(num, stren, conne) {
  stopifnot(
    "num must be a positive integer" = is.numeric(num) && length(num) == 1 && num >= 2,
    "stren must be a positive number" = is.numeric(stren) && stren > 0,
    "conne must be between 0 and 1" = is.numeric(conne) && conne >= 0 && conne <= 1
  )
  Inte <- rnorm(num * num, mean = 0, sd = stren)
  zeroes <- sample(c(rep.int(1, floor(num * num * conne)), rep.int(0, (num * num - floor(num * num * conne)))))
  Inte[which(zeroes == 0)] <- 0
  Inte <- matrix(Inte, ncol = num, nrow = num)
  diag(Inte) <- -1
  return(Inte)
}

#' Internal helper: approximate equality check (replaces dplyr::near)
#' @param x numeric vector
#' @param y numeric vector
#' @param tol tolerance
#' @return logical vector
#' @keywords internal
.near <- function(x, y, tol = .Machine$double.eps^0.5) {
  abs(x - y) < tol
}

#' Internal helper: sample uniformly on the unit sphere
#' Replaces uniformly::runif_on_sphere to avoid the rgl dependency
#' that causes X11/GLX warnings on headless systems.
#' @param n number of points to sample
#' @param d dimension of the sphere
#' @param r radius (default 1)
#' @return n x d matrix of points on the sphere
#' @importFrom stats rnorm
#' @keywords internal
.runif_on_sphere <- function(n, d, r = 1) {
  x <- matrix(rnorm(n * d), nrow = n, ncol = d)
  norms <- sqrt(rowSums(x^2))
  r * x / norms
}

#' function that computes the normalized feasibility from an interaction matrix
#' @import geometry
#' @import mvtnorm
#' @importFrom stats runif
#' @param vertex all the vertexes of the feasibility domain
#' @param raw TRUE: raw omega, FALSE: normalized omega
#' @param nsamples number of sampled points
#' @param method either 'convex_hull' or 'sphere'
#' @param seed random seed for reproducibility (default: 1010). Set to NULL to use current RNG state.
#' @return the normalized feasibility
#' @export
calculate_omega <- function(vertex, raw = FALSE, nsamples = 100,
                            method = "convex_hull", seed = 1010) {
  num <- nrow(vertex)
  vertex <- generate_span_vectors(vertex)

  if (method == "convex_hull") {
    # Scope the RNG so we don't corrupt the global seed
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    if (!is.null(seed)) set.seed(seed)

    vertex <- cbind(
      vertex,
      vertex %*% t(abs(.runif_on_sphere(n = nsamples, d = ncol(vertex), r = 1)))
    )
    if (num < 5) {
      vertex <- generate_span_vectors(vertex) %*% diag(
        runif(
          ncol(vertex),
          (1 - .05 * (num - 2)),
          (1 + .05 * (num - 2))
        )
      )
    } else {
      vertex <- generate_span_vectors(vertex) %*% diag(
        runif(
          ncol(vertex),
          (1 - .05 * (num - 2)),
          (1 + .1 * (num - 2))
        )
      )
    }

    # Restore the global RNG state
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }

    vertex <- cbind(vertex, rep(0, num))

    vol_ori <- (convhulln(t(vertex), output.options = TRUE)$vol)
    vol_ball <- (pi^(num / 2) / gamma(num / 2 + 1))

    omega <- ifelse(raw == FALSE,
      (vol_ori / vol_ball)^(1 / num),
      vol_ori / vol_ball
    )
  }
  if (method == "sphere") {
    d <- pmvnorm(
      lower = rep(0, num),
      upper = rep(Inf, num),
      mean = rep(0, num), sigma = solve(t(vertex) %*% vertex)
    )
    omega <- ifelse(raw == FALSE,
      d[1]^(1 / num),
      d[1]
    )
  }
  omega
}

#' function that normalizes a vector in the L2 norm
#' @param a the original vector
#' @return the normalized vector
#' @export
normalize <- function(a) {
  norm_val <- sqrt(sum(a^2))
  if (norm_val < .Machine$double.eps) {
    warning("Attempting to normalize a near-zero vector")
    return(a)
  }
  a / norm_val
}

#' function that normalizes the spanning vectors of the feasibility domain in the L2 norm
#' @param alpha interaction matrix
#' @return the normalized spanning vectors
#' @export
generate_span_vectors <- function(alpha) {
  apply(alpha, 2, function(x) normalize(x))
}

#' function that computes all the extreme points that belong to original vertexes
#' @param A one interaction matrix
#' @param B another interaction matrix
#' @return all the extreme points that belong to original vertexes
#' @export
inside_vertex_detection <- function(A, B) {
  stopifnot(
    "A and B must have the same dimensions" = all(dim(A) == dim(B)),
    "A must be a square matrix" = nrow(A) == ncol(A)
  )
  SpanA <- generate_span_vectors(A)
  SpanB <- generate_span_vectors(B)
  # to determine whether a vertex of one cone is inside another cone or not.
  inside_detection <- function(Span, vector) {
    lambda <- tryCatch(
      qr.solve(Span, vector),
      error = function(e) rep(-1, ncol(Span))
    )
    if (sum(lambda >= -1e-10) == length(lambda)) {
      return(1)
    } else {
      return(0)
    }
  }
  inside_vertex <- list()
  l <- 1
  for (i in 1:ncol(B)) {
    auxi <- inside_detection(SpanA, SpanB[, i])
    if (auxi == 1) {
      inside_vertex[[l]] <- SpanB[, i]
      l <- l + 1
    }
  }
  for (i in 1:ncol(A)) {
    auxi <- inside_detection(SpanB, SpanA[, i])
    if (auxi == 1) {
      inside_vertex[[l]] <- SpanA[, i]
      l <- l + 1
    }
  }
  return(inside_vertex)
}

#' function that computes all the extreme points generated that are generated by the intersections of the cones
#' @param S one interaction matrix
#' @param M another interaction matrix
#' @importFrom utils combn
#' @importFrom pracma nullspace Rank
#' @return all the extreme points that are generated by the intersections of the cones
#' @export
intersection_vertex_detection <- function(S, M) {
  num <- ncol(S)
  if (num == 2) {
    return(list())
  } else {
    combination_S <- combn(1:ncol(S), 2)
    combination_M <- combn(1:ncol(S), (num - 1))
    Span_S <- generate_span_vectors(S)
    Span_M <- generate_span_vectors(M)

    border_M <- list()
    extreme_point_M <- list()
    for (i in 1:ncol(M)) {
      coeff_matrix <- matrix(0, ncol = num, nrow = num - 1)
      for (j in 1:(num - 1)) {
        coeff_matrix[j, ] <- Span_M[, combination_M[j, i]]
      }
      if (Rank(coeff_matrix) == (num - 1)) {
        border_M[[i]] <- nullspace(coeff_matrix)
      } else {
        # Degenerate face — skip instead of crashing
        border_M[[i]] <- NULL
        next
      }
      extreme_point_M[[i]] <- t(coeff_matrix)[1:(num - 1), 1:(num - 1)]
    }

    inside_face_detection <- function(extreme_point, test_vector) {
      lambda <- tryCatch(
        qr.solve(extreme_point, test_vector),
        error = function(e) rep(-1, ncol(extreme_point))
      )
      if (sum(lambda >= -1e-10) == length(lambda)) {
        return(1)
      } else {
        return(0)
      }
    }

    l <- 1
    intersection_vertex <- list()
    side <- c()
    for (i in 1:ncol(combination_S)) {
      vertex_1 <- Span_S[, combination_S[1, i]]
      vertex_2 <- Span_S[, combination_S[2, i]]
      for (j in seq_along(border_M)) {
        # Skip degenerate faces
        if (is.null(border_M[[j]])) next

        n1 <- sum(vertex_1 * border_M[[j]])
        n2 <- sum(vertex_2 * border_M[[j]])

        auxi <- n1 * n2
        if (auxi < -1e-10) {
          lambda <- n2 / (n2 - n1)
          possible <- lambda * vertex_1 + (1 - lambda) * vertex_2

          ep <- extreme_point_M[[j]]
          det_val <- tryCatch(det(ep), error = function(e) 0)
          if (det_val != 0) {
            auxi2 <- inside_face_detection(ep, possible[1:(num - 1)])
            if (auxi2 == 1) {
              intersection_vertex[[l]] <- possible
              side[l] <- j
              l <- l + 1
            }
          }
        }
      }
    }

    if (length(intersection_vertex) > 0) {
      for (i in 1:length(intersection_vertex)) {
        intersection_vertex[[i]] <- normalize(intersection_vertex[[i]])
      }
    }
  }

  return(intersection_vertex)
}

#' function that computes all the extreme points
#' @param A one interaction matrix
#' @param B another interaction matrix
#' @return all the extreme points that generate the intersection region
#' @export
vertex_detection <- function(A, B) {
  num <- ncol(A)
  inside_vertex <- inside_vertex_detection(A, B)
  intersection_vertex <- intersection_vertex_detection(A, B)

  # combine the two vertex lists
  if (length(inside_vertex) > 0) {
    vertex <- matrix(unlist(inside_vertex), nrow = num, byrow = FALSE)
  } else {
    vertex <- matrix(0, nrow = num, ncol = 2)
  }
  if (length(intersection_vertex) > 0) {
    vertex <- cbind(vertex, matrix(unlist(intersection_vertex), nrow = num, byrow = FALSE))
  }

  # delete the points that are nonzero due to numerical error
  delete_zeroes <- c()
  for (i in 1:ncol(vertex)) {
    if (.near(sum(vertex[, i]^2), 0)) {
      delete_zeroes <- c(delete_zeroes, i)
    }
  }
  if (length(delete_zeroes) > 0) vertex <- vertex[, -delete_zeroes, drop = FALSE]


  # delete the same ones
  if (length(vertex) > num) {
    for (test in 1:ncol(vertex)) {
      vertex[, test] <- normalize(vertex[, test])
    }
    delete_duplicates <- c()
    for (i in 1:(ncol(vertex) - 1)) {
      for (j in (i + 1):ncol(vertex)) {
        if (sum(.near(vertex[, i], vertex[, j])) == nrow(vertex)) {
          delete_duplicates <- c(delete_duplicates, j)
        }
      }
    }
    if (length(delete_duplicates) > 0) vertex <- vertex[, -unique(delete_duplicates), drop = FALSE]
  }
  return(vertex)
}

#' function that computes the overlap of two feasibility domains
#'
#' Three methods are available:
#' \itemize{
#'   \item \code{"sphere"} (default): Uses multivariate normal integration via
#'     \code{pmvnorm}. Constructs the 2n-dimensional constraint matrix
#'     \eqn{C = [-A^{-1}; -B^{-1}]} and computes
#'     \eqn{P(C r > 0)} where \eqn{r \sim N(0, I)}.
#'     This is deterministic, symmetric, and avoids vertex enumeration entirely.
#'   \item \code{"montecarlo"}: Samples uniformly on the unit sphere and checks
#'     membership in both feasibility cones. Good for high dimensions (n > 15).
#'   \item \code{"convex_hull"}: Legacy method using vertex enumeration and Qhull
#'     convex hull computation. Kept for backward compatibility.
#' }
#'
#' @param A one interaction matrix
#' @param B another interaction matrix
#' @param raw TRUE: raw omega, FALSE: normalized omega
#' @param nsamples number of sampled points (for convex_hull and montecarlo methods)
#' @param method one of 'sphere', 'montecarlo', or 'convex_hull'
#' @return the normalized feasibility of the intersection region
#' @export
calculate_omega_overlap <- function(A, B, raw = FALSE, nsamples = 1000,
                                     method = "sphere") {
  stopifnot(
    "A and B must have the same dimensions" = all(dim(A) == dim(B)),
    "A must be a square matrix" = nrow(A) == ncol(A),
    "method must be one of 'sphere', 'montecarlo', or 'convex_hull'" =
      method %in% c("sphere", "montecarlo", "convex_hull")
  )
  num <- nrow(A)

  if (method == "sphere") {
    # Invert both matrices; fall back to montecarlo if singular
    A_inv <- tryCatch(solve(A), error = function(e) NULL)
    B_inv <- tryCatch(solve(B), error = function(e) NULL)

    if (is.null(A_inv) || is.null(B_inv)) {
      # Fall back to Monte Carlo for singular matrices
      return(calculate_omega_overlap(A, B, raw = raw, nsamples = nsamples,
                                      method = "montecarlo"))
    }

    # Construct the 2n x n constraint matrix C = [-A^{-1}; -B^{-1}]
    C <- rbind(-A_inv, -B_inv)

    # Covariance of y = Cr where r ~ N(0, I): Sigma = C C^T
    Sigma <- C %*% t(C)

    # P(y > 0) gives the raw overlap omega
    d <- tryCatch(
      pmvnorm(
        lower = rep(0, 2 * num),
        upper = rep(Inf, 2 * num),
        mean = rep(0, 2 * num),
        sigma = Sigma
      ),
      error = function(e) {
        # If pmvnorm fails (e.g. Sigma not positive definite), fall back
        return(NULL)
      }
    )

    if (is.null(d)) {
      return(calculate_omega_overlap(A, B, raw = raw, nsamples = nsamples,
                                      method = "montecarlo"))
    }

    omega <- ifelse(raw, d[1], d[1]^(1 / num))
    return(omega)
  }

  if (method == "montecarlo") {
    # Sample uniformly on the unit sphere and check membership in both cones
    A_inv <- tryCatch(solve(A), error = function(e) NULL)
    B_inv <- tryCatch(solve(B), error = function(e) NULL)

    if (is.null(A_inv) || is.null(B_inv)) {
      return(0)
    }

    # Sample points uniformly on the unit sphere: nsamples x n matrix
    samples <- .runif_on_sphere(n = nsamples, d = num, r = 1)

    # Vectorized: compute all transformations at once via matrix multiply
    # Na = -A_inv %*% t(samples) -> n x nsamples; each column = one sample
    # Nb = -B_inv %*% t(samples) -> n x nsamples
    Na <- -A_inv %*% t(samples)
    Nb <- -B_inv %*% t(samples)

    # Check if all components > -tol for each sample (column)
    in_A <- colSums(Na >= -1e-10) == num
    in_B <- colSums(Nb >= -1e-10) == num
    in_both <- in_A & in_B

    raw_omega <- mean(in_both)
    omega <- ifelse(raw, raw_omega, raw_omega^(1 / num))
    return(omega)
  }

  if (method == "convex_hull") {
    # Legacy method: vertex enumeration + convex hull
    overlap_vertex <- tryCatch(
      {
        v <- vertex_detection(A, B) %>%
          cbind(vertex_detection(B, A)) %>%
          unique(MARGIN = 2)
        v
      },
      error = function(cond) {
        matrix(0, nrow = num, ncol = 1)
      }
    )

    if (ncol(overlap_vertex) < num || qr(overlap_vertex)$rank < num) {
      volume_overlap <- 0
    } else {
      volume_overlap <- tryCatch(
        {
          calculate_omega(overlap_vertex, raw, nsamples)
        },
        error = function(cond) {
          0
        }
      )
    }

    return(volume_overlap)
  }
}

#' function that computes the constrained feasibility domain size
#'
#' Three methods are available:
#' \itemize{
#'   \item \code{"sphere"} (default): Uses multivariate normal integration. For a
#'     diagonal constraint matrix B, constructs the constraint
#'     \eqn{C = [-A^{-1}; B]} and computes \eqn{P(Cr > 0)}.
#'   \item \code{"montecarlo"}: Samples growth rate vectors and checks feasibility
#'     under the constraint. Good for high dimensions.
#'   \item \code{"legacy"}: Original sampling method kept for backward compatibility.
#' }
#'
#' @param A one interaction matrix
#' @param B the constraint matrix
#' @param raw TRUE: raw omega, FALSE: normalized omega
#' @param nsamples number of sampled points (for montecarlo and legacy methods)
#' @param method one of 'sphere', 'montecarlo', or 'legacy'

#' @return the normalized feasibility of the intersection region
#' @export
calculate_omega_constraint <- function(A, B, raw = FALSE, nsamples,
                                        method = "sphere") {
  stopifnot(
    "A must be a square matrix" = nrow(A) == ncol(A),
    "B must be a square matrix" = nrow(B) == ncol(B),
    "A and B must have the same dimensions" = all(dim(A) == dim(B))
  )
  num <- nrow(A)

  if (method == "sphere") {
    # Use pmvnorm with the 2n-dim constraint C = [-A^{-1}; B]
    A_inv <- tryCatch(solve(A), error = function(e) NULL)

    if (is.null(A_inv)) {
      return(calculate_omega_constraint(A, B, raw = raw, nsamples = nsamples,
                                         method = "montecarlo"))
    }

    C <- rbind(-A_inv, B)
    Sigma <- C %*% t(C)

    d <- tryCatch(
      pmvnorm(
        lower = rep(0, 2 * num),
        upper = rep(Inf, 2 * num),
        mean = rep(0, 2 * num),
        sigma = Sigma
      ),
      error = function(e) NULL
    )

    if (is.null(d)) {
      return(calculate_omega_constraint(A, B, raw = raw, nsamples = nsamples,
                                         method = "montecarlo"))
    }

    omega <- ifelse(raw, d[1], d[1]^(1 / num))
    return(omega)
  }

  if (method == "montecarlo") {
    A_inv <- tryCatch(solve(A), error = function(e) NULL)
    if (is.null(A_inv)) return(0)

    if (missing(nsamples)) nsamples <- max(2^num * 250, 1000)

    # Vectorized sampling and membership check
    samples <- .runif_on_sphere(n = nsamples, d = num, r = 1)
    Na <- -A_inv %*% t(samples)
    Br <- B %*% t(samples)

    in_A <- colSums(Na >= -1e-10) == num
    in_B <- colSums(Br >= -1e-10) == num
    in_both <- in_A & in_B

    raw_omega <- mean(in_both)
    omega <- ifelse(raw, raw_omega, raw_omega^(1 / num))
    return(omega)
  }

  if (method == "legacy") {
    sgn <- -diag(B)

    if (missing(nsamples)) {
      nsamples <- 2^num * 250
    }

    # Pre-compute A inverse once, not per sample
    A_inv <- tryCatch(solve(A), error = function(e) NULL)
    if (is.null(A_inv)) return(0)

    # Vectorized: generate all samples at once
    R <- matrix(rnorm(num * nsamples), nrow = num, ncol = nsamples)
    R <- R / rep(sqrt(colSums(R^2)), each = num)  # normalize columns
    R <- abs(R) * sgn  # apply sign constraint

    # Check feasibility: N* = -A^{-1} r > 0 for all components
    N <- -A_inv %*% R  # n x nsamples
    feasibility <- colSums(N >= 0) == num

    volume_overlap <- ifelse(raw,
      mean(feasibility),
      mean(feasibility)^(1 / num)
    )

    volume_overlap * calculate_omega(B, method = 'sphere')
  }
}
