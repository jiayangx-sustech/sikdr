# Tuning grids, ascent constants, standardization, frames and the Fourier initial estimator.

#' Smoothing grids used by SI-KDR
#'
#' Returns the regression bandwidths \code{g}, density bandwidths \code{h}, ridge
#' \code{lambda} and floor \code{rho} of the paper's grids for a sample of size \code{n}.
#' \code{"final"} is the twenty-tuple grid ranked by the covariance trace, \code{"prelim"}
#' the grid of the preliminary fit and \code{"reference"} the reference rule (SI-KDR-R).
#'
#' @param n Sample size used for fitting.
#' @param type One of \code{"final"}, \code{"prelim"} or \code{"reference"}.
#' @return A list with elements \code{g}, \code{h}, \code{lambda} and \code{rho}.
#' @examples
#' sikdr_grid(600)
#' @export
sikdr_grid <- function(n, type = c("final", "prelim", "reference")) {
  type <- match.arg(type)
  if (type == "reference") {
    return(list(g = n^(-4 / 25) * c(1.2, 1.4, 1.6), h = n^(-1 / 30) * c(0.4, 0.6, 0.8),
                lambda = 0.01 * n^(-7 / 10), rho = 1 / n))
  }
  gm <- if (type == "final") c(0.4, 0.6, 0.8, 1.2) else c(1.2, 1.4, 1.6)
  list(g = 600^(-4 / 25) * gm * (n / 600)^(-1 / 6),
       h = 600^(-1 / 30) * c(0.08, 0.12, 0.18, 0.27, 0.4) * (n / 600)^(-1 / 15),
       lambda = 0.01 * 600^(-7 / 10) * (n / 600)^(-3 / 4), rho = 1 / n)
}

ascent_constants <- function(n) {
  list(R = min(1.5, 2 * n^(-1 / 20)), max_iter = ceiling(log(n + exp(1))^3),
       step0 = 1 / log(n + exp(1)), max_backtracks = 80L, grad_tol = 1 / n)
}

unit_vec <- function(v) v / sqrt(sum(v * v))

# Standardization with divisor n: Z = Sigma^{-1/2}(X - Xbar), y = (Y - Ybar)/sY.
standardize_xy <- function(X, Y) {
  n <- nrow(X)
  q <- ncol(X)
  xbar <- colMeans(X)
  Xc <- sweep(X, 2, xbar)
  S <- crossprod(Xc) / n
  S <- (S + t(S)) / 2
  ybar <- mean(Y)
  sy <- sqrt(mean((Y - ybar)^2))
  if (!all(is.finite(S)) || !is.finite(sy)) stop("non-finite predictor or response moments")
  if (sy <= 0) stop("the response has zero variance")
  e <- eigen(S, symmetric = TRUE)
  lam <- e$values
  V <- e$vectors
  if (min(lam) <= 0 || min(lam) <= q * .Machine$double.eps * max(max(lam), 0)) {
    stop("the predictor covariance is singular; remove collinear predictors or screen first")
  }
  A <- V %*% (t(V) / sqrt(lam))
  Sh <- V %*% (t(V) * sqrt(lam))
  A <- (A + t(A)) / 2
  Sh <- (Sh + t(Sh)) / 2
  list(Z = Xc %*% A, y = (Y - ybar) / sy, xbar = xbar, Sigma = S, A = A, Shalf = Sh,
       ybar = ybar, sy = sy)
}

# Deterministic orthonormal basis of the orthogonal complement of the line [u].
householder_frame <- function(u) {
  u <- unit_vec(u)
  s <- if (u[1] >= 0) 1 else -1
  w <- u
  w[1] <- w[1] + s
  H <- diag(length(u)) - 2 * tcrossprod(w) / sum(w * w)
  H[, -1, drop = FALSE]
}

deterministic_sign <- function(v) {
  k <- which.max(abs(v))
  if (v[k] >= 0) v else -v
}

raw_from_working <- function(gamma, A) unit_vec(drop(crossprod(A, gamma)))

angle_lines <- function(u, v) {
  u <- unit_vec(u)
  v <- unit_vec(v)
  cc <- sum(u * v)
  atan2(sqrt(sum((u - cc * v)^2)), abs(cc))
}

# Fourier initial estimator: average of four matrices (a = 1/q, 1, 3, 9; b = 1), all ordered pairs.
fourier_direction <- function(Z, y) {
  n <- nrow(Z)
  q <- ncol(Z)
  D2 <- matrix(0, n, n)
  for (k in seq_len(q)) D2 <- D2 + outer(Z[, k], Z[, k], "-")^2
  Dy2 <- outer(y, y, "-")^2
  M <- matrix(0, q, q)
  for (a in c(1 / q, 1, 3, 9)) {
    K <- safe_exp(-a * D2 / 2 - Dy2 / 2)
    K1 <- rowSums(K)
    Ma <- a * sum(K1) * diag(q) + ((1 + a)^2 + a^2) * crossprod(Z, K %*% Z) -
      2 * a * (1 + a) * crossprod(Z * K1, Z)
    M <- M + Ma / n^2
  }
  M <- M / 4
  M <- (M + t(M)) / 2
  if (!all(is.finite(M))) stop("the Fourier initial matrix is not finite")
  e <- eigen(M, symmetric = TRUE)
  tie <- q > 1 && abs(e$values[1] - e$values[2]) <= 8 * .Machine$double.eps * max(1, abs(e$values[1]))
  gamma <- if (tie) replace(numeric(q), 1, 1) else deterministic_sign(e$vectors[, 1])
  list(gamma = gamma, tie = tie, values = e$values)
}
