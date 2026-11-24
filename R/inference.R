#' Exact conditional score test for a single-index direction
#'
#' Tests \eqn{H_0: Y \perp X F_\beta \mid X\beta} at a candidate direction \code{beta}, as in
#' Section 5.3 of the paper. The smoothing tuple is chosen at the candidate from \code{(t, Y)}
#' only (regression bandwidth by protected leave-one-out squared error, then density
#' bandwidth by the criterion). Under Gaussian predictors the statistic has exactly the
#' \eqn{F_{p-1, n-p-1}} distribution conditionally on \eqn{(X\beta, Y)}, for any single-index
#' conditional response law.
#'
#' @param x Numeric predictor matrix, or a fitted \code{"sikdr"} object (then \code{y} is
#'   taken from the fit).
#' @param y Response vector (ignored when \code{x} is a fit).
#' @param beta Candidate direction (raw predictor scale), a vector or a matrix with one
#'   candidate per row.
#' @param level Confidence level used for the acceptance decision.
#' @return A data frame with the statistic, p-value and acceptance for each candidate.
#' @examples
#' set.seed(2)
#' n <- 200
#' x <- matrix(rnorm(n * 3), n, 3)
#' y <- sin(2 * drop(x %*% c(1, 1, 0)) / sqrt(2)) + rexp(n) - 1
#' sikdr_test(x, y, rbind(c(1, 1, 0), c(1, 0, 0)))
#' @export
sikdr_test <- function(x, y = NULL, beta, level = 0.95) {
  if (inherits(x, "sikdr")) {
    y <- x$y
    x <- x$x
  }
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  q <- ncol(x)
  dd <- q - 1
  if (n <= q + 1) stop("the test needs n > p + 1")
  B <- if (is.matrix(beta)) beta else matrix(beta, nrow = 1)
  if (ncol(B) != q) stop("beta has the wrong length")
  ys <- (y - mean(y)) / sqrt(mean((y - mean(y))^2))
  fg <- sikdr_grid(n, "final")
  crit <- stats::qf(level, dd, n - dd - 2)
  res <- t(apply(B, 1, function(b) candidate_stat(x, ys, b, fg)))
  data.frame(F = res[, 1], df1 = dd, df2 = n - dd - 2,
             p.value = stats::pf(res[, 1], dd, n - dd - 2, lower.tail = FALSE),
             accepted = res[, 1] <= crit | res[, 2] == 1, g = res[, 3], h = res[, 4])
}

candidate_stat <- function(x, ys, b, fg) {
  n <- nrow(x)
  b <- unit_vec(b)
  xb <- drop(x %*% b)
  xc <- xb - mean(xb)
  s <- sqrt(mean(xc * xc))
  if (!is.finite(s) || s <= 0) return(c(NA_real_, 0, NA_real_, NA_real_))
  tt <- xc / s
  W <- x %*% householder_frame(b)
  Zm <- sweep(W, 2, colMeans(W))
  Zm <- Zm - outer(tt, colMeans(tt * Zm))
  sel <- select_g_then_h(tt, ys, fg$g, fg$h, fg$lambda, fg$rho)
  cc <- kdr_score_coef(tt, ys, sel$g, sel$h, fg$lambda, fg$rho)$c
  u <- cc - mean(cc) - tt * mean(tt * cc)
  un <- sqrt(sum(u * u))
  if (un == 0) return(c(0, 1, sel$g, sel$h))
  z <- drop(crossprod(Zm, u / un))
  A <- crossprod(Zm) - tcrossprod(z)
  dd <- ncol(Zm)
  Fs <- (n - dd - 2) / dd * sum(z * solve(A, z))
  c(Fs, 0, sel$g, sel$h)
}

#' Confidence set for the direction by score inversion
#'
#' Draws candidate directions uniformly from the cap of angular radius
#' \eqn{\log(n + e)/\sqrt{n}} around the estimate, keeps those inside the paper's domain
#' (the cap and the chart of the Fourier initial estimator) and returns the candidates
#' accepted by \code{\link{sikdr_test}}. The set is empty when the inference conditions of the
#' fit fail.
#'
#' @param fit A \code{"sikdr"} object fitted with \code{covariance = TRUE}.
#' @param draws Number of candidate directions.
#' @param level Confidence level.
#' @return A list with the accepted directions (rows), the fraction of accepted draws and the
#'   largest accepted angle (radians) from the estimate.
#' @export
sikdr_confset <- function(fit, draws = 500, level = 0.95) {
  if (!inherits(fit, "sikdr") || is.null(fit$inference_ok)) stop("fit with covariance = TRUE first")
  if (!fit$inference_ok) {
    return(list(accepted = matrix(numeric(0), 0, length(fit$beta)), fraction = 0, max_angle = NA,
                empty_reason = fit$inference_reasons))
  }
  n <- fit$n
  q <- length(fit$beta)
  dd <- q - 1
  a <- log(n + exp(1)) / sqrt(n)
  bh <- fit$beta_internal
  Fr <- householder_frame(bh)
  # uniform draws on the spherical cap of radius a (density proportional to sin^{d-1})
  Ia <- stats::pbeta(sin(a)^2, dd / 2, 0.5)
  theta <- asin(sqrt(stats::qbeta(stats::runif(draws) * Ia, dd / 2, 0.5)))
  dirs <- matrix(stats::rnorm(draws * dd), draws, dd)
  dirs <- dirs / sqrt(rowSums(dirs^2))
  cand <- t(apply(cbind(theta, dirs), 1, function(r) {
    unit_vec(cos(r[1]) * bh + sin(r[1]) * drop(Fr %*% r[-1]))
  }))
  inside <- apply(cand, 1, in_initial_cap, fit = fit)
  keep <- cand[inside, , drop = FALSE]
  ok <- rep(FALSE, nrow(keep))
  if (nrow(keep) > 0) ok <- sikdr_test(fit, beta = keep, level = level)$accepted
  acc <- keep[ok, , drop = FALSE]
  list(accepted = acc, fraction = mean(c(ok, rep(FALSE, sum(!inside)))),
       max_angle = if (nrow(acc)) max(apply(acc, 1, angle_lines, v = bh)) else NA_real_,
       radius = a)
}

in_initial_cap <- function(b, fit) {
  gam <- unit_vec(drop(fit$std$Shalf %*% b))
  c0 <- sum(fit$gamma0 * gam)
  if (c0 == 0) return(FALSE)
  zc <- drop(crossprod(fit$U, gam)) / c0
  sqrt(sum(zc * zc)) <= ascent_constants(fit$n)$R
}
