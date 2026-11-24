#' Single-Index Kernel Density Regression
#'
#' Fits the single-index model \eqn{Y = m(\beta^T X) + \varepsilon} with the error
#' independent of the predictors. At each candidate direction the link is refitted by a
#' protected leave-one-out Nadaraya--Watson smoother and the direction is scored by the
#' leave-one-out kernel log-likelihood of the residuals.
#'
#' With \code{tuning = "trace"} (the default) the twenty smoothing tuples of
#' \code{sikdr_grid(n, "final")} are ranked by the trace of the case-weight covariance
#' estimate at a preliminary fit, and the tuple with the smallest trace among those with
#' positive curvature is used. \code{tuning = "reference"} uses the reference rule (SI-KDR-R).
#' The optimizer is the projected Armijo ascent started from the Fourier initial estimator.
#'
#' @param x Numeric matrix of predictors (n rows, at least two columns).
#' @param y Numeric response vector of length n.
#' @param tuning \code{"trace"} or \code{"reference"}.
#' @param covariance If \code{TRUE}, compute the case-weight covariance at the estimate and
#'   check the conditions used for inference.
#' @return An object of class \code{"sikdr"} with the unit direction \code{beta} (raw
#'   predictor scale, sign fixed so that the largest coordinate is positive), the smoothing
#'   tuple \code{tuple = c(g, h, lambda, rho)}, the covariance \code{V} of
#'   \eqn{\sqrt{n}} times the tangent coordinates in the frame \code{frame}, and diagnostics.
#' @references Xu, J., Jia, Y., Chen, Y. and Chen, X. Finite-smoothing efficiency and
#'   inference in single-index models. arXiv:2611.22384.
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- matrix(rnorm(n * 3), n, 3)
#' e <- ifelse(runif(n) < 0.7, rnorm(n, -0.9, 0.5), rnorm(n, 2.1, 0.3))
#' y <- sin(2 * drop(x %*% c(1, 1, 0)) / sqrt(2)) + e
#' fit <- sikdr(x, y, tuning = "reference", covariance = FALSE)
#' fit$beta
#' @export
sikdr <- function(x, y, tuning = c("trace", "reference"), covariance = TRUE) {
  tuning <- match.arg(tuning)
  cl <- match.call()
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  y <- as.numeric(y)
  if (nrow(x) != length(y)) stop("x and y have different numbers of observations")
  if (ncol(x) < 2) stop("x needs at least two columns")
  if (!all(is.finite(x)) || !all(is.finite(y))) stop("x and y must be finite")
  n <- nrow(x)
  std <- standardize_xy(x, y)
  four <- fourier_direction(std$Z, std$y)
  g0 <- four$gamma
  U <- householder_frame(g0)
  T0 <- drop(std$Z %*% g0)
  diagn <- list(fourier_tie = four$tie)
  if (tuning == "trace") {
    pg <- sikdr_grid(n, "prelim")
    s0 <- select_g_then_h(T0, std$y, pg$g, pg$h, pg$lambda, pg$rho)
    tup_p <- c(s0$g, s0$h, pg$lambda, pg$rho)
    pre <- projected_armijo(std$Z, std$y, g0, U, tup_p)
    if (pre$status == "nonfinite") stop("the preliminary fit failed (non-finite criterion)")
    beta_pre <- raw_from_working(pre$gamma, std$A)
    Fp <- householder_frame(beta_pre)
    fg <- sikdr_grid(n, "final")
    tuples <- as.matrix(expand.grid(g = fg$g, h = fg$h))
    traces <- rep(NA_real_, nrow(tuples))
    pd <- logical(nrow(tuples))
    for (k in seq_len(nrow(tuples))) {
      ij <- ij_covariance(x, y, beta_pre, Fp, c(tuples[k, 1], tuples[k, 2], fg$lambda, fg$rho))
      traces[k] <- ij$trace
      pd[k] <- isTRUE(ij$J_pd) && is.finite(ij$trace)
    }
    if (any(pd)) {
      kbest <- which(pd)[which.min(traces[pd])]
      tup <- c(tuples[kbest, 1], tuples[kbest, 2], fg$lambda, fg$rho)
      fallback <- FALSE
    } else {
      s1 <- select_g_then_h(drop(std$Z %*% pre$gamma), std$y, fg$g, fg$h, fg$lambda, fg$rho)
      tup <- c(s1$g, s1$h, fg$lambda, fg$rho)
      kbest <- which(tuples[, 1] == s1$g & tuples[, 2] == s1$h)[1]
      fallback <- TRUE
    }
    diagn <- c(diagn, list(prelim_tuple = tup_p, beta_prelim = beta_pre, traces = traces,
                           trace_pd = pd, selected_index = kbest, tuning_fallback = fallback))
  } else {
    rg <- sikdr_grid(n, "reference")
    s0 <- select_g_then_h(T0, std$y, rg$g, rg$h, rg$lambda, rg$rho)
    tup <- c(s0$g, s0$h, rg$lambda, rg$rho)
  }
  names(tup) <- c("g", "h", "lambda", "rho")
  fin <- projected_armijo(std$Z, std$y, g0, U, tup)
  if (fin$status == "nonfinite") stop("the final fit failed (non-finite criterion)")
  beta <- raw_from_working(fin$gamma, std$A)
  sgn <- if (beta[which.max(abs(beta))] >= 0) 1 else -1
  out <- list(beta = sgn * beta, beta_internal = beta, tuple = tup, tuning = tuning,
              criterion = fin$F, status = fin$status, iterations = fin$iterations,
              interior = fin$interior, gamma = fin$gamma, gamma0 = g0, U = U, z = fin$z,
              grad = fin$grad, std = std[c("xbar", "A", "Shalf", "ybar", "sy")],
              diagnostics = diagn, x = x, y = y, n = n, call = cl)
  if (!is.null(colnames(x))) names(out$beta) <- colnames(x)
  if (covariance) out <- c(out, inference_conditions(out))
  class(out) <- "sikdr"
  out
}

# Interior solution, small working score, negative-definite working Hessian, and positive-
# definite case-weight J and B at the estimate.
inference_conditions <- function(fit) {
  n <- fit$n
  reasons <- character(0)
  if (!fit$interior) reasons <- c(reasons, "not_interior")
  gnorm <- sqrt(sum(fit$grad^2))
  if (!is.finite(gnorm) || gnorm > n^(-1 / 2) / log(n + exp(1))) {
    reasons <- c(reasons, "working_score_too_large")
  }
  Zs <- fit$std$Z
  if (is.null(Zs)) Zs <- sweep(fit$x, 2, fit$std$xbar) %*% fit$std$A
  ys <- (fit$y - fit$std$ybar) / fit$std$sy
  d <- length(fit$z)
  Hw <- matrix(0, d, d)
  eps <- 1e-5
  for (k in seq_len(d)) {
    e <- replace(numeric(d), k, eps)
    gp <- working_value_grad(fit$z + e, Zs, ys, fit$gamma0, fit$U, fit$tuple)$grad
    gm <- working_value_grad(fit$z - e, Zs, ys, fit$gamma0, fit$U, fit$tuple)$grad
    Hw[, k] <- (gp - gm) / (2 * eps)
  }
  Hw <- (Hw + t(Hw)) / 2
  if (!all(is.finite(Hw)) || max(eigen(Hw, symmetric = TRUE, only.values = TRUE)$values) >= 0) {
    reasons <- c(reasons, "working_curvature_not_negative")
  }
  Fr <- householder_frame(fit$beta_internal)
  ij <- ij_covariance(fit$x, fit$y, fit$beta_internal, Fr, fit$tuple)
  if (!ij$finite) reasons <- c(reasons, "covariance_nonfinite")
  if (!isTRUE(ij$J_pd)) reasons <- c(reasons, "J_not_pd")
  if (!isTRUE(ij$B_pd)) reasons <- c(reasons, "B_not_pd")
  list(inference_ok = length(reasons) == 0, inference_reasons = reasons, frame = Fr, J = ij$J,
       B = ij$B, V = ij$V, working_score_norm = gnorm)
}

#' @export
coef.sikdr <- function(object, ...) object$beta

#' Covariance of the estimated direction
#'
#' Returns the case-weight estimate of the covariance of the unit direction on the raw
#' predictor scale, \eqn{F V F^T / n}, where \eqn{F} spans the tangent space at the estimate.
#' The matrix has rank \eqn{p - 1} because the direction is a unit vector.
#'
#' @param object A fitted \code{"sikdr"} object with \code{covariance = TRUE}.
#' @param ... Unused.
#' @return A p by p covariance matrix.
#' @export
vcov.sikdr <- function(object, ...) {
  if (is.null(object$V)) stop("refit with covariance = TRUE")
  S <- object$frame %*% object$V %*% t(object$frame) / object$n
  dimnames(S) <- list(names(object$beta), names(object$beta))
  S
}

#' @export
print.sikdr <- function(x, digits = 4, ...) {
  cat("Single-Index Kernel Density Regression (", x$tuning, " tuning)\n", sep = "")
  cat("n =", x$n, " p =", length(x$beta), "\n")
  cat("Direction:\n")
  print(round(x$beta, digits))
  cat("Smoothing tuple (g, h, lambda, rho):", format(signif(x$tuple, 4)), "\n")
  if (!is.null(x$inference_ok)) {
    cat("Inference conditions:", if (x$inference_ok) "satisfied" else
      paste("not satisfied:", paste(x$inference_reasons, collapse = ", ")), "\n")
  }
  invisible(x)
}
