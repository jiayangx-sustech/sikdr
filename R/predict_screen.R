#' Predict from a fitted SI-KDR model
#'
#' Gaussian Nadaraya--Watson regression of the training response on the fitted index, with
#' the bandwidth chosen by leave-one-out squared error over
#' \eqn{s_T n^{-1/5} \{0.25, 0.4, 0.6, 0.8, 1, 1.4, 2, 3\}} (\eqn{s_T} is the standard deviation
#' of the training index).
#'
#' @param object A \code{"sikdr"} object.
#' @param newdata Matrix of new predictors; the training predictors when missing.
#' @param ... Unused.
#' @return Predicted responses on the original scale.
#' @export
predict.sikdr <- function(object, newdata, ...) {
  if (missing(newdata)) newdata <- object$x
  newdata <- as.matrix(newdata)
  nw_index_predict(object$x, object$y, newdata, object$beta)$pred
}

#' Kernel regression on a single index
#'
#' Predicts \code{xte} from \code{(xtr, ytr)} by Gaussian Nadaraya--Watson regression on the
#' index \eqn{x^Teta}, with the bandwidth chosen by leave-one-out squared error as in
#' \code{\link{predict.sikdr}}. Useful for comparing directions from different estimators
#' with one prediction rule.
#'
#' @param xtr Training predictors.
#' @param ytr Training response.
#' @param xte Predictors to predict.
#' @param beta Direction defining the index.
#' @return A list with the predictions \code{pred} and the chosen \code{bandwidth}.
#' @export
index_predict <- function(xtr, ytr, xte, beta) {
  nw_index_predict(as.matrix(xtr), as.numeric(ytr), as.matrix(xte), as.numeric(beta))
}

nw_index_predict <- function(xtr, ytr, xte, beta) {
  n <- nrow(xtr)
  ttr <- drop(xtr %*% beta)
  tte <- drop(xte %*% beta)
  bws <- stats::sd(ttr) * n^(-0.2) * c(0.25, 0.4, 0.6, 0.8, 1, 1.4, 2, 3)
  D <- outer(ttr, ttr, "-")
  loo <- vapply(bws, function(b) {
    lw <- -(D * D) / (2 * b * b)
    diag(lw) <- -Inf
    mean((ytr - softmax_rows(lw) %*% ytr)^2)
  }, numeric(1))
  k <- which.min(loo)
  Dte <- outer(tte, ttr, "-")
  list(pred = drop(softmax_rows(-(Dte * Dte) / (2 * bws[k]^2)) %*% ytr), bandwidth = bws[k])
}

softmax_rows <- function(lw) {
  mx <- apply(lw, 1, max)
  e <- safe_exp(lw - mx)
  e / rowSums(e)
}

#' Distance-correlation screening
#'
#' Squared distance correlation (biased V-statistic) between each predictor and the
#' response, as used for independent screening before fitting SI-KDR in high dimensions.
#'
#' @param x Numeric predictor matrix.
#' @param y Response vector.
#' @return A numeric vector of scores, one per column of \code{x}.
#' @export
dc_scores <- function(x, y) {
  x <- as.matrix(x)
  dcen <- function(D) {
    D - rowMeans(D) - rep(colMeans(D), each = nrow(D)) + mean(D)
  }
  B <- dcen(abs(outer(y, y, "-")))
  bb <- mean(B * B)
  vapply(seq_len(ncol(x)), function(j) {
    A <- dcen(abs(outer(x[, j], x[, j], "-")))
    den <- sqrt(mean(A * A) * bb)
    if (den == 0) 0 else mean(A * B) / den
  }, numeric(1))
}

#' SI-KDR after independent screening
#'
#' Splits the sample, ranks the predictors by \code{\link{dc_scores}} on the screening part,
#' keeps the top \code{size} predictors and fits \code{\link{sikdr}} on the estimation part.
#' When the screen retains the support, the paper shows that the covariance comparison and
#' the exact test carry over to the retained model.
#'
#' @param x Numeric predictor matrix (may have more columns than rows).
#' @param y Response vector.
#' @param size Number of predictors to retain (at least 2).
#' @param screen_fraction Fraction of observations used for screening.
#' @param ... Further arguments passed to \code{\link{sikdr}}.
#' @return A list with the screened column indices, the fit on the retained columns and the
#'   direction embedded in all \code{ncol(x)} coordinates.
#' @export
sikdr_screen <- function(x, y, size, screen_fraction = 1 / 3, ...) {
  x <- as.matrix(x)
  n <- nrow(x)
  ns <- floor(screen_fraction * n)
  idx <- sample.int(n, ns)
  sc <- dc_scores(x[idx, , drop = FALSE], y[idx])
  keep <- sort(order(-sc, seq_along(sc))[seq_len(size)])
  fit <- sikdr(x[-idx, keep, drop = FALSE], y[-idx], ...)
  beta <- numeric(ncol(x))
  beta[keep] <- fit$beta
  list(selected = keep, scores = sc, fit = fit, beta = beta, screen_rows = idx)
}
