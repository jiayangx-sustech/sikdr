# Projected Armijo ascent on the working chart z (||z|| <= R) around the Fourier direction.

chart_gamma <- function(z, g0, U) (g0 + drop(U %*% z)) / sqrt(1 + sum(z * z))

chart_jacobian <- function(z, g0, U) {
  s <- 1 + sum(z * z)
  U / sqrt(s) - tcrossprod(g0 + drop(U %*% z), z) / s^1.5
}

working_value <- function(z, Z, y, g0, U, tup) {
  kdr_value(drop(Z %*% chart_gamma(z, g0, U)), y, tup[1], tup[2], tup[3], tup[4])
}

working_value_grad <- function(z, Z, y, g0, U, tup) {
  sc <- kdr_score_coef(drop(Z %*% chart_gamma(z, g0, U)), y, tup[1], tup[2], tup[3], tup[4])
  grad_gamma <- drop(crossprod(Z, sc$c)) / nrow(Z)
  list(F = sc$L, grad = drop(crossprod(chart_jacobian(z, g0, U), grad_gamma)))
}

project_ball <- function(x, R) {
  nx <- sqrt(sum(x * x))
  if (nx <= R) x else x * (R / nx)
}

projected_armijo <- function(Z, y, g0, U, tup) {
  n <- nrow(Z)
  C <- ascent_constants(n)
  z <- numeric(ncol(U))
  vg <- working_value_grad(z, Z, y, g0, U, tup)
  Fz <- vg$F
  gz <- vg$grad
  it <- 0L
  status <- "max_iter"
  if (!is.finite(Fz) || !all(is.finite(gz))) {
    return(list(z = z, F = Fz, grad = gz, iterations = 0L, status = "nonfinite", R = C$R,
                interior = TRUE))
  }
  while (it < C$max_iter) {
    if (sqrt(sum(gz * gz)) <= C$grad_tol) {
      status <- "grad_tol"
      break
    }
    it <- it + 1L
    step <- C$step0
    accepted <- FALSE
    for (trial in seq_len(C$max_backtracks)) {
      cand <- z + step * gz
      znew <- project_ball(cand, C$R)
      dz <- znew - z
      Fn <- working_value(znew, Z, y, g0, U, tup)
      if (is.finite(Fn) && Fn >= Fz + 0.5 * sum(gz * dz)) {
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    if (!accepted) {
      status <- "armijo_exhausted"
      break
    }
    if (sqrt(sum(dz * dz)) == 0) {
      status <- "projected_stationary"
      break
    }
    z <- znew
    vg <- working_value_grad(z, Z, y, g0, U, tup)
    Fz <- vg$F
    gz <- vg$grad
    if (!is.finite(Fz) || !all(is.finite(gz))) {
      status <- "nonfinite"
      break
    }
  }
  if (status == "max_iter" && sqrt(sum(gz * gz)) <= C$grad_tol) status <- "grad_tol"
  list(z = z, F = Fz, grad = gz, iterations = it, status = status, R = C$R,
       interior = sqrt(sum(z * z)) < C$R, gamma = chart_gamma(z, g0, U))
}
