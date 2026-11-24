# Case-weight (infinitesimal jackknife) covariance of SI-KDR in a raw tangent chart.
#
# The weighted criterion of the paper (Section 5.1) is
#   L_w = sum_i p_i log(d_i + rho),  p = w / sum(w),
# with the weighted index t_i = (x_i - P_w x)/s_T(w), x = X beta(v), the weighted response,
# the weighted ridge smoother and the weighted density estimate.  cw_pass() evaluates L_1,
# its gradient in v and its gradient with respect to the weights p (reverse mode, written out
# by hand).  Every operation is analytic in v, so a complex step in v returns the exact
# Hessian in v and the exact mixed derivatives eta_i = n d/dw_i grad_v L_w.

cw_pass <- function(v, X, Y, beta, Fr, g, h, lambda, rho) {
  n <- nrow(X)
  p <- rep(1 / n, n)
  nv <- 1 + sum(v * v)
  b <- (beta + drop(Fr %*% v)) / sqrt(nv)
  x <- drop(X %*% b)
  xbar <- sum(p * x)
  xc <- x - xbar
  s2 <- sum(p * xc * xc)
  sT <- sqrt(s2)
  tt <- xc / sT
  ybar <- sum(p * Y)
  yc <- Y - ybar
  sy2 <- sum(p * yc * yc)
  sy <- sqrt(sy2)
  yy <- yc / sy
  # weighted ridge smoother
  D <- outer(tt, tt, "-")
  K <- gauss_kernel(D, g)
  diag(K) <- 0
  Kp <- (-D / (g * g)) * K
  num <- drop(K %*% (p * yy))
  den <- drop(K %*% p) + (1 - p) * lambda
  m <- num / den
  r <- yy - m
  # weighted density estimate at the residuals
  Dr <- outer(r, r, "-")
  H <- gauss_kernel(Dr, h)
  diag(H) <- 0
  Hp <- (-Dr / (h * h)) * H
  dd <- drop(H %*% p) / (1 - p)
  ld <- log(dd + rho)
  L <- sum(p * ld)
  # reverse pass
  u <- p / (dd + rho)
  ut <- u / (1 - p)
  gp <- ld + drop(crossprod(H, ut)) + u * dd / (1 - p)
  ar <- ut * drop(Hp %*% p) + p * drop(Hp %*% ut)
  ay <- ar
  am <- -ar
  bn <- am / den
  bd <- -am * m / den
  KTbn <- drop(crossprod(K, bn))
  gp <- gp + yy * KTbn + drop(crossprod(K, bd)) - lambda * bd
  ay <- ay + p * KTbn
  at <- bn * drop(Kp %*% (p * yy)) + bd * drop(Kp %*% p) + p * yy * drop(Kp %*% bn) +
    p * drop(Kp %*% bd)
  sat <- sum(at)
  satt <- sum(at * tt)
  gp <- gp - (x / sT) * sat - (xc * xc / (2 * s2)) * satt
  gp <- gp - (Y / sy) * sum(ay) - (yc * yc / (2 * sy2)) * sum(ay * yy)
  # gradient in v at fixed weights
  ax <- at / sT - p * sat / sT - p * xc * satt / s2
  gb <- drop(crossprod(X, ax))
  Jb <- Fr / sqrt(nv) - tcrossprod(beta + drop(Fr %*% v), v) / nv^1.5
  list(L = L, grad_v = drop(crossprod(Jb, gb)), grad_p = gp)
}

# Real pass at v = 0 that keeps every intermediate needed by the tangent pass.
cw_base <- function(X, Y, beta, g, h, lambda, rho) {
  n <- nrow(X)
  p <- rep(1 / n, n)
  x <- drop(X %*% beta)
  xc <- x - sum(p * x)
  s2 <- sum(p * xc * xc)
  sT <- sqrt(s2)
  tt <- xc / sT
  yc <- Y - sum(p * Y)
  sy2 <- sum(p * yc * yc)
  sy <- sqrt(sy2)
  yy <- yc / sy
  D <- outer(tt, tt, "-")
  K <- gauss_kernel(D, g)
  diag(K) <- 0
  Kp <- (-D / (g * g)) * K
  num <- drop(K %*% (p * yy))
  den <- drop(K %*% p) + (1 - p) * lambda
  m <- num / den
  r <- yy - m
  Dr <- outer(r, r, "-")
  H <- gauss_kernel(Dr, h)
  diag(H) <- 0
  Hp <- (-Dr / (h * h)) * H
  dd <- drop(H %*% p) / (1 - p)
  ld <- log(dd + rho)
  u <- p / (dd + rho)
  ut <- u / (1 - p)
  am <- -(ut * drop(Hp %*% p) + p * drop(Hp %*% ut))
  bn <- am / den
  bd <- -am * m / den
  at <- bn * drop(Kp %*% (p * yy)) + bd * drop(Kp %*% p) + p * yy * drop(Kp %*% bn) +
    p * drop(Kp %*% bd)
  list(n = n, p = p, x = x, xc = xc, s2 = s2, sT = sT, tt = tt, Y = Y, yc = yc, sy2 = sy2, sy = sy,
       yy = yy, D = D, K = K, Kp = Kp, den = den, m = m, Dr = Dr, H = H, Hp = Hp, dd = dd, u = u,
       ut = ut, am = am, bn = bn, bd = bd, at = at, KTbn = drop(crossprod(K, bn)),
       sat = sum(at), satt = sum(at * tt), g = g, h = h, lambda = lambda, rho = rho,
       gb = drop(crossprod(X, at / sT - p * sum(at) / sT - p * xc * sum(at * tt) / s2)))
}

# Forward-mode tangent of cw_pass along v = e_k at v = 0 (same quantities as a complex step,
# in real arithmetic and reusing the kernels of the base pass).
cw_tangent <- function(B, X, beta, fk) {
  p <- B$p
  g <- B$g
  h <- B$h
  xd <- drop(X %*% fk)
  xcd <- xd - sum(p * xd)
  s2d <- 2 * sum(p * B$xc * xcd)
  sTd <- s2d / (2 * B$sT)
  td <- xcd / B$sT - B$tt * sTd / B$sT
  Dd <- outer(td, td, "-")
  Kd <- B$Kp * Dd
  Kpd <- (-Dd / (g * g)) * B$K + (-B$D / (g * g)) * Kd
  py <- p * B$yy
  numd <- drop(Kd %*% py)
  dend <- drop(Kd %*% p)
  md <- (numd - B$m * dend) / B$den
  rd <- -md
  Drd <- outer(rd, rd, "-")
  Hd <- B$Hp * Drd
  Hpd <- (-Drd / (h * h)) * B$H + (-B$Dr / (h * h)) * Hd
  ddd <- drop(Hd %*% p) / (1 - p)
  ldd <- ddd / (B$dd + B$rho)
  ud <- -B$u * ddd / (B$dd + B$rho)
  utd <- ud / (1 - p)
  gpd <- ldd + drop(crossprod(Hd, B$ut)) + drop(crossprod(B$H, utd)) +
    (ud * B$dd + B$u * ddd) / (1 - p)
  ard <- utd * drop(B$Hp %*% p) + B$ut * drop(Hpd %*% p) +
    p * (drop(Hpd %*% B$ut) + drop(B$Hp %*% utd))
  amd <- -ard
  bnd <- (amd - B$bn * dend) / B$den
  bdd <- -(amd * B$m + B$am * md - B$am * B$m * dend / B$den) / B$den
  KTbnd <- drop(crossprod(Kd, B$bn)) + drop(crossprod(B$K, bnd))
  gpd <- gpd + B$yy * KTbnd + drop(crossprod(Kd, B$bd)) + drop(crossprod(B$K, bdd)) -
    B$lambda * bdd
  ayd <- ard + p * KTbnd
  atd <- bnd * drop(B$Kp %*% py) + B$bn * drop(Kpd %*% py) + bdd * drop(B$Kp %*% p) +
    B$bd * drop(Kpd %*% p) + py * (drop(Kpd %*% B$bn) + drop(B$Kp %*% bnd)) +
    p * (drop(Kpd %*% B$bd) + drop(B$Kp %*% bdd))
  satd <- sum(atd)
  sattd <- sum(atd * B$tt + B$at * td)
  gpd <- gpd - (xd / B$sT - B$x * sTd / B$s2) * B$sat - (B$x / B$sT) * satd -
    (B$xc * xcd / B$s2 - B$xc^2 * s2d / (2 * B$s2^2)) * B$satt - (B$xc^2 / (2 * B$s2)) * sattd
  gpd <- gpd - (B$Y / B$sy) * sum(ayd) - (B$yc^2 / (2 * B$sy2)) * sum(ayd * B$yy)
  axd <- atd / B$sT - B$at * sTd / B$s2 - p * (satd / B$sT - B$sat * sTd / B$s2) -
    p * (xcd * B$satt / B$s2 + B$xc * sattd / B$s2 - B$xc * B$satt * s2d / B$s2^2)
  list(gbd = drop(crossprod(X, axd)), gpd = gpd)
}

# J = -Hess_v L_1, B = n^{-1} sum (eta - etabar)(eta - etabar)^T, V = J^{-1} B J^{-T}.
# method = "tangent" (default, real forward mode) or "complex" (complex step through cw_pass).
ij_covariance <- function(X, Y, beta, Fr, tup, method = c("tangent", "complex")) {
  method <- match.arg(method)
  n <- nrow(X)
  d <- ncol(Fr)
  Hs <- matrix(0, d, d)
  eta <- matrix(0, n, d)
  if (method == "tangent") {
    B0 <- cw_base(X, Y, beta, tup[1], tup[2], tup[3], tup[4])
    for (k in seq_len(d)) {
      td <- cw_tangent(B0, X, beta, Fr[, k])
      Hs[, k] <- drop(crossprod(Fr, td$gbd)) - replace(numeric(d), k, sum(beta * B0$gb))
      eta[, k] <- td$gpd - mean(td$gpd)
    }
    base <- list(grad_v = drop(crossprod(Fr, B0$gb)), L = sum(B0$p * log(B0$dd + B0$rho)))
  } else {
    hstep <- 1e-20
    base <- cw_pass(numeric(d), X, Y, beta, Fr, tup[1], tup[2], tup[3], tup[4])
    for (k in seq_len(d)) {
      v <- complex(real = numeric(d), imaginary = replace(numeric(d), k, hstep))
      res <- cw_pass(v, X, Y, beta, Fr, tup[1], tup[2], tup[3], tup[4])
      Hs[, k] <- Im(res$grad_v) / hstep
      de <- Im(res$grad_p) / hstep
      eta[, k] <- de - mean(de)
    }
  }
  J <- -(Hs + t(Hs)) / 2
  B <- crossprod(eta) / n
  out <- list(S = base$grad_v, L = base$L, J = J, B = B, eta = eta, finite = all(is.finite(c(J, B))),
              J_pd = FALSE, B_pd = FALSE, V = NULL, trace = NA_real_)
  if (out$finite) {
    out$J_pd <- min(eigen(J, symmetric = TRUE, only.values = TRUE)$values) > 0
    out$B_pd <- min(eigen((B + t(B)) / 2, symmetric = TRUE, only.values = TRUE)$values) > 0
    if (out$J_pd) {
      Ji <- solve(J)
      V <- Ji %*% B %*% t(Ji)
      out$V <- (V + t(V)) / 2
      out$trace <- sum(diag(out$V))
    }
  }
  out
}
