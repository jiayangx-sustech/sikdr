# Protected double leave-one-out criterion of SI-KDR and its analytic score coefficient.
# All kernels are normalized Gaussian kernels; diagonals are removed exactly (set to zero).

# Exponents below -700 are set to -Inf: those kernel values (< 1e-304) are negligible against the
# ridge and the floor, and skipping them avoids very slow subnormal arithmetic. The comparison
# uses the real part, so the kernel stays usable with a complex step.
safe_exp <- function(e) {
  low <- Re(e) < -700
  if (any(low)) e[low] <- -Inf
  exp(e)
}

gauss_kernel <- function(u, b) safe_exp(-(u * u) / (2 * b * b)) / (sqrt(2 * pi) * b)

# Leave-one-out Nadaraya-Watson fit with ridge (n - 1) * lambda in the denominator.
nw_protected <- function(tt, y, g, lambda) {
  n <- length(tt)
  D <- outer(tt, tt, "-")
  Kg <- gauss_kernel(D, g)
  diag(Kg) <- 0
  den <- rowSums(Kg) + (n - 1) * lambda
  m <- drop(Kg %*% y) / den
  list(m = m, D = D, Kg = Kg, den = den)
}

# Criterion L_n = n^{-1} sum log(d_i + rho) together with the pieces needed for the score.
kdr_parts <- function(tt, y, g, h, lambda, rho) {
  n <- length(tt)
  nw <- nw_protected(tt, y, g, lambda)
  r <- y - nw$m
  Dr <- outer(r, r, "-")
  Kh <- gauss_kernel(Dr, h)
  diag(Kh) <- 0
  d <- rowSums(Kh) / (n - 1)
  c(nw, list(r = r, Dr = Dr, Kh = Kh, d = d, L = mean(log(d + rho))))
}

kdr_value <- function(tt, y, g, h, lambda, rho) kdr_parts(tt, y, g, h, lambda, rho)$L

loo_mse <- function(tt, y, g, lambda) mean((y - nw_protected(tt, y, g, lambda)$m)^2)

# Score coefficient c with n dL_n/dT_i = c_i (paper, Section 4.1):
# c_i = -a_i A_i + sum_j a_j Q_ji (y_i - m_j).
kdr_score_coef <- function(tt, y, g, h, lambda, rho) {
  P <- kdr_parts(tt, y, g, h, lambda, rho)
  n <- length(tt)
  Q <- (-P$D / (g * g)) * P$Kg / P$den
  A <- drop(Q %*% y) - P$m * rowSums(Q)
  inv <- 1 / (P$d + rho)
  Khp <- (-P$Dr / (h * h)) * P$Kh
  a <- (rowSums(Khp) * inv + drop(Khp %*% inv)) / (n - 1)
  QTa <- drop(crossprod(Q, a))
  QTam <- drop(crossprod(Q, a * P$m))
  list(c = -a * A + y * QTa - QTam, L = P$L, a = a)
}

# g by protected leave-one-out squared error, then h maximizing L_n at that g (first index on ties).
select_g_then_h <- function(tt, y, gs, hs, lambda, rho) {
  mse <- vapply(gs, function(g) loo_mse(tt, y, g, lambda), numeric(1))
  ig <- which.min(mse)
  Ls <- vapply(hs, function(h) kdr_value(tt, y, gs[ig], h, lambda, rho), numeric(1))
  ih <- which.max(Ls)
  list(ig = ig, ih = ih, g = gs[ig], h = hs[ih], mse = mse, L = Ls)
}
