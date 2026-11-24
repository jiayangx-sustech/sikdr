# Agreement with the Python implementation used for the paper's numerical work
# (reference values were produced in float64 on the CPU from the same data).

check_case <- function(file) {
  ref <- readRDS(test_path(file))
  X <- ref$X
  Y <- ref$Y
  n <- nrow(X)
  std <- sikdr:::standardize_xy(X, Y)
  fo <- sikdr:::fourier_direction(std$Z, std$y)
  expect_equal(fo$gamma, ref$fourier_gamma, tolerance = 1e-10)
  T0 <- drop(std$Z %*% ref$fourier_gamma)
  sc <- sikdr:::kdr_score_coef(T0, std$y, 0.3, 0.2, 0.001, 1 / n)
  expect_equal(sc$L, ref$L, tolerance = 1e-12)
  expect_equal(sc$c, ref$c, tolerance = 1e-10)
  expect_equal(sikdr:::loo_mse(T0, std$y, 0.3, 0.001), ref$loo, tolerance = 1e-12)
  for (m in c("tangent", "complex")) {
    ij <- sikdr:::ij_covariance(X, Y, ref$beta0, ref$ij_F, ref$ij_tuple, method = m)
    expect_equal(ij$S, ref$ij_S, tolerance = 1e-10)
    expect_equal(ij$J, ref$ij_J, tolerance = 1e-10)
    expect_equal(ij$B, ref$ij_B, tolerance = 1e-10)
  }
  fit <- sikdr(X, Y, tuning = "trace")
  expect_equal(unname(fit$tuple), ref$kdr_tuple, tolerance = 1e-12)
  expect_identical(as.integer(fit$diagnostics$selected_index), as.integer(ref$kdr_selected))
  ok <- is.finite(ref$kdr_traces)
  expect_equal(fit$diagnostics$traces[ok], ref$kdr_traces[ok], tolerance = 1e-8)
  expect_equal(fit$beta_internal, ref$kdr_beta, tolerance = 1e-8)
  expect_true(fit$inference_ok)
  expect_equal(fit$V, ref$kdr_V, tolerance = 1e-8)
  fr <- sikdr(X, Y, tuning = "reference", covariance = FALSE)
  expect_equal(fr$beta_internal, ref$kdrr_beta, tolerance = 1e-8)
  ts <- sikdr_test(X, Y, rbind(ref$beta0, ref$kdr_beta))
  expect_equal(ts$F, ref$cs_F, tolerance = 1e-8)
  expect_equal(dc_scores(X, Y), ref$dcsis, tolerance = 1e-12)
  half <- seq_len(n / 2)
  pr <- sikdr:::nw_index_predict(X[half, ], Y[half], X[-half, ], ref$kdr_beta)
  expect_equal(pr$bandwidth, ref$pred_bw, tolerance = 1e-12)
  expect_equal(mean((Y[-half] - pr$pred)^2), ref$pred_err, tolerance = 1e-10)
}

test_that("SI-KDR matches the reference implementation (n = 200, p = 3)", {
  check_case("ref_c1.rds")
})

test_that("SI-KDR matches the reference implementation (n = 300, p = 5)", {
  skip_on_cran()
  check_case("ref_c2.rds")
})
