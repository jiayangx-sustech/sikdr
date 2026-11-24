test_that("grids follow the paper's rates", {
  g <- sikdr_grid(600)
  expect_length(g$g, 4)
  expect_length(g$h, 5)
  expect_equal(g$rho, 1 / 600)
  expect_equal(sikdr_grid(600, "reference")$lambda, 0.01 * 600^(-7 / 10))
})

test_that("fit, prediction, covariance and test work end to end", {
  set.seed(3)
  n <- 150
  x <- matrix(rnorm(n * 3), n, 3)
  e <- ifelse(runif(n) < 0.7, rnorm(n, -0.9, 0.5), rnorm(n, 2.1, 0.3))
  b0 <- c(1, 1, 0) / sqrt(2)
  y <- sin(2 * drop(x %*% b0)) + e
  fit <- sikdr(x, y)
  expect_s3_class(fit, "sikdr")
  expect_equal(sum(fit$beta^2), 1, tolerance = 1e-12)
  expect_lt(abs(abs(sum(fit$beta * b0)) - 1), 0.05)
  expect_length(predict(fit, x[1:5, ]), 5)
  S <- vcov(fit)
  expect_equal(dim(S), c(3, 3))
  expect_lt(max(abs(S %*% fit$beta)), 1e-10)
  tt <- sikdr_test(fit, beta = b0)
  expect_true(tt$p.value > 0 && tt$p.value <= 1)
  expect_output(print(fit), "Direction")
})

test_that("screening keeps the requested number of predictors", {
  set.seed(4)
  n <- 120
  x <- matrix(rnorm(n * 30), n, 30)
  y <- sin(drop(x[, 1:2] %*% c(1, 1)) / sqrt(2)) + rexp(n) - 1
  s <- sikdr_screen(x, y, size = 4, tuning = "reference", covariance = FALSE)
  expect_length(s$selected, 4)
  expect_length(s$beta, 30)
})
