# sikdr: Single-Index Kernel Density Regression

R implementation of SI-KDR from

> Jiayang Xu, Yusen Jia, Yangzhou Chen and Xin Chen. *Finite-Smoothing Efficiency and
> Inference in Single-Index Models.* [arXiv:2611.22384](https://arxiv.org/abs/2611.22384).

In the single-index model `Y = m(beta' X) + e` with the error independent of the predictors,
methods based on the conditional mean (such as MAVE) ignore the error distribution. When the
error is skewed or multimodal, that distribution also carries information about `beta`.
SI-KDR refits the link at every candidate direction and scores the direction by the
leave-one-out kernel log-likelihood of the residuals.

The package provides

* `sikdr()`: the estimator, with bandwidths ranked by the trace of a case-weight covariance
  estimate (`tuning = "trace"`) or the reference rule (`tuning = "reference"`);
* `vcov()`: the case-weight covariance of the estimated direction;
* `sikdr_test()` and `sikdr_confset()`: a conditional score test that is exact under Gaussian
  predictors, and the confidence set obtained by inverting it;
* `sikdr_screen()` and `dc_scores()`: distance-correlation screening for more predictors than
  observations;
* `predict()`: kernel regression on the fitted index;
* `energy`: the UCI Energy Efficiency data used in the paper's benchmark.

The numerical results agree with the Python implementation used for the paper to about
1e-12 (see `tests/testthat/test-reference.R`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("jiayangx-sustech/sikdr")
```

## Example

```r
library(sikdr)
set.seed(1)
n <- 300
b0 <- c(1, 1, 0) / sqrt(2)
x <- matrix(rnorm(n * 3), n, 3)
e <- ifelse(runif(n) < 0.7, rnorm(n, -0.9, 0.5), rnorm(n, 2.1, 0.3))
y <- sin(2 * drop(x %*% b0)) + e

fit <- sikdr(x, y)
fit
vcov(fit)
sikdr_test(fit, beta = b0)
```

See `vignette("sikdr")` for confidence sets, screening and the Energy Efficiency benchmark.

## License

GPL (>= 3)
