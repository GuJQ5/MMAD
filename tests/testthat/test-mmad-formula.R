## End-to-end Phase 5 tests: pass a formula straight to mmad() and verify
## the same fit comes out as when we hand-construct the mmad_expr.

test_that("mmad() accepts a formula and matches direct mmad_expr construction (multinomial)", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr_direct <- 12 * log(0.5 * v1 + 0.5 * v2) +
                 15 * log((2/3) * v1 + (1/3) * v2) +
                 9  * log((1/3) * v1 + (2/3) * v2) -
                 6 * v1 - 6 * v2

  fit_direct <- mmad(expr_direct, init = c(4, 2),
                     tol = 1e-6, max_iter = 2000)

  fit_form <- mmad(
    ~ 12 * log(0.5 * theta[1] + 0.5 * theta[2]) +
      15 * log((2/3) * theta[1] + (1/3) * theta[2]) +
       9 * log((1/3) * theta[1] + (2/3) * theta[2]) -
       6 * theta[1] - 6 * theta[2],
    init = c(4, 2), tol = 1e-6, max_iter = 2000)

  expect_true(fit_form$converged)
  expect_equal(fit_form$estimate, fit_direct$estimate, tolerance = 1e-5)
})

test_that("named-parameter formula yields a named estimate vector", {
  fit <- mmad(
    ~ 12 * log(0.5 * a + 0.5 * b) +
      15 * log((2/3) * a + (1/3) * b) +
       9 * log((1/3) * a + (2/3) * b) -
       6 * a - 6 * b,
    init = c(a = 4, b = 2), tol = 1e-6, max_iter = 2000)

  expect_true(fit$converged)
  expect_equal(names(fit$estimate), c("a", "b"))
})

test_that("data argument supplies external values referenced in the formula", {
  # Toy linear model: maximize -0.5 * sum((y - x * theta)^2). Pre-aggregated
  # to: -0.5 * (sum(y^2) - 2 * sum(x*y) * theta + sum(x^2) * theta^2).
  x <- c(1, 2, 3); y <- c(2, 3.9, 6.1)
  fit <- mmad(
    ~ -0.5 * (sum(y^2) - 2 * sum(x * y) * theta[1] + sum(x^2) * theta[1]^2),
    init = c(theta = 0),
    data = list(x = x, y = y),
    tol = 1e-8)

  expect_true(fit$converged)
  beta_hat_ols <- sum(x * y) / sum(x^2)        # closed form
  expect_equal(unname(fit$estimate), beta_hat_ols, tolerance = 1e-6)
})

test_that("Poisson regression via sum() and X %*% theta matches glm()", {
  set.seed(42)
  n <- 100; p <- 3
  X <- matrix(rnorm(n * p, sd = 0.5), nrow = n)
  beta_true <- c(0.3, -0.2, 0.5)
  y <- rpois(n, exp(X %*% beta_true))

  fit <- mmad(
    ~ sum(y * (X %*% theta) - exp(X %*% theta)),
    init = rep(0, p),
    data = list(X = X, y = y),
    tol = 1e-6, max_iter = 500)

  expect_true(fit$converged)

  glm_fit <- glm(y ~ X - 1, family = poisson)
  expect_equal(unname(fit$estimate), unname(coef(glm_fit)), tolerance = 1e-3)
})
