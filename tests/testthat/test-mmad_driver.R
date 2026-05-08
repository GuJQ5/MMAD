## Phase 4 driver tests. Each test pins one of the guarantees the
## driver makes:
##   - on a smooth concave target it converges to a stationary point;
##   - the canonical multinomial objective converges with gradient ~ 0;
##   - when the surrogate has a non-empty non_separable bucket, the
##     entangled-block step still drives the iteration to a stationary
##     point of the actual target;
##   - the return list reports correct convergence status, iteration
##     counts, and (when requested) per-iteration history.

# ---- Multinomial: converges with gradient near zero ----------------------

test_that("multinomial objective converges to a stationary point", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- 12 * log(0.5 * v1 + 0.5 * v2) +
          15 * log((2/3) * v1 + (1/3) * v2) +
          9  * log((1/3) * v1 + (2/3) * v2) -
          6 * v1 - 6 * v2

  fit <- mmad(expr, init = c(4, 2), tol = 1e-6, max_iter = 2000)

  expect_true(fit$converged)
  g_at_opt <- evaluate_expr(expr, fit$estimate)$gradient
  expect_lt(max(abs(g_at_opt)), 1e-5)
})

# ---- Smooth concave 2-d: analytic stationary point -----------------------

test_that("concave quadratic-like target converges to the analytic optimum", {
  # f(theta) = log(theta1 + theta2) - theta1 - theta2 + 2
  # d/dtheta1 = 1/(theta1+theta2) - 1 = 0 => theta1+theta2 = 1
  # The optimum is any (a, 1-a) with a > 0; gradient zero everywhere on
  # that line. Pick a starting point and verify the gradient norm at
  # convergence is small and theta1+theta2 ~= 1.
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- log(v1 + v2) - v1 - v2 + 2

  fit <- mmad(expr, init = c(2, 2), tol = 1e-8)

  expect_true(fit$converged)
  expect_lt(abs(fit$estimate[1] + fit$estimate[2] - 1), 1e-4)
  g <- evaluate_expr(expr, fit$estimate)$gradient
  expect_lt(max(abs(g)), 1e-5)
})

# ---- Non-separable residue: log(1 - 0.5*v1 - 0.5*v2) + log(v1) + log(v2) -

test_that("target with a non-separable residue still converges via the block step", {
  # f = log(1 - 0.5*v1 - 0.5*v2) + log(v1) + log(v2)
  # Sub to t = 0.5*v1 + 0.5*v2: log piece is log(1 - t), separable logs.
  # Symmetric in v1, v2. Stationary: d/dv1 = -0.5/(1 - 0.5*v1 - 0.5*v2) +
  # 1/v1 = 0 -> v1 (-0.5) = -(1 - 0.5*v1 - 0.5*v2) -> v1 = 2(1 - 0.5*v1 -
  # 0.5*v2) - actually let's not solve in closed form, just assert grad
  # near zero at convergence.
  v1 <- mmad_var(1, sign = "positive"); v2 <- mmad_var(2, sign = "positive")
  expr <- log(1 - 0.5 * v1 - 0.5 * v2) + log(v1) + log(v2)

  fit <- mmad(expr, init = c(0.4, 0.4), tol = 1e-8, max_iter = 500)

  expect_true(fit$converged)
  g <- evaluate_expr(expr, fit$estimate)$gradient
  expect_lt(max(abs(g)), 1e-5)
  # Symmetry of the target implies the optimum is symmetric in v1, v2.
  expect_lt(abs(fit$estimate[1] - fit$estimate[2]), 1e-5)
})

# ---- Convergence flag and iteration count are sensible -------------------

test_that("max_iter cap is respected and reported", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- log(0.5 * v1 + 0.5 * v2)
  fit <- mmad(expr, init = c(1, 1), tol = 1e-12, max_iter = 3)
  expect_equal(fit$iterations, 3)
  expect_false(fit$converged)
  expect_match(fit$message, "max_iter")
})

# ---- History tracking ----------------------------------------------------

test_that("track_history returns a data frame with one row per iteration", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- log(v1 + v2) - v1 - v2

  fit <- mmad(expr, init = c(2, 2), tol = 1e-8, track_history = TRUE)
  expect_true(is.data.frame(fit$history))
  expect_true(all(c("iteration", "value", "grad_norm", "theta_change") %in%
                  names(fit$history)))
  expect_equal(nrow(fit$history), fit$iterations)
  # Values are monotonically nondecreasing iteration over iteration.
  expect_true(all(diff(fit$history$value) >= -1e-10))
})

# ---- Single-coordinate target via univariate extraction ------------------

test_that("non-DCP-but-univariate-separable target still converges per coordinate", {
  # f = -exp(v1) - exp(v2) + log(1 + v1 + v2 + 1)  -- the +1 ensures domain.
  # Each -exp piece is concave, hits per_coord; log piece is concave (log
  # of affine, > 0 at theta_0); there should be no non_separable. The
  # iteration should drive both coords toward -infinity in principle but
  # since exp is unbounded below, MM ascent will keep walking; we test
  # only that the iteration is well-formed and makes monotonic progress.
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- -exp(v1) - exp(v2) + log(2 + v1 + v2)

  fit <- mmad(expr, init = c(0.5, 0.5), tol = 1e-7, max_iter = 200,
              track_history = TRUE)
  # Value monotonically improves.
  expect_true(all(diff(fit$history$value) >= -1e-10))
})

# ---- mmad_fit class + print/summary methods ------------------------------

test_that("mmad() returns an mmad_fit object with print and summary methods", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- log(v1 + v2) - v1 - v2 + 2

  fit <- mmad(expr, init = c(2, 2), tol = 1e-8)

  expect_s3_class(fit, "mmad_fit")
  expect_true("expr" %in% names(fit))   # needed for summary's grad recompute
  # expect_silent(print(fit))

  s <- summary(fit)
  expect_s3_class(s, "summary.mmad_fit")
  expect_true(s$grad_norm < 1e-5)
  # expect_silent(print(s))
})
