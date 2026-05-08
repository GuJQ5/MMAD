## Phase 3 tests for minorize_at(). The contract has four pieces:
##   - tangency at theta_0: surrogate value AND gradient match the target;
##   - lower bound at perturbed theta: surrogate(theta) <= target(theta);
##   - separability is maximised: univariate sub-expressions are extracted
##     into per_coord; multivariate Jensen-able pieces become per_coord
##     Jensen surrogates; multivariate convex pieces become per_coord
##     hyperplanes; remaining multi-coord residue lands in non_separable;
##   - graceful degradation: when nothing further can be minorized, the
##     surrogate is still returned (with non_separable populated) so that
##     value/gradient/Hessian computations remain available.

surrogate_value <- function(surr, theta) {
  evaluate_expr(surr$expr, theta)$value
}
surrogate_gradient <- function(surr, theta) {
  evaluate_expr(surr$expr, theta)$gradient
}

# ---- Affine target: surrogate equals target everywhere -------------------

test_that("affine target produces an identity surrogate", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr    <- 3 * v1 - 2 * v2 + 5
  theta_0 <- c(1, 1)
  surr    <- minorize_at(expr, theta_0)

  expect_true(surr$is_separable)
  for (theta in list(c(1, 1), c(0.5, 2), c(-1, 3))) {
    expect_equal(surrogate_value(surr, theta),
                 evaluate_expr(expr, theta)$value,
                 tolerance = 1e-12)
  }
})

# ---- Univariate extraction -----------------------------------------------

test_that("univariate non-affine sub-expressions are extracted as-is", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- exp(v1) + log(v2) + 3 * v1 - v2 + 7
  theta_0 <- c(0.5, 2)
  surr <- minorize_at(expr, theta_0)

  expect_true(surr$is_separable)
  expect_null(surr$non_separable)
  expect_false(is.null(surr$per_coord[[1]]))
  expect_false(is.null(surr$per_coord[[2]]))

  # Tangency and exact equality everywhere (since the per_coord pieces
  # are the original 1-d functions: no minorization was applied).
  for (theta in list(c(0.5, 2), c(1, 1.5), c(0.2, 3))) {
    expect_equal(surrogate_value(surr, theta),
                 evaluate_expr(expr, theta)$value,
                 tolerance = 1e-10)
  }
})

test_that("non-DCP but per-coordinate-separable expression is fully separated", {
  # exp(v1) - exp(v2): convex - convex = unknown curvature, but each piece
  # depends on a single theta -- both are extracted into per_coord.
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- exp(v1) - exp(v2)
  theta_0 <- c(0.3, 0.4)
  surr <- minorize_at(expr, theta_0)

  expect_true(surr$is_separable)
  expect_null(surr$non_separable)

  # Surrogate equals target exactly (no Jensen / hyperplane was applied).
  for (theta in list(theta_0, c(0.6, 0.1), c(-0.2, 0.5))) {
    expect_equal(surrogate_value(surr, theta),
                 evaluate_expr(expr, theta)$value,
                 tolerance = 1e-10)
  }
})

# ---- Convex multivariate: tangent + lower bound --------------------------

test_that("convex multivariate sub-expression is minorized by its supporting hyperplane", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- exp(v1 + 0.5 * v2)              # multivariate -> hyperplane
  theta_0 <- c(0.3, -0.2)
  surr <- minorize_at(expr, theta_0)

  expect_true(surr$is_separable)

  ev0 <- evaluate_expr(expr, theta_0)
  expect_equal(surrogate_value(surr, theta_0), ev0$value, tolerance = 1e-10)
  expect_equal(surrogate_gradient(surr, theta_0), ev0$gradient,
               tolerance = 1e-10)

  set.seed(1)
  for (k in 1:20) {
    theta <- theta_0 + rnorm(2, sd = 0.1)
    expect_lte(surrogate_value(surr, theta),
               evaluate_expr(expr, theta)$value + 1e-9)
  }
})

# ---- Concave multivariate via Jensen: tangent + lower bound --------------

test_that("concave-of-affine multivariate is minorized by Jensen", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr    <- log(0.5 * v1 + 0.5 * v2)
  theta_0 <- c(2, 1)
  surr    <- minorize_at(expr, theta_0)

  expect_true(surr$is_separable)

  ev0 <- evaluate_expr(expr, theta_0)
  expect_equal(surrogate_value(surr, theta_0),    ev0$value,    tolerance = 1e-10)
  expect_equal(surrogate_gradient(surr, theta_0), ev0$gradient, tolerance = 1e-10)

  set.seed(2)
  for (k in 1:20) {
    theta <- pmax(theta_0 + rnorm(2, sd = 0.2), 0.05)
    expect_lte(surrogate_value(surr, theta),
               evaluate_expr(expr, theta)$value + 1e-9)
  }
})

# ---- Canonical multinomial objective: tangent + LB ----------------------

test_that("multinomial objective: tangent at theta_0 and lower-bounds at perturbed thetas", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- 12 * log(0.5 * v1 + 0.5 * v2) +
          15 * log((2/3) * v1 + (1/3) * v2) +
          9  * log((1/3) * v1 + (2/3) * v2) -
          6 * v1 - 6 * v2
  theta_0 <- c(4, 2)
  surr <- minorize_at(expr, theta_0)

  expect_true(surr$is_separable)

  ev0 <- evaluate_expr(expr, theta_0)
  expect_equal(surrogate_value(surr, theta_0),    ev0$value,    tolerance = 1e-10)
  expect_equal(surrogate_gradient(surr, theta_0), ev0$gradient, tolerance = 1e-10)

  set.seed(3)
  for (k in 1:30) {
    theta <- pmax(theta_0 + rnorm(2, sd = 0.3), 0.05)
    expect_lte(surrogate_value(surr, theta),
               evaluate_expr(expr, theta)$value + 1e-9)
  }
})

# ---- Soft failure: log(1 - positive_affine) ------------------------------

test_that("log(1 - positive_affine) returns a usable surrogate via the non-separable bucket", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- log(1 - (1/3) * v1 - (1/3) * v2)
  theta_0 <- c(0.5, 0.5)
  surr <- minorize_at(expr, theta_0)

  expect_false(surr$is_separable)
  expect_false(is.null(surr$non_separable))

  # The full surrogate equals the target at theta_0 (and everywhere,
  # since the only piece was added to non_separable as-is).
  ev0 <- evaluate_expr(expr, theta_0)
  expect_equal(surrogate_value(surr, theta_0),    ev0$value,    tolerance = 1e-10)
  expect_equal(surrogate_gradient(surr, theta_0), ev0$gradient, tolerance = 1e-10)

  for (theta in list(c(0.5, 0.5), c(0.6, 0.4), c(0.2, 0.7))) {
    expect_equal(surrogate_value(surr, theta),
                 evaluate_expr(expr, theta)$value,
                 tolerance = 1e-10)
  }
})

# ---- Concave-of-concave: handled by nested Jensen ------------------------

test_that("concave-of-concave is minorized via the nested Jensen trick", {
  # log(log(v1) + log(v2) + 5): the inner is concave (sum of concaves
  # plus a positive constant) and the outer log is concave + nondecreasing,
  # so the nested-Jensen path should produce a fully separable surrogate.
  v_pos <- mmad_var(1, sign = "positive")
  v2    <- mmad_var(2, sign = "positive")
  expr  <- log(log(v_pos) + log(v2) + 5)
  theta_0 <- c(2, 2)
  surr <- minorize_at(expr, theta_0)

  expect_true(surr$is_separable)
  expect_null(surr$non_separable)
  expect_false(is.null(surr$per_coord[[1]]))
  expect_false(is.null(surr$per_coord[[2]]))

  # Tangency at theta_0: value AND gradient match the target.
  ev0 <- evaluate_expr(expr, theta_0)
  expect_equal(surrogate_value(surr, theta_0),    ev0$value,    tolerance = 1e-10)
  expect_equal(surrogate_gradient(surr, theta_0), ev0$gradient, tolerance = 1e-10)

  # Lower bound at perturbed thetas (in the surrogate's domain, which is
  # tighter than the original's: each per_coord requires log(theta_j) > 0
  # i.e. theta_j > 1).
  set.seed(11)
  for (k in 1:20) {
    theta <- pmax(theta_0 + rnorm(2, sd = 0.2), 1.05)
    expect_lte(surrogate_value(surr, theta),
               evaluate_expr(expr, theta)$value + 1e-9)
  }
})

# ---- Mixed: separable + non-separable in one objective -------------------

test_that("an objective with both separable and non-separable parts yields a usable mixed surrogate", {
  # exp(v1) (univariate, separable) + log(1 - 0.5*v1 - 0.5*v2) (non-separable).
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- exp(v1) + log(1 - 0.5 * v1 - 0.5 * v2)
  theta_0 <- c(0.3, 0.3)
  surr <- minorize_at(expr, theta_0)

  expect_false(surr$is_separable)
  # exp(v1) extracted into per_coord[[1]] (only depends on theta_1):
  expect_false(is.null(surr$per_coord[[1]]))
  # The non-separable bucket carries the log(1 - ...) piece.
  expect_false(is.null(surr$non_separable))

  ev0 <- evaluate_expr(expr, theta_0)
  expect_equal(surrogate_value(surr, theta_0),    ev0$value,    tolerance = 1e-10)
  expect_equal(surrogate_gradient(surr, theta_0), ev0$gradient, tolerance = 1e-10)
})
