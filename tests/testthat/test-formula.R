## Phase 5 parser tests. Each test pins one piece of the documented
## grammar by parsing a formula and comparing the parsed tree with a
## hand-built mmad_expr that we already trust through Phase 1.

# Helper: build an evaluator that compares two mmad_expr trees by
# evaluating them at a few thetas and checking value and gradient agree.
expect_equivalent_expr <- function(parsed, hand_built, thetas) {
  for (theta in thetas) {
    a <- evaluate_expr(parsed,     theta)
    b <- evaluate_expr(hand_built, theta)
    expect_equal(a$value,    b$value,    tolerance = 1e-12)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-12)
  }
}

# ---- theta[i] indexing ---------------------------------------------------

test_that("theta[i] indexing parses to the right mmad_var", {
  parsed     <- as_mmad_expr(~ theta[1] + 2 * theta[2], init = c(0, 0))
  hand_built <- mmad_var(1) + 2 * mmad_var(2)
  expect_equivalent_expr(parsed, hand_built,
                         list(c(0, 0), c(1, 2), c(-1, 0.5)))
})

# ---- Named parameters ----------------------------------------------------

test_that("named parameters resolve via names(init)", {
  parsed <- as_mmad_expr(~ alpha + 2 * beta,
                         init = c(alpha = 0, beta = 0))
  hand_built <- mmad_var(1) + 2 * mmad_var(2)
  expect_equivalent_expr(parsed, hand_built,
                         list(c(0, 0), c(1, 2)))
})

test_that("indexed and named parameter syntax can be mixed", {
  parsed <- as_mmad_expr(~ alpha + theta[2],
                         init = c(alpha = 0, beta = 0))
  hand_built <- mmad_var(1) + mmad_var(2)
  expect_equivalent_expr(parsed, hand_built, list(c(1, 2), c(3, 1)))
})

# ---- Theta-free sub-expressions are evaluated eagerly --------------------

test_that("theta-free sub-expressions are evaluated and frozen as constants", {
  X <- c(1, 2, 3, 4)            # in this test's frame
  parsed <- as_mmad_expr(~ sum(X) * theta[1] - mean(X) * theta[2],
                         init = c(0, 0))
  hand_built <- 10 * mmad_var(1) - 2.5 * mmad_var(2)
  expect_equivalent_expr(parsed, hand_built, list(c(0, 0), c(1, 1)))
})

test_that("the data argument supplies values for theta-free symbols", {
  parsed <- as_mmad_expr(~ a * theta[1] + b * theta[2],
                         init = c(0, 0),
                         data = list(a = 3, b = -2))
  hand_built <- 3 * mmad_var(1) - 2 * mmad_var(2)
  expect_equivalent_expr(parsed, hand_built, list(c(0, 0), c(2, 5)))
})

# ---- log(1 - x) parses as expected, no special syntax needed -------------

test_that("log(1 - affine) parses to log of an additive structure", {
  parsed <- as_mmad_expr(~ log(1 - 0.5 * theta[1] - 0.5 * theta[2]),
                         init = c(0, 0))
  hand_built <- log(1 - 0.5 * mmad_var(1) - 0.5 * mmad_var(2))
  expect_equivalent_expr(parsed, hand_built,
                         list(c(0.2, 0.2), c(0.5, 0.4)))
})

# ---- Math primitives -----------------------------------------------------

test_that("log/exp/sqrt all parse correctly", {
  parsed <- as_mmad_expr(~ log(theta[1]) + exp(theta[2]) + sqrt(theta[1]),
                         init = c(1, 0))
  hand_built <- log(mmad_var(1)) + exp(mmad_var(2)) +
    sqrt(mmad_var(1))
  expect_equivalent_expr(parsed, hand_built,
                         list(c(2, 0.5), c(0.7, -0.3)))
})

# ---- Multinomial via formula (the canonical reviewer-level example) -----

test_that("multinomial objective parses and evaluates identically to direct construction", {
  parsed <- as_mmad_expr(
    ~ 12 * log(0.5 * theta[1] + 0.5 * theta[2]) +
      15 * log((2/3) * theta[1] + (1/3) * theta[2]) +
       9 * log((1/3) * theta[1] + (2/3) * theta[2]) -
       6 * theta[1] - 6 * theta[2],
    init = c(0, 0))

  v1 <- mmad_var(1); v2 <- mmad_var(2)
  hand_built <- 12 * log(0.5 * v1 + 0.5 * v2) +
                15 * log((2/3) * v1 + (1/3) * v2) +
                9  * log((1/3) * v1 + (2/3) * v2) -
                6 * v1 - 6 * v2

  expect_equivalent_expr(parsed, hand_built,
                         list(c(4, 2), c(1, 1), c(0.7, 0.3)))
})

# ---- Refusal cases -------------------------------------------------------

test_that("expr * expr is refused with a clear message", {
  expect_error(
    as_mmad_expr(~ theta[1] * theta[2], init = c(1, 1)),
    "non-DCP")
})

test_that("numeric / expr is refused", {
  expect_error(
    as_mmad_expr(~ 1 / theta[1], init = c(1)),
    "division by an expression")
})

test_that("expression exponent is refused", {
  expect_error(
    as_mmad_expr(~ theta[1] ^ theta[2], init = c(1, 1)),
    "exponent")
})

test_that("unsupported function names produce a helpful error", {
  expect_error(
    as_mmad_expr(~ tan(theta[1]), init = c(1)),
    "Unsupported expression")
})

test_that("an unresolvable symbol produces an error pointing at the name", {
  expect_error(
    as_mmad_expr(~ undefined_thing + theta[1], init = c(1)),
    "Cannot resolve symbol 'undefined_thing'")
})

test_that("non-integer theta index is rejected", {
  expect_error(
    as_mmad_expr(~ theta[1.5], init = c(1, 1)),
    "positive integer")
})

# ---- Regression: bare numerics coerce in Ops.mmad_expr -------------------
# Phase 5 redefined as_mmad_expr() as an S3 generic; for a while it had
# no .numeric method, which broke `1 - mmad_expr` style hand-built trees.
# Pin the behaviour here.

test_that("bare numeric scalars coerce on either side of Ops.mmad_expr", {
  hand_built <- log(1 - 0.5 * mmad_var(1) - 0.5 * mmad_var(2))
  parsed     <- as_mmad_expr(~ log(1 - 0.5 * theta[1] - 0.5 * theta[2]),
                             init = c(0, 0))
  for (theta in list(c(0.2, 0.2), c(0.5, 0.4), c(0.1, 0.7))) {
    a <- evaluate_expr(hand_built, theta)
    b <- evaluate_expr(parsed,     theta)
    expect_equal(a$value,    b$value,    tolerance = 1e-12)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-12)
  }
})

test_that("as_mmad_expr() dispatches correctly across all input kinds", {
  expect_s3_class(as_mmad_expr(3.5),                        "mmad_const")
  expect_s3_class(as_mmad_expr(mmad_var(1)),                "mmad_var")
  expect_s3_class(as_mmad_expr(~ theta[1] + 2, init = 0),   "mmad_expr")
  expect_error(   as_mmad_expr("not a numeric"),            "Cannot coerce")
  expect_error(   as_mmad_expr(c(1, 2)),                    "scalar numerics")
})

# ---- sum() over observations and X %*% theta -----------------------------

test_that("theta-free sum() collapses to a numeric constant (regression for Phase 5.5)", {
  # Phase 5.5 added row-wise expansion for sum() with theta-dependent
  # inners. Theta-free sums must still collapse to a scalar so they can
  # be used as coefficients via `sum(X) * theta[1]` etc.
  X <- c(1, 2, 3, 4)
  parsed <- as_mmad_expr(~ sum(X) * theta[1] - mean(X) * theta[2],
                         init = c(0, 0))
  hand_built <- 10 * mmad_var(1) - 2.5 * mmad_var(2)
  for (theta in list(c(0, 0), c(1, 1), c(2, -1))) {
    a <- evaluate_expr(parsed,     theta)
    b <- evaluate_expr(hand_built, theta)
    expect_equal(a$value,    b$value,    tolerance = 1e-12)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-12)
  }
})

test_that("sum() over a vector data reference expands row-by-row", {
  X <- c(1, 2, 3, 4)
  parsed <- as_mmad_expr(~ sum(X * theta[1]),
                         init = c(0),
                         data = list(X = X))
  hand_built <- 1 * mmad_var(1) + 2 * mmad_var(1) +
                3 * mmad_var(1) + 4 * mmad_var(1)
  for (theta in list(c(1), c(0.5), c(-2))) {
    a <- evaluate_expr(parsed,     theta)
    b <- evaluate_expr(hand_built, theta)
    expect_equal(a$value,    b$value,    tolerance = 1e-12)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-12)
  }
})

test_that("X %*% theta inside sum() expands to row-i affine combos", {
  X <- matrix(c(1, 2, 3, 4), nrow = 2)   # 2x2
  parsed <- as_mmad_expr(~ sum(X %*% theta),
                         init = c(0, 0),
                         data = list(X = X))
  # Expected: (X[1,1]*theta[1] + X[1,2]*theta[2]) + (X[2,1]*theta[1] + X[2,2]*theta[2])
  #         = (1*theta[1] + 3*theta[2]) + (2*theta[1] + 4*theta[2])
  #         = 3 * theta[1] + 7 * theta[2]
  hand_built <- 3 * mmad_var(1) + 7 * mmad_var(2)
  for (theta in list(c(1, 1), c(2, -1), c(0.5, 0.5))) {
    a <- evaluate_expr(parsed,     theta)
    b <- evaluate_expr(hand_built, theta)
    expect_equal(a$value,    b$value,    tolerance = 1e-12)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-12)
  }
})

test_that("Poisson-style sum: y * (X %*% theta) - exp(X %*% theta) parses", {
  set.seed(1)
  n <- 4; p <- 2
  X <- matrix(rnorm(n * p), nrow = n)
  y <- c(1, 2, 0, 3)
  parsed <- as_mmad_expr(
    ~ sum(y * (X %*% theta) - exp(X %*% theta)),
    init = c(0, 0),
    data = list(X = X, y = y))

  # Hand-built: row-by-row sum.
  hand_built <- mmad_const(0)
  for (i in seq_len(n)) {
    affine_i <- X[i, 1] * mmad_var(1) + X[i, 2] * mmad_var(2)
    hand_built <- hand_built + y[i] * affine_i - exp(affine_i)
  }

  for (theta in list(c(0, 0), c(0.5, -0.3), c(1, 1))) {
    a <- evaluate_expr(parsed,     theta)
    b <- evaluate_expr(hand_built, theta)
    expect_equal(a$value,    b$value,    tolerance = 1e-12)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-12)
  }
})

# ---- Refusal cases for the new grammar -----------------------------------

test_that("X %*% theta outside sum() is refused with a helpful message", {
  X <- matrix(1:4, nrow = 2)
  expect_error(
    as_mmad_expr(~ X %*% theta, init = c(0, 0), data = list(X = X)),
    "only supported inside sum")
})

test_that("bare 'theta' symbol outside %*% is refused", {
  expect_error(
    as_mmad_expr(~ theta + 1, init = c(0)),
    "Bare 'theta'")
})

test_that("X %*% theta with mismatched dimensions errors", {
  X <- matrix(1:6, nrow = 2)   # 2 rows, 3 cols
  expect_error(
    as_mmad_expr(~ sum(X %*% theta), init = c(0, 0), data = list(X = X)),
    "must be a numeric matrix")
})

test_that("inconsistent vector lengths inside sum() error", {
  X <- c(1, 2, 3); Y <- c(1, 2, 3, 4)
  expect_error(
    as_mmad_expr(~ sum(X * theta[1] + Y * theta[2]),
                 init = c(0, 0),
                 data = list(X = X, Y = Y)),
    "inconsistent vector lengths")
})
