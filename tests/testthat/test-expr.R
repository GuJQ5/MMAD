## Phase 1 tests: the expression tree, the atom registry, and the evaluator.
## We pin correctness in two independent ways:
##   1. closed-form check on a simple operator-built expression;
##   2. central-difference numerical derivatives on a mixed expression.

# ---- Operator-built expression: closed-form derivatives ------------------

test_that("operator-built expressions evaluate to closed-form derivatives", {
  th1 <- mmad_var(1)
  th2 <- mmad_var(2)
  # f(theta) = log(theta[1] + 2 * theta[2]) - 0.5 * exp(theta[1])
  expr <- log(th1 + 2 * th2) - 0.5 * exp(th1)

  expect_s3_class(expr, "mmad_expr")

  theta <- c(0.5, 0.7)
  res <- evaluate_expr(expr, theta)

  s   <- theta[1] + 2 * theta[2]
  e1  <- exp(theta[1])
  v_expected <- log(s) - 0.5 * e1
  g_expected <- c(1 / s - 0.5 * e1,
                  2 / s)
  H_expected <- matrix(c(
    -1 / s^2 - 0.5 * e1,  -2 / s^2,
    -2 / s^2,             -4 / s^2
  ), nrow = 2, byrow = TRUE)

  expect_equal(res$value,    v_expected, tolerance = 1e-12)
  expect_equal(res$gradient, g_expected, tolerance = 1e-12)
  expect_equal(res$hessian,  H_expected, tolerance = 1e-12)
})

# ---- Numerical-derivative check on a mixed expression ---------------------

test_that("evaluate_expr derivatives agree with central differences", {
  th1 <- mmad_var(1); th2 <- mmad_var(2)
  expr <- log(0.3 * th1 + 0.7 * th2) + (th1 + th2)^2 + 0.25 * exp(th1)

  theta <- c(1.2, 0.8)
  res <- evaluate_expr(expr, theta)

  f <- function(t) evaluate_expr(expr, t)$value
  eps <- 1e-5
  num_grad <- numeric(2)
  for (i in 1:2) {
    tp <- theta; tp[i] <- tp[i] + eps
    tm <- theta; tm[i] <- tm[i] - eps
    num_grad[i] <- (f(tp) - f(tm)) / (2 * eps)
  }
  expect_equal(res$gradient, num_grad, tolerance = 1e-6)

  num_hess <- matrix(0, 2, 2)
  for (i in 1:2) for (j in 1:2) {
    tpp <- theta; tpp[i] <- tpp[i] + eps; tpp[j] <- tpp[j] + eps
    tpm <- theta; tpm[i] <- tpm[i] + eps; tpm[j] <- tpm[j] - eps
    tmp <- theta; tmp[i] <- tmp[i] - eps; tmp[j] <- tmp[j] + eps
    tmm <- theta; tmm[i] <- tmm[i] - eps; tmm[j] <- tmm[j] - eps
    num_hess[i, j] <- (f(tpp) - f(tpm) - f(tmp) + f(tmm)) / (4 * eps^2)
  }
  expect_equal(res$hessian, num_hess, tolerance = 1e-4)
})

# ---- Construction-time validation ----------------------------------------

test_that("constructor input validation rejects bad arguments", {
  expect_error(mmad_var(0),     "positive integer")
  expect_error(mmad_var(1.5),   "positive integer")
  expect_error(mmad_const(NA),  "finite numeric scalar")
  expect_error(mmad_const(c(1, 2)), "finite numeric scalar")
  expect_error(mmad_call("not_a_real_atom", list(mmad_var(1))),
               "Unknown atom")
})

test_that("Ops generic forbids non-DCP combinations at build time", {
  th1 <- mmad_var(1); th2 <- mmad_var(2)
  expect_error(th1 * th2, "non-DCP")
  expect_error(1 / th1,   "division by a numeric")
  expect_error(th1 ^ th2, "numeric scalar")
})

# ---- simplify_expr: algebraic identity rewrites ---------------------------

test_that("simplify_expr leaves leaves unchanged", {
  v1 <- mmad_var(1)
  cn <- mmad_const(3)
  expect_identical(simplify_expr(v1), v1)
  expect_identical(simplify_expr(cn), cn)
})

test_that("E3: simplify_expr rewrites log(exp(h)) to h", {
  v1   <- mmad_var(1)
  expr <- log(exp(v1))
  s    <- simplify_expr(expr)
  # Result should be exactly the theta[1] node.
  expect_true(inherits(s, "mmad_var"))
  expect_equal(s$index, 1L)
  # Value and gradient are preserved.
  res_orig <- evaluate_expr(expr, c(2))
  res_simp <- evaluate_expr(s,    c(2))
  expect_equal(res_simp$value,    res_orig$value,    tolerance = 1e-12)
  expect_equal(res_simp$gradient, res_orig$gradient, tolerance = 1e-12)
})

test_that("E3: simplify_expr rewrites log(exp(h)) for a non-trivial h", {
  # h = log(t1) + log(t2); log(exp(h)) should simplify to h.
  v1 <- mmad_var(1, sign = "positive")
  v2 <- mmad_var(2, sign = "positive")
  h    <- log(v1) + log(v2)
  expr <- log(exp(h))
  s    <- simplify_expr(expr)
  # s should be an add node (the original h), not a log node.
  expect_true(inherits(s, "mmad_call"))
  expect_equal(s$op, "add")
  # Numerical value check.
  theta <- c(1.5, 2.0)
  expect_equal(evaluate_expr(s, theta)$value,
               evaluate_expr(expr, theta)$value, tolerance = 1e-12)
})

test_that("E3: nested log(exp(log(exp(h)))) collapses to h in one pass", {
  v1   <- mmad_var(1)
  expr <- log(exp(log(exp(v1))))
  s    <- simplify_expr(expr)
  expect_true(inherits(s, "mmad_var"))
  expect_equal(s$index, 1L)
})

test_that("E4: simplify_expr rewrites exp(log(h)) to h when h > 0", {
  v1   <- mmad_var(1, sign = "positive")
  expr <- exp(log(v1))
  s    <- simplify_expr(expr)
  expect_true(inherits(s, "mmad_var"))
  expect_equal(s$index, 1L)
  theta <- c(3.0)
  expect_equal(evaluate_expr(s, theta)$value,
               evaluate_expr(expr, theta)$value, tolerance = 1e-12)
})

test_that("E4: simplify_expr does NOT rewrite exp(log(h)) when h sign unknown", {
  v1   <- mmad_var(1)   # sign = "unknown"
  expr <- exp(log(v1))
  s    <- simplify_expr(expr)
  # Should be left as exp(log(theta[1])), i.e. an mmad_call with op "exp".
  expect_true(inherits(s, "mmad_call"))
  expect_equal(s$op, "exp")
})

test_that("E5: simplify_expr rewrites log(h^c) to c*log(h) when h > 0", {
  v1   <- mmad_var(1, sign = "positive")
  # log(v1^3) should become scale(3, log(v1)).
  expr <- log(v1 ^ 3)
  s    <- simplify_expr(expr)
  expect_true(inherits(s, "mmad_call"))
  expect_equal(s$op, "scale")
  expect_equal(s$params$c, 3)
  # The inner of scale should be log(v1).
  inner <- s$args[[1L]]
  expect_equal(inner$op, "log")
  # Numerical value check at a positive point.
  theta <- c(2.0)
  expect_equal(evaluate_expr(s, theta)$value,
               evaluate_expr(expr, theta)$value, tolerance = 1e-12)
  expect_equal(evaluate_expr(s, theta)$gradient,
               evaluate_expr(expr, theta)$gradient, tolerance = 1e-10)
})

test_that("E5: simplify_expr rewrites log(h^(-1)) to -log(h) when h > 0", {
  v1   <- mmad_var(1, sign = "positive")
  expr <- log(v1 ^ (-1))
  s    <- simplify_expr(expr)
  # scale(-1, log(v1)) -- the smart constructor might leave it as scale or neg.
  theta <- c(2.5)
  expect_equal(evaluate_expr(s, theta)$value,
               evaluate_expr(expr, theta)$value, tolerance = 1e-12)
})

test_that("E5: simplify_expr does NOT rewrite log(h^c) when h sign unknown", {
  v1   <- mmad_var(1)   # sign = "unknown"
  expr <- log(v1 ^ 2)
  s    <- simplify_expr(expr)
  # Should remain log(pow(...)).
  expect_true(inherits(s, "mmad_call"))
  expect_equal(s$op, "log")
  expect_equal(s$args[[1L]]$op, "pow")
})

test_that("simplify_expr rewrites children before the parent (bottom-up)", {
  # log(exp(exp(log(v1_pos)))) -- inner exp(log(v1)) -> v1 first,
  # then log(exp(v1)) -> v1 at the outer level.
  v1   <- mmad_var(1, sign = "positive")
  expr <- log(exp(exp(log(v1))))
  s    <- simplify_expr(expr)
  expect_true(inherits(s, "mmad_var"))
  expect_equal(s$index, 1L)
})

test_that("simplify_expr preserves values and gradients for an unchanged tree", {
  # An expression with no applicable rules should be returned unmodified in value.
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- log(0.4 * v1 + 0.6 * v2) + 2 * exp(v1)
  s    <- simplify_expr(expr)
  theta <- c(1.0, 2.0)
  res_orig <- evaluate_expr(expr, theta)
  res_simp <- evaluate_expr(s,    theta)
  expect_equal(res_simp$value,    res_orig$value,    tolerance = 1e-12)
  expect_equal(res_simp$gradient, res_orig$gradient, tolerance = 1e-12)
})
