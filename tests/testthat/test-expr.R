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
