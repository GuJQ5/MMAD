## Phase 2 corpus: each entry pairs an expression with its expected
## DCP-inferred curvature. The tests cover the rules that matter most
## for MMAD's intended use cases: affine combinations, log of affine
## (multinomial-style), exp of affine, even-integer powers, signed-input
## sqrt/recip, and a handful of canonical DCP false-negatives that we
## explicitly want reported as "unknown" rather than misclassified.

# Small helper to keep each row of the corpus visually compact.
expect_curvature <- function(expr, expected) {
  expect_equal(curvature(expr), expected,
               info = paste0("expression: ", format_mmad_expr(expr)))
}
expect_sign <- function(expr, expected) {
  expect_equal(sign_of(expr), expected,
               info = paste0("expression: ", format_mmad_expr(expr)))
}

# ---- Affine corpus --------------------------------------------------------

test_that("leaves and affine combinations are recognised as affine", {
  expect_curvature(mmad_var(1),                                  "affine")
  expect_curvature(mmad_const(5),                                "affine")
  expect_curvature(mmad_const(-3),                               "affine")
  expect_curvature(mmad_const(0),                                "affine")
  expect_curvature(mmad_var(1) + mmad_var(2),                    "affine")
  expect_curvature(mmad_var(1) - mmad_var(2),                    "affine")
  expect_curvature(2 * mmad_var(1) - 3 * mmad_var(2) + 1,        "affine")
  expect_curvature(-mmad_var(1),                                 "affine")
  expect_curvature(mmad_var(1) ^ 1,                              "affine")
  expect_curvature(mmad_var(1) ^ 0,                              "affine")
})

# ---- Convex corpus --------------------------------------------------------

test_that("exp/even-power compositions are recognised as convex", {
  expect_curvature(exp(mmad_var(1)),                             "convex")
  expect_curvature(exp(mmad_var(1) + 2 * mmad_var(2)),           "convex")
  expect_curvature(mmad_var(1) ^ 2,                              "convex")
  expect_curvature((mmad_var(1) + mmad_var(2)) ^ 2,              "convex")
  expect_curvature(mmad_var(1) ^ 4,                              "convex")
  expect_curvature(exp(mmad_var(1)) + exp(mmad_var(2)),          "convex")
  expect_curvature(exp(mmad_var(1)) + mmad_var(1) ^ 2,           "convex")
  expect_curvature(2 * exp(mmad_var(1)),                         "convex")
})

# ---- Concave corpus -------------------------------------------------------

test_that("log/neg-exp compositions are recognised as concave", {
  expect_curvature(log(mmad_var(1) + mmad_var(2)),               "concave")
  expect_curvature(log(mmad_var(1)) + log(mmad_var(2)),          "concave")
  expect_curvature(-exp(mmad_var(1)),                            "concave")
  expect_curvature(log(mmad_var(1)) - exp(mmad_var(2)),          "concave")
  expect_curvature(2 * log(mmad_var(1) + mmad_var(2)),           "concave")
})

test_that("the canonical multinomial-style objective is concave", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- 12 * log(0.5 * v1 + 0.5 * v2) +
          15 * log((2/3) * v1 + (1/3) * v2) +
          9  * log((1/3) * v1 + (2/3) * v2) -
          6 * v1 - 6 * v2
  expect_curvature(expr, "concave")
})

# ---- Sign-driven cases ----------------------------------------------------

test_that("pow with declared-sign argument resolves curvature", {
  v_pos <- mmad_var(1, sign = "positive")
  v_nn  <- mmad_var(1, sign = "nonneg")

  expect_curvature(sqrt(v_nn),                                   "concave")
  expect_curvature(v_nn ^ 0.5,                                   "concave")
  expect_curvature(v_pos ^ (-1),                                 "convex")
  expect_curvature(v_nn ^ 3,                                     "convex")
  expect_curvature(v_pos ^ 1.5,                                  "convex")
  expect_curvature(log(v_pos),                                   "concave")
})

test_that("pow without sign info on the argument falls back to unknown", {
  expect_curvature(sqrt(mmad_var(1)),                            "unknown")
  expect_curvature(mmad_var(1) ^ 0.5,                            "unknown")
  expect_curvature(mmad_var(1) ^ 3,                              "unknown")
  expect_curvature(mmad_var(1) ^ (-1),                           "unknown")
  expect_curvature(mmad_var(1) ^ 1.5,                            "unknown")
})

# ---- Genuinely unknown / DCP false-negatives ------------------------------

test_that("non-DCP combinations report unknown", {
  # convex - convex: indefinite in general
  expect_curvature(exp(mmad_var(1)) - exp(mmad_var(2)),          "unknown")
  # concave + convex: indefinite
  expect_curvature(log(mmad_var(1)) + exp(mmad_var(2)),          "unknown")
})

test_that("known DCP false-negatives report unknown rather than guessing", {
  # log(exp(x)) is x in math, but the DCP rule cannot see through.
  expect_curvature(log(exp(mmad_var(1))),                        "affine")
  # exp(log(x)) is x for x > 0; same story.
  expect_curvature(exp(log(mmad_var(1, sign = "positive"))),     "affine")
})

# ---- Sign inference -------------------------------------------------------

test_that("sign_of() propagates correctly through affine and atom layers", {
  expect_sign(mmad_const(5),                                     "positive")
  expect_sign(mmad_const(-3),                                    "negative")
  expect_sign(mmad_const(0),                                     "zero")
  expect_sign(exp(mmad_var(1)),                                  "positive")
  expect_sign(-exp(mmad_var(1)),                                 "negative")
  expect_sign(mmad_var(1) ^ 2,                                   "nonneg")
  expect_sign(log(mmad_var(1)),                                  "unknown")
  # positive + positive = positive
  expect_sign(exp(mmad_var(1)) + mmad_const(1),                  "positive")
  # negative scaling flips sign
  expect_sign(-2 * exp(mmad_var(1)),                             "negative")
})

# ---- is_dcp() convenience predicate ---------------------------------------

test_that("is_dcp() agrees with curvature()", {
  expect_true(is_dcp(mmad_var(1) + mmad_var(2)))
  expect_true(is_dcp(exp(mmad_var(1))))
  expect_true(is_dcp(log(mmad_var(1))))
  expect_false(is_dcp(exp(mmad_var(1)) - exp(mmad_var(2))))
  expect_false(is_dcp(sqrt(mmad_var(1))))
})
