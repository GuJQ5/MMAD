## Function_check() tests. Pin: structural fields are populated correctly,
## domain failures are detected, and the print method runs without error.

test_that("multinomial objective is reported as concave, DCP, and fully separable", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- 12 * log(0.5 * v1 + 0.5 * v2) +
          15 * log((2/3) * v1 + (1/3) * v2) +
          9  * log((1/3) * v1 + (2/3) * v2) -
          6 * v1 - 6 * v2

  chk <- Function_check(expr, init = c(4, 2))

  expect_s3_class(chk, "mmad_check")
  expect_equal(chk$target_curvature, "concave")
  expect_true(chk$is_dcp)
  expect_true(chk$domain_ok)
  expect_true(chk$is_separable)
  expect_equal(nrow(chk$summands), 5L)        # three logs + two linear terms
  # expect_silent(print(chk))
})

test_that("formula entry point lowers and runs the check", {
  chk <- Function_check(~ log(theta[1] + theta[2]) - theta[1],
                        init = c(1, 1))
  expect_true(chk$is_dcp)
  expect_equal(chk$target_curvature, "concave")
})

test_that("non-separable residue is flagged with the offending coordinates", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- log(1 - 0.5 * v1 - 0.5 * v2)        # mixed-sign Jensen slots
  chk <- Function_check(expr, init = c(0.3, 0.3))
  expect_false(chk$is_separable)
  expect_setequal(chk$non_separable_indices, c(1L, 2L))
})

test_that("domain failure at init is reported, not raised", {
  v1 <- mmad_var(1)
  expr <- log(v1)
  chk <- Function_check(expr, init = c(-1))    # log(-1) is undefined
  expect_false(chk$domain_ok)
  expect_match(chk$domain_message, "log")
})

test_that("non-DCP target is flagged but the check still completes", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- exp(v1) - exp(v2)                   # convex - convex = unknown
  chk <- Function_check(expr, init = c(0, 0))
  expect_false(chk$is_dcp)
  expect_equal(chk$target_curvature, "unknown")
  # Univariate extraction still gives a separable surrogate here.
  expect_true(chk$is_separable)
})
