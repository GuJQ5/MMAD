## Function_check() tests.
##
## Coverage:
##   - Structural fields are populated correctly (is_minorizable, expr_tree,
##     summands$minorizable column added in the tree-output revision).
##   - Domain failures are detected and don't raise errors.
##   - The print method runs without error and emits the tree header.
##   - Expression tree: correct node labels, curvature, sign, minorizability.
##   - Extended rules reflected in the tree:
##       E1 additive Jensen  -- log(1 + exp(theta))
##       E2 hyperplane-concave-inner -- -log(1 + exp(t1) + exp(t2) + exp(t3))
##       E3 log(exp(h)) = h
##       E4 exp(log(h)) = h
##       E5 log(h^c) = c * log(h)

# ── Existing tests (updated for new fields) ─────────────────────────────────

test_that("multinomial objective is reported as concave, DCP, and fully separable", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- 12 * log(0.5 * v1 + 0.5 * v2) +
          15 * log((2/3) * v1 + (1/3) * v2) +
           9 * log((1/3) * v1 + (2/3) * v2) -
           6 * v1 - 6 * v2

  chk <- Function_check(expr, init = c(4, 2))

  expect_s3_class(chk, "mmad_check")
  expect_equal(chk$target_curvature, "concave")
  expect_true(chk$is_dcp)
  expect_true(chk$is_minorizable)        # new field: mirrors is_separable
  expect_true(chk$domain_ok)
  expect_true(chk$is_separable)
  expect_equal(nrow(chk$summands), 5L)   # three logs + two linear terms
  # All summands should be minorizable, and curvature/sign columns are absent.
  expect_true(all(chk$summands$minorizable == "yes"))
  expect_false("curvature" %in% names(chk$summands))
  expect_false("sign"      %in% names(chk$summands))
  # expr_tree is present and is a list with the expected fields.
  expect_type(chk$expr_tree, "list")
  expect_named(chk$expr_tree,
               c("label", "curvature", "sign", "minorizable", "children"),
               ignore.order = TRUE)
})

test_that("formula entry point lowers and runs the check", {
  chk <- Function_check(~ log(theta[1] + theta[2]) - theta[1],
                        init = c(1, 1))
  expect_true(chk$is_dcp)
  expect_equal(chk$target_curvature, "concave")
  expect_true(chk$is_minorizable)
  expect_type(chk$expr_tree, "list")
})

test_that("non-separable residue is flagged with the offending coordinates", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- log(1 - 0.5 * v1 - 0.5 * v2)   # mixed-sign Jensen slots
  chk <- Function_check(expr, init = c(0.3, 0.3))
  expect_false(chk$is_separable)
  expect_false(chk$is_minorizable)
  expect_setequal(chk$non_separable_indices, c(1L, 2L))
  # The top-level summand for the log node should be "no".
  expect_true(any(chk$summands$minorizable == "no"))
})

test_that("domain failure at init is reported, not raised", {
  v1 <- mmad_var(1)
  expr <- log(v1)
  chk <- Function_check(expr, init = c(-1))   # log(-1) is undefined
  expect_false(chk$domain_ok)
  expect_match(chk$domain_message, "log")
  # is_separable is NA when domain check fails.
  expect_true(is.na(chk$is_separable))
  expect_false(chk$is_minorizable)
  # Tree is still built (with minorizable = "unknown" for non-trivial nodes).
  expect_type(chk$expr_tree, "list")
})

test_that("non-DCP target is flagged but the check still completes", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- exp(v1) - exp(v2)               # convex - convex = unknown
  chk <- Function_check(expr, init = c(0, 0))
  expect_false(chk$is_dcp)
  expect_equal(chk$target_curvature, "unknown")
  # Univariate extraction still gives a separable surrogate here.
  expect_true(chk$is_separable)
  expect_true(chk$is_minorizable)
})

# ── Print method ─────────────────────────────────────────────────────────────

test_that("print.mmad_check does not show tree by default", {
  v1  <- mmad_var(1)
  chk <- Function_check(log(v1), init = c(1))
  out <- capture.output(print(chk))
  expect_false(any(grepl("Expression tree", out)))
})

test_that("print.mmad_check shows tree when tree = TRUE at call time", {
  v1  <- mmad_var(1)
  chk <- Function_check(log(v1), init = c(1))
  out <- capture.output(print(chk, tree = TRUE))
  expect_true(any(grepl("Expression tree", out)))
  expect_true(any(grepl("minorizable", out, ignore.case = TRUE)))
})

test_that("print.mmad_check shows tree when tree = TRUE at construction time", {
  v1  <- mmad_var(1)
  chk <- Function_check(log(v1), init = c(1), tree = TRUE)
  out <- capture.output(print(chk))
  expect_true(any(grepl("Expression tree", out)))
})

# ── Expression-tree structure ─────────────────────────────────────────────────

test_that("build_expr_tree annotates a leaf theta node correctly", {
  v1 <- mmad_var(1)
  tree <- build_expr_tree(v1, init_vec = c(2))
  expect_equal(tree$label,       "theta[1]")
  expect_equal(tree$curvature,   "affine")
  expect_equal(tree$minorizable, "trivial")
  expect_length(tree$children,  0L)
})

test_that("build_expr_tree annotates a constant node correctly", {
  cn <- mmad_const(3.5)
  tree <- build_expr_tree(cn, init_vec = c(1))
  expect_equal(tree$curvature,   "affine")
  expect_equal(tree$sign,        "positive")
  expect_equal(tree$minorizable, "trivial")
  expect_length(tree$children,  0L)
})

test_that("build_expr_tree: log(theta[1]) is concave, univariate -> trivial", {
  v1   <- mmad_var(1, sign = "positive")
  tree <- build_expr_tree(log(v1), init_vec = c(1))
  expect_equal(tree$curvature,   "concave")
  expect_equal(tree$minorizable, "trivial")   # univariate
  expect_length(tree$children,  1L)
  expect_equal(tree$children[[1L]]$label, "theta[1]")
})

test_that("build_expr_tree: affine sum has two children and 'trivial' flag", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  tree <- build_expr_tree(v1 + v2, init_vec = c(1, 1))
  expect_equal(tree$curvature,   "affine")
  expect_equal(tree$minorizable, "trivial")
  # Smart constructor flattens add so we get a 2-child add node.
  expect_equal(tree$label, "(+) add")
  expect_length(tree$children, 2L)
})

test_that("build_expr_tree: scale node label includes coefficient", {
  v1   <- mmad_var(1)
  tree <- build_expr_tree(3 * v1, init_vec = c(1))
  expect_match(tree$label, "scale")
  expect_match(tree$label, "3")
})

test_that("build_expr_tree: pow node label includes exponent", {
  v1   <- mmad_var(1, sign = "positive")
  tree <- build_expr_tree(v1 ^ 2, init_vec = c(1))
  expect_match(tree$label, "pow")
  expect_match(tree$label, "2")
})

test_that("build_expr_tree: log(t1 + t2) is concave and minorizable", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  tree <- build_expr_tree(log(v1 + v2), init_vec = c(1, 1))
  expect_equal(tree$curvature,   "concave")
  expect_equal(tree$minorizable, "yes")
  # Child of log is the add node.
  add_node <- tree$children[[1L]]
  expect_equal(add_node$label,     "(+) add")
  expect_equal(add_node$curvature, "affine")
})

test_that("build_expr_tree: non-separable log(1 - t1 - t2) is flagged 'no'", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  tree <- build_expr_tree(log(1 - 0.5 * v1 - 0.5 * v2),
                          init_vec = c(0.3, 0.3))
  expect_equal(tree$curvature,   "concave")
  expect_equal(tree$minorizable, "no")
})

# ── Extended rules visible in the tree ───────────────────────────────────────

test_that("E1 additive Jensen: log(1 + exp(t1)) is concave and minorizable", {
  # Standard DCP cannot prove curvature (convex inner), but E1 applies.
  v1   <- mmad_var(1)
  expr <- log(mmad_const(1) + exp(v1))
  chk  <- Function_check(expr, init = c(1))

  expect_equal(chk$target_curvature, "concave")
  expect_true(chk$is_minorizable)

  # Root of tree is the log node.
  tree <- chk$expr_tree
  expect_equal(tree$label,       "log")
  expect_equal(tree$curvature,   "concave")
  expect_equal(tree$minorizable, "trivial")   # univariate => trivial

  # Child is the add node with two children: const(1) and exp(t1).
  add_node <- tree$children[[1L]]
  expect_equal(add_node$label, "(+) add")
  labels <- vapply(add_node$children, `[[`, character(1), "label")
  expect_true(any(labels == "1"))
  expect_true(any(labels == "exp"))
})

test_that("E1 additive Jensen: log(1+exp(t1)+exp(t2)) is concave and minorizable (multivariate)", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- log(mmad_const(1) + exp(v1) + exp(v2))
  chk  <- Function_check(expr, init = c(1, 1))

  expect_equal(chk$target_curvature, "concave")
  expect_true(chk$is_minorizable)

  tree <- chk$expr_tree
  expect_equal(tree$label,       "log")
  expect_equal(tree$curvature,   "concave")
  expect_equal(tree$minorizable, "yes")
})

test_that("E2 hyperplane-concave-inner: -log(1+exp(t1)+exp(t2)+exp(t3)) is minorizable", {
  v1 <- mmad_var(1); v2 <- mmad_var(2); v3 <- mmad_var(3)
  expr <- -60 * log(mmad_const(1) + exp(v1) + exp(v2) + exp(v3))
  chk  <- Function_check(expr, init = c(0, 0, 0))

  # The whole expression is convex (neg-log of convex inner).
  expect_equal(chk$target_curvature, "convex")
  expect_true(chk$is_minorizable)

  # The root is the scale node; its child should be the log (convex, minorizable).
  tree  <- chk$expr_tree
  # After scale/neg unwrapping, expect a scale or neg node at the top.
  expect_true(grepl("scale|neg", tree$label))
  expect_equal(tree$curvature,   "convex")
  expect_equal(tree$minorizable, "yes")
})

test_that("E3 log(exp(h)): curvature equals that of h", {
  v1   <- mmad_var(1)
  expr <- log(exp(v1))    # should reduce to affine (= v1)
  chk  <- Function_check(expr, init = c(1))

  expect_equal(chk$target_curvature, "affine")
  expect_true(chk$is_minorizable)

  tree <- chk$expr_tree
  expect_equal(tree$label,     "log")
  expect_equal(tree$curvature, "affine")
})

test_that("E3 log(exp(h)) with concave h: curvature is concave and minorizable", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  # h = log(t1) + log(t2) is concave; log(exp(h)) = h is concave.
  # simplify_expr() rewrites log(exp(h)) -> h before minorize_at() runs,
  # so the minorization engine sees log(t1) + log(t2) directly and
  # decomposes it into two univariate terms.
  inner <- log(mmad_var(1, sign = "positive")) +
           log(mmad_var(2, sign = "positive"))
  expr  <- log(exp(inner))
  chk   <- Function_check(expr, init = c(1, 1))

  expect_equal(chk$target_curvature, "concave")
  expect_true(chk$is_minorizable)
})

test_that("E4 exp(log(h)): curvature equals that of h when h > 0", {
  v1   <- mmad_var(1, sign = "positive")
  expr <- exp(log(v1))    # should reduce to affine (= v1, since v1 > 0)
  chk  <- Function_check(expr, init = c(1))

  expect_equal(chk$target_curvature, "affine")
  expect_true(chk$is_minorizable)

  tree <- chk$expr_tree
  expect_equal(tree$label,     "exp")
  expect_equal(tree$curvature, "affine")
})

test_that("E5 log(h^c) with c > 0 and h > 0: curvature is concave", {
  v1   <- mmad_var(1, sign = "positive")
  # log(v1^2) = 2 * log(v1), which is concave.
  expr <- log(v1 ^ 2)
  chk  <- Function_check(expr, init = c(1))

  expect_equal(chk$target_curvature, "concave")
  expect_true(chk$is_minorizable)   # univariate

  tree <- chk$expr_tree
  expect_equal(tree$label,     "log")
  expect_equal(tree$curvature, "concave")
})

test_that("E5 log(h^c) with c < 0 and h > 0: curvature is convex", {
  v1   <- mmad_var(1, sign = "positive")
  # log(v1^(-1)) = -log(v1), which is convex.
  expr <- log(v1 ^ (-1))
  chk  <- Function_check(expr, init = c(1))

  expect_equal(chk$target_curvature, "convex")
})

# ── summands$minorizable column ──────────────────────────────────────────────

test_that("summands data.frame carries a minorizable column with expected values", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- log(0.5 * v1 + 0.5 * v2) - v1
  chk  <- Function_check(expr, init = c(2, 2))

  expect_true("minorizable" %in% names(chk$summands))
  expect_false("curvature" %in% names(chk$summands))
  expect_false("sign"      %in% names(chk$summands))
  expect_true(all(chk$summands$minorizable %in%
                  c("yes", "no", "trivial", "error", "unknown")))
  # log summand: multivariate concave => "yes".
  # -v1 summand: univariate so minorize_at returns is_separable=TRUE => "yes".
  # summands$minorizable uses the dry-run path ("yes"/"no"), never "trivial".
  minor_vals <- chk$summands$minorizable
  expect_true(all(minor_vals == "yes"))
})

test_that("is_minorizable is FALSE when the whole expression is not separable", {
  v1 <- mmad_var(1); v2 <- mmad_var(2)
  expr <- log(1 - 0.5 * v1 - 0.5 * v2)
  chk  <- Function_check(expr, init = c(0.3, 0.3))
  expect_false(chk$is_minorizable)
  expect_false(chk$is_separable)
})
