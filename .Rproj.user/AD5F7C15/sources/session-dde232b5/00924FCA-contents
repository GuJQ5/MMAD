#' MMAD: Minorization-Maximization via Assembly-Decomposition Technology
#'
#' Tools for maximizing a target function via the MM algorithm, built
#' around three layers:
#'
#' * a symbolic expression tree (`mmad_expr`) constructed via the
#'   formula interface or operator overloads on [mmad_var()];
#' * disciplined-convex-programming inference of curvature, sign, and
#'   monotonicity ([curvature()], [sign_of()], [is_dcp()]);
#' * a surrogate construction ([minorize_at()]) that uses only Jensen's
#'   inequality and the supporting hyperplane, plus a top-level driver
#'   ([mmad()]) that maximizes the surrogate via block-coordinate Newton.
#'
#' The formula interface accepts `theta[i]` indexing or named parameters,
#' eager-evaluates theta-free sub-expressions in `data`, and recognizes
#' `sum()` over observations together with `X %*% theta`, so standard
#' statistical likelihoods can be written in one line. See [mmad()] for
#' the headline example and [Function_check()] for a pre-flight diagnostic.
#'
#' @keywords internal
"_PACKAGE"
