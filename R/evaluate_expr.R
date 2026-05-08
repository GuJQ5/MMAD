## Forward-mode evaluator for `mmad_expr`.
##
## Recurses through the expression tree. At each node it returns
##   list(value = scalar, gradient = numeric(p), hessian = matrix(p, p))
## where derivatives are taken with respect to theta (not the atom's local
## arguments). Atom-local derivatives come from the registry; the chain
## rule is applied at every internal node.
##
## Chain rule used at an atom application g(a_1, ..., a_k) where each
## a_j(theta) has gradient g_j and Hessian H_j w.r.t. theta and the atom
## reports local gradient lg[j] and local Hessian lH[j, m]:
##
##   d g / d theta_i           = sum_j lg[j] * g_j[i]
##   d^2 g / d theta_i d theta_l
##                             = sum_{j,m} lH[j, m] * g_j[i] * g_m[l]
##                             + sum_j lg[j] * H_j[i, l]
##
## This is plain forward-mode AD to second order. It is O(k * p^2) per
## node; for the small parameter dimensions used in MM applications that
## is plenty fast.

#' Evaluate an `mmad_expr` at a parameter vector
#'
#' Computes the value, gradient, and Hessian of the expression at `theta`,
#' all with respect to theta. The gradient is a numeric vector of length
#' `length(theta)`; the Hessian is a `length(theta)` by `length(theta)`
#' numeric matrix.
#'
#' This is the Phase 1 replacement for `Function_evaluation()` -- it works
#' on the new expression tree. The legacy function continues to work on
#' the old list representation; `legacy_to_expr()` bridges between them.
#'
#' @param expr  An object of class `mmad_expr`.
#' @param theta A numeric vector of parameter values.
#'
#' @return A list with components `value`, `gradient`, and `hessian`.
#' @examples
#' expr <- log(mmad_var(1) + mmad_var(2))
#' evaluate_expr(expr, c(1, 1))
#' @export
evaluate_expr <- function(expr, theta) {
  if (!inherits(expr, "mmad_expr")) {
    stop("evaluate_expr(): expr must be an mmad_expr object.")
  }
  if (!is.numeric(theta) || length(theta) == 0L) {
    stop("evaluate_expr(): theta must be a non-empty numeric vector.")
  }
  theta <- as.numeric(theta)
  p <- length(theta)
  eval_mmad_node(expr, theta, p)
}

# Internal recursion. Kept separate from the public entry point so that the
# arity / type checks above run only once per top-level call.
eval_mmad_node <- function(node, theta, p) {
  if (inherits(node, "mmad_var")) {
    i <- node$index
    if (i > p) {
      stop(sprintf("theta[%d] referenced but theta has length %d.", i, p))
    }
    g <- numeric(p)
    g[i] <- 1
    return(list(value    = theta[i],
                gradient = g,
                hessian  = matrix(0, p, p)))
  }

  if (inherits(node, "mmad_const")) {
    return(list(value    = node$value,
                gradient = numeric(p),
                hessian  = matrix(0, p, p)))
  }

  if (inherits(node, "mmad_call")) {
    children <- lapply(node$args, eval_mmad_node, theta = theta, p = p)
    arg_vals <- vapply(children, function(ch) ch$value, numeric(1))
    atom     <- mmad_atom(node$op)

    V  <- atom$value(arg_vals, node$params)
    lg <- atom$grad(arg_vals, node$params)
    lH <- atom$hess(arg_vals, node$params)
    k  <- length(arg_vals)

    if (length(lg) != k) {
      stop(sprintf("Atom '%s' grad() returned length %d, expected %d.",
                   node$op, length(lg), k))
    }
    if (!is.matrix(lH) || nrow(lH) != k || ncol(lH) != k) {
      stop(sprintf("Atom '%s' hess() returned a non-%dx%d matrix.",
                   node$op, k, k))
    }

    # First-order chain rule.
    g_tot <- numeric(p)
    for (j in seq_len(k)) {
      if (lg[j] != 0) {
        g_tot <- g_tot + lg[j] * children[[j]]$gradient
      }
    }

    # Second-order chain rule.
    H_tot <- matrix(0, p, p)
    for (j in seq_len(k)) {
      if (lg[j] != 0) {
        H_tot <- H_tot + lg[j] * children[[j]]$hessian
      }
    }
    for (j in seq_len(k)) {
      for (m in seq_len(k)) {
        h_jm <- lH[j, m]
        if (h_jm != 0) {
          H_tot <- H_tot +
            h_jm * tcrossprod(children[[j]]$gradient, children[[m]]$gradient)
        }
      }
    }

    return(list(value = V, gradient = g_tot, hessian = H_tot))
  }

  stop("evaluate_expr(): unrecognised expression node of class ",
       paste(class(node), collapse = "/"))
}
