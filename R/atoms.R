## Atom registry for the mmad_expr framework.
##
## An atom is a primitive operator/function in the expression language. Each
## atom carries metadata for two purposes:
##   1. Evaluation: value(), grad(), hess() in the atom's own argument space
##      (not theta space). The Phase 1 evaluator combines these via the
##      chain rule to produce derivatives w.r.t. theta.
##   2. DCP curvature inference (Phase 2): dcp_info(args_sign, params)
##      returns the atom's intrinsic curvature, per-argument monotonicity,
##      and output sign as a function of the children's signs and the
##      atom's parameters. The composition rule lives in R/dcp.R.

.mmad_atoms <- new.env(parent = emptyenv())

#' Register a primitive atom in the `mmad_expr` framework
#'
#' @param name         Atom name (a single character string).
#' @param arity        Integer arity, or `NA_integer_` for variadic atoms.
#' @param value,grad,hess Functions computing value, gradient, and Hessian
#'   with respect to the atom's *own arguments* (not theta). Each takes
#'   `(args, params)`, where `args` is a numeric vector of argument values
#'   and `params` is the atom's parameter list. `value()` returns a scalar,
#'   `grad()` returns a numeric vector of length `length(args)`, and
#'   `hess()` returns a `length(args) x length(args)` numeric matrix.
#' @param dcp_info     Function `(args_sign, params)` returning a list with
#'   components `curvature` (one of `"affine"`, `"convex"`, `"concave"`,
#'   `"unknown"`), `monotonicity` (a character vector of length equal to
#'   the arity, with entries `"nondecreasing"`, `"nonincreasing"`, or
#'   `"unknown"`), and `sign` (one of `"positive"`, `"nonneg"`,
#'   `"negative"`, `"nonpos"`, `"zero"`, `"unknown"`). Used by
#'   [curvature()] and the Phase 3 minorization engine.
#'
#' @return Invisibly `NULL`. Called for its side effect on the registry.
#' @examples
#' \dontrun{
#' # Register a custom 'tanh' atom (advanced; rarely needed):
#' register_atom("tanh",
#'   arity = 1L,
#'   value = function(args, params) tanh(args[1]),
#'   grad  = function(args, params) 1 - tanh(args[1])^2,
#'   hess  = function(args, params) matrix(-2 * tanh(args[1]) * (1 - tanh(args[1])^2), 1, 1))
#' }
#' @keywords internal
#' @export
register_atom <- function(name, arity, value, grad, hess, dcp_info = NULL) {
  if (!is.character(name) || length(name) != 1L) {
    stop("register_atom(name, ...): name must be a single character string.")
  }
  .mmad_atoms[[name]] <- list(
    arity    = arity,
    value    = value,
    grad     = grad,
    hess     = hess,
    dcp_info = dcp_info
  )
  invisible(NULL)
}

#' Names of currently registered atoms
#' @return A character vector.
#' @examples
#' mmad_atom_names()
#' @keywords internal
#' @export
mmad_atom_names <- function() {
  ls(envir = .mmad_atoms, sorted = TRUE)
}

#' Look up an atom by name (internal)
#' @keywords internal
#' @noRd
mmad_atom <- function(name) {
  a <- .mmad_atoms[[name]]
  if (is.null(a)) {
    stop(sprintf("Unknown atom '%s'.", name))
  }
  a
}

# ---- Built-in atoms -------------------------------------------------------
# The dcp_info functions reference helpers (is_pos, is_nonneg, flip_sign,
# scale_sign, reduce_signs) that live in R/dcp.R. R's lazy binding means
# this is fine: those helpers exist by the time any dcp_info is invoked,
# even though atoms.R is sourced before dcp.R alphabetically.

# Variadic addition: f(x_1, ..., x_k) = sum_j x_j.
register_atom("add",
  arity = NA_integer_,
  value = function(args, params) sum(args),
  grad  = function(args, params) rep(1, length(args)),
  hess  = function(args, params) matrix(0, length(args), length(args)),
  dcp_info = function(args_sign, params) {
    list(curvature    = "affine",
         monotonicity = rep("nondecreasing", length(args_sign)),
         sign         = reduce_signs(args_sign))
  }
)

# Unary negation: f(x) = -x.
register_atom("neg",
  arity = 1L,
  value = function(args, params) -args[1],
  grad  = function(args, params) -1,
  hess  = function(args, params) matrix(0, 1, 1),
  dcp_info = function(args_sign, params) {
    list(curvature    = "affine",
         monotonicity = "nonincreasing",
         sign         = flip_sign(args_sign[1]))
  }
)

# Unary scaling by a numeric parameter c: f(x; c) = c * x.
register_atom("scale",
  arity = 1L,
  value = function(args, params) params$c * args[1],
  grad  = function(args, params) params$c,
  hess  = function(args, params) matrix(0, 1, 1),
  dcp_info = function(args_sign, params) {
    cv <- params$c
    list(curvature    = "affine",
         monotonicity = if (cv >= 0) "nondecreasing" else "nonincreasing",
         sign         = scale_sign(cv, args_sign[1]))
  }
)

# Natural logarithm: f(x) = log(x), domain x > 0.
register_atom("log",
  arity = 1L,
  value = function(args, params) {
    x <- args[1]
    if (!is.finite(x) || x <= 0) {
      stop(sprintf("log(): argument must be positive (got %g).", x))
    }
    log(x)
  },
  grad  = function(args, params) 1 / args[1],
  hess  = function(args, params) matrix(-1 / args[1]^2, 1, 1),
  dcp_info = function(args_sign, params) {
    list(curvature    = "concave",
         monotonicity = "nondecreasing",
         sign         = "unknown")
  }
)

# Exponential: f(x) = exp(x).
register_atom("exp",
  arity = 1L,
  value = function(args, params) exp(args[1]),
  grad  = function(args, params) exp(args[1]),
  hess  = function(args, params) matrix(exp(args[1]), 1, 1),
  dcp_info = function(args_sign, params) {
    list(curvature    = "convex",
         monotonicity = "nondecreasing",
         sign         = "positive")
  }
)

# Power: f(x; c) = x^c. Curvature depends jointly on the exponent c and
# the sign of x. Decision table:
#
#   c == 0:                             affine, positive (constant 1)
#   c == 1:                             affine, sign passes through
#   c is positive even integer (>=2):   convex on R; mono = sign of x;
#                                         sign nonneg
#   c in (0, 1):                        concave + nondecreasing only when
#                                         x >= 0; otherwise unknown
#   c > 1, not even-int:                convex + nondecreasing only when
#                                         x >= 0; otherwise unknown
#   c < 0:                              convex + nonincreasing only when
#                                         x > 0; otherwise unknown
#
# The conservative "unknown" branches are deliberate: declaring a curvature
# we can't prove would feed garbage into the Phase 3 minorization engine.
register_atom("pow",
  arity = 1L,
  value = function(args, params) {
    x  <- args[1]
    cc <- params$c
    if (x == 0 && cc < 0) {
      stop(sprintf("(.)^%g: zero base with negative exponent.", cc))
    }
    if (x < 0 && cc != round(cc)) {
      stop(sprintf("(.)^%g: non-integer exponent applied to negative base %g.",
                   cc, x))
    }
    x^cc
  },
  grad = function(args, params) {
    params$c * args[1]^(params$c - 1)
  },
  hess = function(args, params) {
    matrix(params$c * (params$c - 1) * args[1]^(params$c - 2), 1, 1)
  },
  dcp_info = function(args_sign, params) {
    cc <- params$c
    s  <- args_sign[1]
    is_int <- (cc == round(cc))

    if (cc == 0) {
      return(list(curvature = "affine",
                  monotonicity = "nondecreasing",
                  sign = "positive"))
    }
    if (cc == 1) {
      return(list(curvature = "affine",
                  monotonicity = "nondecreasing",
                  sign = s))
    }
    if (is_int && cc > 1 && cc %% 2 == 0) {
      mono <- if (is_nonneg(s)) "nondecreasing"
              else if (is_nonpos(s)) "nonincreasing"
              else "unknown"
      return(list(curvature = "convex",
                  monotonicity = mono,
                  sign = "nonneg"))
    }
    if (cc > 0 && cc < 1) {
      if (is_nonneg(s)) {
        return(list(curvature = "concave",
                    monotonicity = "nondecreasing",
                    sign = "nonneg"))
      }
      return(list(curvature = "unknown",
                  monotonicity = "unknown",
                  sign = "unknown"))
    }
    if (cc > 1) {
      if (is_nonneg(s)) {
        return(list(curvature = "convex",
                    monotonicity = "nondecreasing",
                    sign = "nonneg"))
      }
      return(list(curvature = "unknown",
                  monotonicity = "unknown",
                  sign = "unknown"))
    }
    # cc < 0
    if (is_pos(s)) {
      return(list(curvature = "convex",
                  monotonicity = "nonincreasing",
                  sign = "positive"))
    }
    list(curvature = "unknown",
         monotonicity = "unknown",
         sign = "unknown")
  }
)
