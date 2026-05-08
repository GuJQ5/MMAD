## Phase 1 of the v3 redesign: a symbolic expression tree that will replace
## the bespoke nested-list representation. No DCP curvature inference, no
## minorization rewriting -- just a typed tree, operator sugar for users,
## and enough structure for an evaluator and a legacy-list adapter to be
## built on top.
##
## Three node kinds, all sharing the parent class `mmad_expr`:
##   * mmad_var(i)        -- the i-th parameter, theta[i]
##   * mmad_const(value)  -- a numeric constant
##   * mmad_call(op, args, params) -- application of a registered atom
##
## Users typically don't construct nodes directly; the Ops/Math group
## generics let them write `2 * theta1 + log(theta1 + theta2)` and have
## the tree built for them.

#' Construct a parameter reference `theta[i]`
#'
#' @param i Positive integer index into the parameter vector.
#' @param sign Optional declaration of the sign of `theta[i]`. One of
#'   `"unknown"` (default), `"positive"`, `"nonneg"`, `"negative"`,
#'   `"nonpos"`, or `"zero"`. The DCP curvature inference (Phase 2) uses
#'   this to refine curvature for sign-dependent atoms such as `pow`.
#'
#' @return An object of class `mmad_expr` (subclass `mmad_var`).
#' @examples
#' mmad_var(1)
#' mmad_var(2, sign = "positive")
#' @export
mmad_var <- function(i, sign = "unknown") {
  if (!is.numeric(i) || length(i) != 1L || !is.finite(i) ||
      i <= 0 || i != round(i)) {
    stop("mmad_var(i): i must be a positive integer scalar.")
  }
  valid_signs <- c("unknown", "positive", "nonneg",
                   "negative", "nonpos", "zero")
  if (!is.character(sign) || length(sign) != 1L || !sign %in% valid_signs) {
    stop("mmad_var(i, sign): sign must be one of: ",
         paste(valid_signs, collapse = ", "), ".")
  }
  structure(list(index = as.integer(i), sign = sign),
            class = c("mmad_var", "mmad_expr"))
}

#' Construct a numeric-constant expression node
#'
#' @param value A finite numeric scalar.
#'
#' @return An object of class `mmad_expr` (subclass `mmad_const`).
#' @examples
#' mmad_const(2.5)
#' @export
mmad_const <- function(value) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("mmad_const(value): value must be a finite numeric scalar.")
  }
  structure(list(value = as.numeric(value)),
            class = c("mmad_const", "mmad_expr"))
}

#' Construct an atom-application expression node
#'
#' Internal constructor; not normally called by users. Validates that the
#' atom name is registered (see [mmad_atom_names()]) and that all child
#' expressions are themselves `mmad_expr`.
#'
#' @param op    Name of a registered atom.
#' @param args  List of `mmad_expr` children.
#' @param params Named list of additional numeric atom-parameters
#'   (e.g. `list(c = 0.5)` for the exponent of a power atom).
#'
#' @return An object of class `mmad_expr` (subclass `mmad_call`).
#' @keywords internal
#' @export
mmad_call <- function(op, args, params = list()) {
  if (!is.character(op) || length(op) != 1L) {
    stop("mmad_call(op, ...): op must be a single character string.")
  }
  if (!is.list(args) ||
      !all(vapply(args, inherits, logical(1), "mmad_expr"))) {
    stop("mmad_call(op, args): args must be a list of mmad_expr objects.")
  }
  if (is.null(.mmad_atoms[[op]])) {
    stop(sprintf("Unknown atom '%s'. Registered atoms: %s.",
                 op, paste(mmad_atom_names(), collapse = ", ")))
  }
  structure(list(op = op, args = args, params = params),
            class = c("mmad_call", "mmad_expr"))
}

# NOTE: the `as_mmad_expr` coercer used by Ops.mmad_expr lives in
# R/parse_formula.R as an S3 generic with methods for `formula`,
# `mmad_expr`, `numeric`, and `default`. The methods together cover both
# the formula-parsing entry point and the bare-numeric-to-mmad_const
# coercion that the Ops/Math group generics rely on.

# ---- Internal smart-constructors ------------------------------------------
# These do light algebraic simplification at build time: flatten nested
# additions, fold scale-of-scale, drop zero terms, etc. The simplifications
# are only those that don't change semantics and don't depend on curvature.

add_expr <- function(args) {
  flat <- list()
  for (a in args) {
    if (inherits(a, "mmad_call") && a$op == "add") {
      flat <- c(flat, a$args)
    } else if (inherits(a, "mmad_const") && a$value == 0) {
      # drop additive zero
    } else {
      flat <- c(flat, list(a))
    }
  }
  if (length(flat) == 0L) return(mmad_const(0))
  if (length(flat) == 1L) return(flat[[1L]])
  mmad_call("add", flat)
}

neg_expr <- function(x) {
  if (inherits(x, "mmad_const")) return(mmad_const(-x$value))
  if (inherits(x, "mmad_call") && x$op == "neg") return(x$args[[1L]])
  mmad_call("neg", list(x))
}

scale_expr <- function(c, x) {
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c)) {
    stop("Scalar coefficient required when scaling an mmad_expr.")
  }
  if (c == 0) return(mmad_const(0))
  if (c == 1) return(x)
  if (inherits(x, "mmad_const")) return(mmad_const(c * x$value))
  if (inherits(x, "mmad_call") && x$op == "scale") {
    return(scale_expr(c * x$params$c, x$args[[1L]]))
  }
  if (inherits(x, "mmad_call") && x$op == "neg") {
    return(scale_expr(-c, x$args[[1L]]))
  }
  mmad_call("scale", list(x), params = list(c = as.numeric(c)))
}

pow_expr <- function(x, c) {
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c)) {
    stop("Power exponent must be a finite numeric scalar.")
  }
  if (c == 0) return(mmad_const(1))
  if (c == 1) return(x)
  mmad_call("pow", list(x), params = list(c = as.numeric(c)))
}

# ---- Group generics for ergonomic syntax ----------------------------------

#' Arithmetic operators for `mmad_expr`
#'
#' Lets users write expressions naturally: `2 * theta1 - log(theta1 + theta2)`
#' becomes a properly-typed expression tree.
#'
#' Restrictions enforced at construction time:
#' * `expr * expr` is rejected (one side must be a numeric scalar) -- products
#'   of two parameter-dependent expressions are non-DCP.
#' * `numeric / expr` is rejected -- division by an expression is non-DCP.
#' * `expr ^ expr` is rejected -- the exponent must be a numeric scalar.
#'
#' @param e1,e2 An `mmad_expr` and/or numeric scalar.
#' @return An `mmad_expr`.
#' @export
Ops.mmad_expr <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "-") return(neg_expr(e1))
    if (.Generic == "+") return(e1)
    stop(sprintf("Unary operator '%s' is not supported on mmad_expr.",
                 .Generic))
  }
  switch(.Generic,
    "+" = add_expr(list(as_mmad_expr(e1), as_mmad_expr(e2))),
    "-" = add_expr(list(as_mmad_expr(e1),
                        neg_expr(as_mmad_expr(e2)))),
    "*" = {
      if (is.numeric(e1) && length(e1) == 1L) {
        scale_expr(e1, as_mmad_expr(e2))
      } else if (is.numeric(e2) && length(e2) == 1L) {
        scale_expr(e2, as_mmad_expr(e1))
      } else {
        stop("mmad_expr * mmad_expr is not supported: products of two ",
             "parameter-dependent expressions are non-DCP.")
      }
    },
    "/" = {
      if (!is.numeric(e2) || length(e2) != 1L) {
        stop("mmad_expr only supports division by a numeric scalar.")
      }
      scale_expr(1 / e2, as_mmad_expr(e1))
    },
    "^" = {
      if (!is.numeric(e2) || length(e2) != 1L) {
        stop("mmad_expr exponent must be a numeric scalar.")
      }
      pow_expr(as_mmad_expr(e1), e2)
    },
    stop(sprintf("Operator '%s' is not supported on mmad_expr.", .Generic))
  )
}

#' Mathematical functions for `mmad_expr`
#'
#' `log()`, `exp()`, and `sqrt()` (the latter as `(.)^(1/2)`).
#'
#' @param x An `mmad_expr`.
#' @param ... Unused.
#' @return An `mmad_expr`.
#' @export
Math.mmad_expr <- function(x, ...) {
  switch(.Generic,
    "log"  = mmad_call("log", list(x)),
    "exp"  = mmad_call("exp", list(x)),
    "sqrt" = pow_expr(x, 0.5),
    stop(sprintf("Math function '%s' is not supported on mmad_expr.",
                 .Generic))
  )
}

# ---- Pretty printing -------------------------------------------------------

#' @export
print.mmad_expr <- function(x, ...) {
  cat("<mmad_expr>\n  ", format_mmad_expr(x), "\n", sep = "")
  invisible(x)
}

# Small helper: turn a node into a one-line readable form. Used by print().
format_mmad_expr <- function(x) {
  if (inherits(x, "mmad_var"))   return(sprintf("theta[%d]", x$index))
  if (inherits(x, "mmad_const")) return(format(x$value))
  if (inherits(x, "mmad_call")) {
    inner <- vapply(x$args, format_mmad_expr, character(1))
    switch(x$op,
      "add"   = paste(inner, collapse = " + "),
      "neg"   = sprintf("-(%s)", inner[1]),
      "scale" = sprintf("%s*(%s)", format(x$params$c), inner[1]),
      "pow"   = sprintf("(%s)^%s", inner[1], format(x$params$c)),
      sprintf("%s(%s)", x$op, paste(inner, collapse = ", "))
    )
  } else {
    "<unrecognised>"
  }
}
