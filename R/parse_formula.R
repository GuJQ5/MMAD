## Phase 5: formula-based interface.
##
## Lowers a one-sided R formula like
##
##   ~ 12 * log(0.5 * theta[1] + 0.5 * theta[2]) - 6 * theta[1]
##
## into the same `mmad_expr` tree the user could have built by hand via
## the operator overloads. Two parameter conventions are supported:
##
##   * Indexed via `theta[i]` (matches the math notation in the paper);
##     `i` may be a literal integer or any expression that evaluates to
##     a positive integer in the formula's environment / `data`.
##
##   * Named parameters (e.g. `alpha`, `beta`) when `init` carries
##     names. The name -> index map is `names(init)`.
##
## Both can be mixed: `~ alpha + theta[2]` with `init = c(alpha = 1,
## beta = 0.5)` resolves both to the same tree as `~ theta[1] +
## theta[2]`.
##
## Theta-free sub-expressions are evaluated eagerly in the formula's
## environment (or `data`, if supplied) and frozen as numeric constants.
## `sum(X)`, `mean(Y)`, `crossprod(...)`, anything -- as long as it does
## not reference a parameter.
##
## `log(1 - x)` needs no special syntax: the AST walker maps the literal
## R expression `1 - x` to `add(const(1), neg(...))` via the existing
## smart constructors, so `log(1 - 0.5 * theta[1] - 0.5 * theta[2])`
## parses to exactly the tree the user would build by hand.

#' Convert a formula or expression into an `mmad_expr`
#'
#' @param x      A one-sided formula (`~ ...`) or an existing `mmad_expr`.
#' @param init   Numeric vector of initial parameter values. If named,
#'   the names supply the parameter vocabulary; otherwise parameters
#'   must be referenced via `theta[i]`.
#' @param data   Optional list / data frame in which to evaluate
#'   theta-free sub-expressions. Default: the formula's environment.
#'
#' @return An `mmad_expr`.
#' @examples
#' expr <- as_mmad_expr(~ log(theta[1] + theta[2]), init = c(0, 0))
#' evaluate_expr(expr, c(1, 2))$value
#' @export
as_mmad_expr <- function(x, init = NULL, data = NULL) {
  UseMethod("as_mmad_expr")
}

#' @export
as_mmad_expr.mmad_expr <- function(x, init = NULL, data = NULL) x

#' @export
as_mmad_expr.formula <- function(x, init = NULL, data = NULL) {
  if (length(x) != 2L) {
    stop("as_mmad_expr.formula(): formula must be one-sided (use `~ ...`).")
  }
  theta_names <- if (!is.null(init) &&
                     !is.null(names(init)) &&
                     all(nzchar(names(init)))) {
    names(init)
  } else {
    character(0L)
  }
  data_env <- if (is.null(data)) {
    environment(x)
  } else {
    list2env(as.list(data), parent = environment(x))
  }
  p <- if (is.null(init)) NA_integer_ else as.integer(length(init))
  parse_formula_node(x[[2L]], theta_names, data_env, p)
}

#' @export
as_mmad_expr.numeric <- function(x, init = NULL, data = NULL) {
  if (length(x) != 1L || !is.finite(x)) {
    stop("as_mmad_expr.numeric: only finite scalar numerics can be coerced.")
  }
  mmad_const(as.numeric(x))
}

#' @export
as_mmad_expr.default <- function(x, init = NULL, data = NULL) {
  stop("Cannot coerce object of class ",
       paste(class(x), collapse = "/"), " to mmad_expr.")
}

# ---- Internal walker -----------------------------------------------------

# Recursively map an R AST node to an mmad_expr. Errors point at the
# specific sub-expression that couldn't be handled.
parse_formula_node <- function(node, theta_names, data_env, p = NA_integer_) {
  # Numeric literal.
  if (is.numeric(node) && length(node) == 1L) {
    return(mmad_const(as.numeric(node)))
  }

  # Symbol: parameter or eager-evaluated constant.
  if (is.symbol(node)) {
    name <- as.character(node)
    if (name == "theta") {
      stop("Bare 'theta' is not a valid formula term outside `X %*% theta` ",
           "inside sum(). Use `theta[i]` for the i-th parameter, or named ",
           "parameters when `init` carries names.")
    }
    if (length(theta_names) > 0L && name %in% theta_names) {
      return(mmad_var(match(name, theta_names)))
    }
    val <- tryCatch(eval(node, envir = data_env),
                    error = function(e) NULL)
    if (is.numeric(val) && length(val) == 1L && is.finite(val)) {
      return(mmad_const(as.numeric(val)))
    }
    stop(sprintf(
      "Cannot resolve symbol '%s' in formula. Either declare it via init = c(%s = ...) or supply it through `data`.",
      name, name))
  }

  if (!is.call(node)) {
    stop(sprintf("Unrecognised AST node in formula: %s",
                 deparse(node, width.cutoff = 80L)))
  }

  op_node <- node[[1L]]
  if (!is.symbol(op_node)) {
    stop(sprintf("Unsupported call expression: %s",
                 deparse(node, width.cutoff = 80L)))
  }
  op   <- as.character(op_node)
  args <- as.list(node[-1L])

  # Parens.
  if (op == "(") {
    return(parse_formula_node(args[[1L]], theta_names, data_env, p))
  }

  # sum() over observations: row-wise expansion.
  if (op == "sum" && length(args) == 1L) {
    return(parse_sum_call(args[[1L]], theta_names, data_env, p))
  }

  # theta[i] indexing.
  if (op == "[" && length(args) == 2L &&
      is.symbol(args[[1L]]) && as.character(args[[1L]]) == "theta") {
    idx_val <- tryCatch(eval(args[[2L]], envir = data_env),
                        error = function(e) NULL)
    if (is.numeric(idx_val) && length(idx_val) == 1L &&
        is.finite(idx_val) && idx_val == round(idx_val) && idx_val > 0) {
      return(mmad_var(as.integer(idx_val)))
    }
    stop(sprintf(
      "theta[%s]: index must evaluate to a positive integer, got %s.",
      deparse(args[[2L]]), deparse(idx_val)))
  }

  # Arithmetic ops.
  if (op %in% c("+", "-", "*", "/", "^")) {
    if (length(args) == 1L) {
      a <- parse_formula_node(args[[1L]], theta_names, data_env, p)
      if (op == "-") return(-a)
      if (op == "+") return( a)
      stop(sprintf("Unary `%s` is not supported in formula.", op))
    }
    if (length(args) == 2L) {
      return(parse_binop(op, args[[1L]], args[[2L]],
                         theta_names, data_env, node, p))
    }
    stop(sprintf("Operator `%s` requires one or two arguments; got %d.",
                 op, length(args)))
  }

  # Math functions.
  if (op %in% c("log", "exp", "sqrt") && length(args) == 1L) {
    a <- parse_formula_node(args[[1L]], theta_names, data_env, p)
    if (op == "log")  return(log(a))
    if (op == "exp")  return(exp(a))
    if (op == "sqrt") return(sqrt(a))
  }

  # `%*%` outside of a sum() is not supported in v1.
  if (op == "%*%") {
    stop(sprintf(
      "`%s`: matrix-vector products with theta are only supported inside sum(). Wrap the expression in sum(...) to expand it row-by-row.",
      deparse(node, width.cutoff = 80L)))
  }

  # Theta-free fallback: eager evaluation.
  if (is_theta_free(node, theta_names)) {
    val <- tryCatch(eval(node, envir = data_env), error = function(e) NULL)
    if (is.numeric(val) && length(val) == 1L && is.finite(val)) {
      return(mmad_const(as.numeric(val)))
    }
  }

  stop(sprintf(
    "Unsupported expression in formula: %s. Supported operators: + - * / ^; supported functions: log, exp, sqrt, sum (row-wise); supported references: theta[i], named parameters, theta-free sub-expressions; supported reduction: sum() with X %%*%% theta inside.",
    deparse(node, width.cutoff = 80L)))
}

# Binary operator handling. Strategy:
#   1. If both sides are theta-free, attempt to eagerly evaluate the whole
#      sub-expression as a numeric constant; if that fails, fall through
#      so per-side parsing surfaces the right error message.
#   2. Otherwise parse each side recursively (which itself handles the
#      "theta-free but unresolvable symbol" case with a proper message).
#      Then apply the operator-specific DCP rules: + and - work for any
#      pair; *, /, ^ require a constant on at least one (or specifically
#      the right-hand) side, with operator-specific error messages.
parse_binop <- function(op, a_node, b_node, theta_names, data_env,
                        full_node, p = NA_integer_) {
  a_free <- is_theta_free(a_node, theta_names)
  b_free <- is_theta_free(b_node, theta_names)

  if (a_free && b_free) {
    val <- tryCatch(eval(full_node, envir = data_env),
                    error = function(e) NULL)
    if (is.numeric(val) && length(val) == 1L && is.finite(val)) {
      return(mmad_const(as.numeric(val)))
    }
    # Fall through; the per-side parse below will produce a precise error.
  }

  a <- parse_formula_node(a_node, theta_names, data_env, p)
  b <- parse_formula_node(b_node, theta_names, data_env, p)

  a_val <- if (inherits(a, "mmad_const")) a$value else NULL
  b_val <- if (inherits(b, "mmad_const")) b$value else NULL

  if (op == "+") {
    if (!is.null(a_val)) return(a_val + b)
    if (!is.null(b_val)) return(a + b_val)
    return(a + b)
  }
  if (op == "-") {
    if (!is.null(a_val)) return(a_val - b)
    if (!is.null(b_val)) return(a - b_val)
    return(a - b)
  }
  if (op == "*") {
    if (!is.null(a_val)) return(a_val * b)
    if (!is.null(b_val)) return(a * b_val)
    stop(sprintf(
      "`%s`: products of two parameter-dependent expressions are non-DCP. One side must be a constant.",
      deparse(full_node, width.cutoff = 80L)))
  }
  if (op == "/") {
    if (!is.null(b_val)) return(a / b_val)
    stop(sprintf(
      "`%s`: division by an expression is non-DCP. The divisor must be a constant.",
      deparse(full_node, width.cutoff = 80L)))
  }
  if (op == "^") {
    if (!is.null(b_val)) return(a ^ b_val)
    stop(sprintf(
      "`%s`: exponent must be a numeric constant; expression-exponent is non-DCP.",
      deparse(full_node, width.cutoff = 80L)))
  }
  stop(sprintf("Operator `%s` is not supported in formula.", op))
}

# ---- sum() over observations and X %*% theta -----------------------------
#
# When the parser sees sum(<inner>) where <inner> references vector data
# or contains an X %*% theta product, it expands the inner expression
# row-by-row: for each i = 1..n it builds a scalar AST by substituting
#
#   * vector data symbols `v` (length n) with their i-th element v[i],
#   * `X %*% theta` with the affine combination
#       X[i, 1] * theta[1] + X[i, 2] * theta[2] + ... + X[i, p] * theta[p],
#
# parses each substituted AST with the standard parser, and assembles the
# n scalar mmad_exprs into a single sum via mmad_call("add", ...).
#
# v1 restrictions:
#   * `%*%` is only valid inside sum(); outside, it raises an error.
#   * The right-hand side of `%*%` must be the literal symbol `theta`.
#   * The left-hand side of `%*%` must be a symbol that resolves to a
#     numeric matrix in `data` whose row count equals the row count of
#     any other vector data referenced in the same sum().
#   * The user must have specified `init` so that p = length(init) is
#     known.
#
# These constraints cover the canonical statistical patterns: Poisson,
# logistic, multinomial logit, Cox partial likelihood, etc.

parse_sum_call <- function(inner_node, theta_names, data_env, p) {
  # Theta-free sum() is just a numeric reduction -- compute it once and
  # freeze as a constant. This handles cases like `sum(X) * theta[1]`
  # where `sum(X)` is a coefficient computed from data.
  if (is_theta_free(inner_node, theta_names)) {
    val <- tryCatch(eval(call("sum", inner_node), envir = data_env),
                    error = function(e) NULL)
    if (is.numeric(val) && length(val) == 1L && is.finite(val)) {
      return(mmad_const(as.numeric(val)))
    }
    stop(sprintf(
      "Theta-free sum(%s) did not evaluate to a finite numeric scalar.",
      deparse(inner_node, width.cutoff = 80L)))
  }

  n <- detect_row_dim(inner_node, theta_names, data_env)
  if (is.null(n)) {
    # Theta-dependent but no vector data: sum() is a no-op on a scalar.
    return(parse_formula_node(inner_node, theta_names, data_env, p))
  }
  if (is.na(p)) {
    stop("sum() with X %*% theta requires init to be supplied so the ",
         "parameter dimension can be determined.")
  }
  terms <- vector("list", n)
  for (i in seq_len(n)) {
    sub_ast <- expand_row(inner_node, theta_names, data_env, p, i, n)
    terms[[i]] <- parse_formula_node(sub_ast, theta_names, data_env, p)
  }
  # Build a single n-ary add directly to avoid O(n^2) flattening cost.
  if (length(terms) == 1L) return(terms[[1L]])
  mmad_call("add", terms)
}

# Walk an AST collecting row-dimension hints from any vector data
# references and any X %*% theta calls. Returns a single integer n if
# all collected hints agree, NULL if no hints exist, and stops with an
# error if the hints conflict.
detect_row_dim <- function(node, theta_names, data_env) {
  ns <- collect_row_dims(node, theta_names, data_env)
  if (length(ns) == 0L) return(NULL)
  uns <- unique(ns)
  if (length(uns) > 1L) {
    stop(sprintf(
      "sum(): inconsistent vector lengths inside the sum (found %s).",
      paste(uns, collapse = ", ")))
  }
  uns[1]
}

collect_row_dims <- function(node, theta_names, data_env) {
  if (is.symbol(node)) {
    name <- as.character(node)
    if (name == "theta") return(integer(0))
    if (name %in% theta_names) return(integer(0))
    val <- tryCatch(eval(node, envir = data_env), error = function(e) NULL)
    if (is.matrix(val))                           return(nrow(val))
    if (is.numeric(val) && length(val) > 1L)      return(length(val))
    return(integer(0))
  }
  if (!is.call(node)) return(integer(0))
  op <- if (is.symbol(node[[1L]])) as.character(node[[1L]]) else ""

  # theta[i]: do not recurse into the index expression as a row-dim source.
  if (op == "[" && length(node) >= 2L &&
      is.symbol(node[[2L]]) && as.character(node[[2L]]) == "theta") {
    return(integer(0))
  }
  ns <- integer(0)
  if (length(node) >= 2L) {
    for (i in 2L:length(node)) {
      ns <- c(ns, collect_row_dims(node[[i]], theta_names, data_env))
    }
  }
  ns
}

# Build a row-i scalar AST by substituting vector data references with
# their i-th element and `X %*% theta` calls with the affine combo
# X[i, 1] * theta[1] + ... + X[i, p] * theta[p].
expand_row <- function(node, theta_names, data_env, p, i, n) {
  if (is.symbol(node)) {
    name <- as.character(node)
    if (name == "theta")            return(node)   # handled inside %*%
    if (name %in% theta_names)      return(node)
    val <- tryCatch(eval(node, envir = data_env), error = function(e) NULL)
    if (is.numeric(val) && length(val) == n) {
      return(val[i])
    }
    if (is.matrix(val) && nrow(val) == n) {
      # A bare matrix reference outside %*% is not meaningful row-by-row;
      # leave it for the parser to error on.
      return(node)
    }
    # Scalar in data, or unknown symbol: leave as-is for the parser.
    return(node)
  }

  if (!is.call(node)) return(node)
  op <- if (is.symbol(node[[1L]])) as.character(node[[1L]]) else ""

  # theta[i]: leave indexing as-is.
  if (op == "[" && length(node) >= 2L &&
      is.symbol(node[[2L]]) && as.character(node[[2L]]) == "theta") {
    return(node)
  }

  # X %*% theta -> sum_j X[i, j] * theta[j].
  if (op == "%*%" && length(node) == 3L) {
    a_node <- node[[2L]]; b_node <- node[[3L]]
    if (is.symbol(b_node) && as.character(b_node) == "theta" &&
        is.symbol(a_node)) {
      mat <- tryCatch(eval(a_node, envir = data_env),
                      error = function(e) NULL)
      if (is.matrix(mat) && nrow(mat) == n && ncol(mat) == p) {
        if (p < 1L) stop("`X %*% theta`: parameter dimension must be >= 1.")
        if (p == 1L) {
          return(bquote(.(mat[i, 1L]) * theta[1L]))
        }
        terms <- lapply(seq_len(p), function(k) {
          bquote(.(mat[i, k]) * theta[.(k)])
        })
        return(Reduce(function(x, y) bquote(.(x) + .(y)), terms))
      }
      stop(sprintf(
        "`%s %%*%% theta`: left operand must be a numeric matrix with %d columns and %d rows.",
        deparse(a_node), p, n))
    }
    stop(sprintf(
      "Unsupported `%%*%%` form in sum(): %s. Only `X %%*%% theta` (matrix times the parameter vector) is supported.",
      deparse(node, width.cutoff = 80L)))
  }

  # General recursion: rebuild call with substituted args.
  new_args <- lapply(as.list(node[-1L]), expand_row,
                     theta_names = theta_names,
                     data_env = data_env, p = p, i = i, n = n)
  as.call(c(list(node[[1L]]), new_args))
}

# Detect whether an AST sub-expression contains any reference to a
# parameter (a named symbol in `theta_names`, the bare symbol `theta`
# used as a parameter vector inside %*%, or `theta[...]`).
is_theta_free <- function(node, theta_names) {
  if (is.symbol(node)) {
    name <- as.character(node)
    if (name == "theta") return(FALSE)
    return(!(name %in% theta_names))
  }
  if (is.numeric(node) || is.character(node) || is.logical(node)) {
    return(TRUE)
  }
  if (is.call(node)) {
    op <- if (is.symbol(node[[1L]])) as.character(node[[1L]]) else ""
    if (op == "[" && length(node) >= 2L &&
        is.symbol(node[[2L]]) && as.character(node[[2L]]) == "theta") {
      return(FALSE)
    }
    if (length(node) >= 2L) {
      for (i in 2L:length(node)) {
        if (!is_theta_free(node[[i]], theta_names)) return(FALSE)
      }
    }
    return(TRUE)
  }
  TRUE
}
