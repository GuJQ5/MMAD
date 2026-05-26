## Diagnostic helper: a one-call sanity check on a target expression
## before launching mmad(). Reports overall curvature, the per-summand
## breakdown, domain feasibility at the initial point, and what
## minorize_at() will produce (fully separable vs. with a non_separable
## residue).
##
## NOTE on DCP vs. minorizability: a summand (or the whole expression)
## can be fully minorizable without satisfying the strict DCP composition
## rules. For example, log(1 + exp(theta)) has "unknown" DCP curvature
## but is handled by the additive-Jensen rule, and exp(log(theta1) +
## log(theta2)) is handled by the hyperplane-on-concave-inner rule. The
## report therefore separates the two concepts:
##   is_dcp        -- whether the whole expression is DCP-valid (strict)
##   is_separable  -- whether minorize_at() produces a fully separable
##                    surrogate (the operationally relevant criterion)
## The per-summand table likewise carries both a `curvature` column
## (DCP inference result) and a `minorizable` column (outcome of a
## dry-run minorize_at() on that summand alone at `init`).
##
## The print() method also renders a full expression tree in which every
## node is annotated with [curvature | sign | minorizable?]. Internal
## structural nodes (add/neg/scale) are labelled with their role;
## atom-call nodes show the atom name and any parameters; leaf nodes
## show theta[i] or the constant value.

#' Diagnose an `mmad_expr` (or formula) before optimization
#'
#' Runs curvature inference, top-level summand decomposition, domain
#' check at `init`, and a dry-run of [minorize_at()]. Returns a
#' structured report with a `print()` method intended as a quick health
#' check before invoking [mmad()].
#'
#' The `is_dcp` flag reflects strict DCP composition rules and may be
#' `FALSE` even for expressions that are fully minorizable by MMAD's
#' extended rules (e.g. `log(1 + exp(theta))` or
#' `exp(log(theta1) + log(theta2))`). The operationally relevant
#' criterion is `is_separable`: if that is `TRUE`, [mmad()] will
#' optimise the expression successfully.
#'
#' The `print()` method renders a full expression tree in which every
#' node is annotated with `[curvature | sign | minorizable?]`.
#'
#' @param expr A formula or `mmad_expr` representing the target.
#' @param init Numeric vector of initial parameter values (named or not).
#' @param data Optional `list`/`data.frame` of theta-free symbols
#'   referenced in `expr` (only relevant when `expr` is a formula).
#'
#' @param tree Logical; if `TRUE`, `print()` will also render the full
#'   annotated expression tree showing the minorizability of every node.
#'   Defaults to `FALSE`.  The tree data is always built and stored in
#'   `$expr_tree` regardless of this flag; the flag only affects printing.
#'
#' @return A list with class `"mmad_check"` and components
#'   `target_curvature`, `target_sign`, `is_dcp`, `is_minorizable`,
#'   `summands`, `domain_ok`, `domain_message`, `is_separable`,
#'   `non_separable_indices`, `expr_tree`, `show_tree`.
#' @examples
#' chk <- Function_check(~ log(theta[1] + theta[2]) - theta[1],
#'                       init = c(1, 1))
#' print(chk)
#' print(chk, tree = TRUE)   # show the expression tree
#' @export
Function_check <- function(expr, init, data = NULL, tree = FALSE) {
  if (inherits(expr, "formula")) {
    expr <- as_mmad_expr(expr, init = init, data = data)
  }
  if (!inherits(expr, "mmad_expr")) {
    stop("Function_check(): expr must be a formula or an mmad_expr.")
  }
  init_vec <- as.numeric(init)

  # Whole-expression diagnostics.
  target_curv <- curvature(expr)
  target_sign <- sign_of(expr)
  dcp_overall <- target_curv %in% c("affine", "convex", "concave")

  # Per-summand decomposition (top-level additive structure only).
  raw <- collect_summands(expr)

  # Domain check: just try evaluating at init.
  domain_ok      <- TRUE
  domain_message <- "ok"
  ev <- tryCatch(evaluate_expr(expr, init_vec),
                 error = function(e) e)
  if (inherits(ev, "error")) {
    domain_ok      <- FALSE
    domain_message <- conditionMessage(ev)
  }

  # Dry-run minorization of the whole expression.
  if (domain_ok) {
    surr <- tryCatch(minorize_at(expr, init_vec),
                     error = function(e) e)
    if (inherits(surr, "error")) {
      is_sep <- NA
      ns_idx <- NA_integer_
      domain_message <- paste0(
        "minorize_at() raised an error: ", conditionMessage(surr))
    } else {
      is_sep <- surr$is_separable
      ns_idx <- if (!is.null(surr$non_separable))
        theta_indices(surr$non_separable) else integer(0)
    }
  } else {
    is_sep <- NA
    ns_idx <- NA_integer_
  }

  # is_minorizable: TRUE iff the whole expression is separably minorizable.
  # This is the operationally relevant flag: it can be TRUE even when
  # is_dcp is FALSE (e.g. log(1+exp(theta)), exp(log(t1)+log(t2))).
  is_minor <- isTRUE(is_sep)

  # Per-summand table. For each summand we record:
  #   coef, theta_indices, minorizable.
  # The `minorizable` column is obtained by a dry-run minorize_at() on
  # each summand's sub-expression scaled by its coefficient. This allows
  # the user to pinpoint exactly which summands resist minorization.
  if (length(raw) == 0L) {
    summands <- data.frame(
      coef          = numeric(0),
      theta_indices = character(0),
      minorizable   = character(0),
      stringsAsFactors = FALSE)
  } else {
    summand_minor <- vapply(raw, function(r) {
      if (!domain_ok) return("unknown")
      # Build the scaled summand as an mmad_expr and dry-run minorize_at.
      scaled <- scale_expr(r$coef, r$expr)
      res <- tryCatch(minorize_at(scaled, init_vec), error = function(e) e)
      if (inherits(res, "error"))  return("error")
      if (res$is_separable)        return("yes")
      return("no")
    }, character(1))

    summands <- data.frame(
      coef          = vapply(raw, function(r) r$coef,          numeric(1)),
      theta_indices = vapply(raw, function(r) {
        idx <- theta_indices(r$expr)
        if (length(idx) == 0L) "(none)"
        else paste(idx, collapse = ",")
      }, character(1)),
      minorizable   = summand_minor,
      stringsAsFactors = FALSE)
  }

  # Build the full annotated expression tree.
  expr_tree <- build_expr_tree(expr, init_vec, domain_ok)

  structure(list(
    target_curvature      = target_curv,
    target_sign           = target_sign,
    is_dcp                = dcp_overall,
    is_minorizable        = is_minor,
    summands              = summands,
    domain_ok             = domain_ok,
    domain_message        = domain_message,
    is_separable          = is_sep,
    non_separable_indices = ns_idx,
    expr_tree             = expr_tree,
    show_tree             = isTRUE(tree)
  ), class = "mmad_check")
}

# Walk the top-level additive structure (add / neg / scale only) and
# emit a list of (coef, sub_expr) pairs where each sub_expr is the leaf
# (any non-additive node). Used by Function_check for per-summand
# diagnostics.
collect_summands <- function(expr) {
  out <- list()
  walk <- function(node, c) {
    if (inherits(node, "mmad_call")) {
      if (node$op == "add") {
        for (a in node$args) walk(a, c)
        return(invisible(NULL))
      }
      if (node$op == "neg") {
        walk(node$args[[1L]], -c)
        return(invisible(NULL))
      }
      if (node$op == "scale") {
        walk(node$args[[1L]], c * node$params$c)
        return(invisible(NULL))
      }
    }
    out[[length(out) + 1L]] <<- list(coef = c, expr = node)
    invisible(NULL)
  }
  walk(expr, 1)
  out
}

#' @export
print.mmad_check <- function(x, tree = x$show_tree, ...) {
  cat("MMAD function check\n")
  cat(sprintf("  Curvature:      %s\n", x$target_curvature))
  cat(sprintf("  Sign:           %s\n", x$target_sign))
  cat(sprintf("  DCP-valid:      %s\n", if (x$is_dcp) "yes" else "no"))
  cat(sprintf("  Minorizable:    %s\n",
              if (is.na(x$is_separable)) "unknown"
              else if (x$is_minorizable) "yes (fully separable surrogate)"
              else "no (non-separable residue; see non_separable_indices)"))
  cat(sprintf("  Domain at init: %s\n",
              if (x$domain_ok) "ok"
              else paste0("FAIL (", x$domain_message, ")")))
  if (!is.na(x$is_separable) && !x$is_separable) {
    cat(sprintf("  Non-sep. theta: %s\n",
                paste(x$non_separable_indices, collapse = ", ")))
  }
  if (nrow(x$summands) > 0L) {
    cat(sprintf("  Top-level summands (%d):\n", nrow(x$summands)))
    print(x$summands, row.names = FALSE)
  }
  if (isTRUE(tree)) {
    cat("\nExpression tree  [minorizable?]\n")
    print_expr_tree(x$expr_tree)
  }
  invisible(x)
}

# ---- Expression-tree annotation -------------------------------------------

# Recursively annotate every node of an mmad_expr with:
#   label       - short human-readable description of the node itself
#   curvature   - result of infer_dcp()$curvature
#   sign        - result of infer_dcp()$sign
#   minorizable - "yes" / "no" / "trivial" / "unknown" (dry-run minorize_at)
#   children    - list of annotated child trees (empty for leaves)
#
# `init_vec`  is needed for the minorizability dry-run.
# `domain_ok` suppresses the dry-run when the domain check already failed.
build_expr_tree <- function(expr, init_vec, domain_ok = TRUE) {
  dcp   <- infer_dcp(expr)
  curv  <- dcp$curvature
  sgn   <- dcp$sign

  # Minorizability: dry-run minorize_at on the sub-expression.
  # For affine or univariate nodes we skip the dry-run (trivially handled).
  idx <- theta_indices(expr)
  if (curv == "affine" || length(idx) <= 1L) {
    minor <- "trivial"
  } else if (!domain_ok) {
    minor <- "unknown"
  } else {
    res <- tryCatch(minorize_at(expr, init_vec), error = function(e) e)
    if (inherits(res, "error")) {
      minor <- "error"
    } else if (res$is_separable) {
      minor <- "yes"
    } else {
      minor <- "no"
    }
  }

  # Build the node label and recurse into children.
  if (inherits(expr, "mmad_var")) {
    label    <- sprintf("theta[%d]", expr$index)
    children <- list()
  } else if (inherits(expr, "mmad_const")) {
    label    <- format(expr$value)
    children <- list()
  } else if (inherits(expr, "mmad_call")) {
    label <- switch(expr$op,
      "add"   = "(+) add",
      "neg"   = "(-) neg",
      "scale" = sprintf("(*) scale  c=%s", format(expr$params$c)),
      "pow"   = sprintf("pow  c=%s", format(expr$params$c)),
      expr$op   # all other atoms: log, exp, etc.
    )
    children <- lapply(expr$args, build_expr_tree,
                       init_vec = init_vec, domain_ok = domain_ok)
  } else {
    label    <- "<unknown>"
    children <- list()
  }

  list(label       = label,
       curvature   = curv,
       sign        = sgn,
       minorizable = minor,
       children    = children)
}

# Render an annotated tree produced by build_expr_tree().
# Uses plain ASCII box-drawing: the last child gets "`--", others "|--",
# and their subtrees are indented with "|   " / "    " respectively.
print_expr_tree <- function(node, prefix = "", is_last = TRUE) {
  connector <- if (is_last) "`-- " else "|-- "
  annotation <- sprintf("[%s]", node$minorizable)
  cat(prefix, connector, node$label, "  ", annotation, "\n", sep = "")

  child_prefix <- paste0(prefix, if (is_last) "    " else "|   ")
  n <- length(node$children)
  for (i in seq_len(n)) {
    print_expr_tree(node$children[[i]],
                    prefix  = child_prefix,
                    is_last = (i == n))
  }
  invisible(NULL)
}
