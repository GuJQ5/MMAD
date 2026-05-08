## Diagnostic helper: a one-call sanity check on a target expression
## before launching mmad(). Reports overall curvature, the per-summand
## breakdown, domain feasibility at the initial point, and what
## minorize_at() will produce (fully separable vs. with a non_separable
## residue).

#' Diagnose an `mmad_expr` (or formula) before optimization
#'
#' Runs curvature inference, top-level summand decomposition, domain
#' check at `init`, and a dry-run of [minorize_at()]. Returns a
#' structured report with a `print()` method intended as a quick health
#' check before invoking [mmad()].
#'
#' @param expr A formula or `mmad_expr` representing the target.
#' @param init Numeric vector of initial parameter values (named or not).
#' @param data Optional `list`/`data.frame` of theta-free symbols
#'   referenced in `expr` (only relevant when `expr` is a formula).
#'
#' @return A list with class `"mmad_check"` and components
#'   `target_curvature`, `target_sign`, `is_dcp`, `summands`,
#'   `domain_ok`, `domain_message`, `is_separable`,
#'   `non_separable_indices`.
#' @examples
#' chk <- Function_check(~ log(theta[1] + theta[2]) - theta[1],
#'                       init = c(1, 1))
#' print(chk)
#' @export
Function_check <- function(expr, init, data = NULL) {
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
  if (length(raw) == 0L) {
    summands <- data.frame(
      coef         = numeric(0),
      curvature    = character(0),
      sign         = character(0),
      theta_indices = character(0),
      stringsAsFactors = FALSE)
  } else {
    summands <- data.frame(
      coef         = vapply(raw, function(r) r$coef, numeric(1)),
      curvature    = vapply(raw, function(r) curvature(r$expr), character(1)),
      sign         = vapply(raw, function(r) sign_of(r$expr),  character(1)),
      theta_indices = vapply(raw, function(r) {
        idx <- theta_indices(r$expr)
        if (length(idx) == 0L) "(none)"
        else paste(idx, collapse = ",")
      }, character(1)),
      stringsAsFactors = FALSE)
  }

  # Domain check: just try evaluating at init.
  domain_ok      <- TRUE
  domain_message <- "ok"
  ev <- tryCatch(evaluate_expr(expr, init_vec),
                 error = function(e) e)
  if (inherits(ev, "error")) {
    domain_ok      <- FALSE
    domain_message <- conditionMessage(ev)
  }

  # Dry-run minorization.
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

  structure(list(
    target_curvature      = target_curv,
    target_sign           = target_sign,
    is_dcp                = dcp_overall,
    summands              = summands,
    domain_ok             = domain_ok,
    domain_message        = domain_message,
    is_separable          = is_sep,
    non_separable_indices = ns_idx
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
print.mmad_check <- function(x, ...) {
  cat("MMAD function check\n")
  cat(sprintf("  Curvature:     %s\n", x$target_curvature))
  cat(sprintf("  Sign:          %s\n", x$target_sign))
  cat(sprintf("  DCP-valid:     %s\n", if (x$is_dcp) "yes" else "no"))
  cat(sprintf("  Domain at init: %s\n",
              if (x$domain_ok) "ok" else paste0("FAIL (", x$domain_message, ")")))
  if (!is.na(x$is_separable)) {
    cat(sprintf("  Surrogate:     %s\n",
                if (x$is_separable) "fully separable"
                else sprintf("non-separable on theta indices %s",
                             paste(x$non_separable_indices, collapse = ", "))))
  }
  if (nrow(x$summands) > 0L) {
    cat(sprintf("  Top-level summands (%d):\n", nrow(x$summands)))
    print(x$summands, row.names = FALSE)
  }
  invisible(x)
}
