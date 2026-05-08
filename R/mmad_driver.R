## Phase 4: top-level MM driver. Maximizes an mmad_expr target function
## by repeatedly:
##   1. building a separable-as-possible surrogate at the current iterate
##      via minorize_at();
##   2. taking a 1-d Newton step on each per_coord[[j]] whose coordinate j
##      is not entangled by the non_separable bucket;
##   3. taking a multivariate Newton step on the entangled block, where
##      the block surrogate is sum_{j in entangled} per_coord[[j]] +
##      non_separable;
##   4. line-searching on the actual target to enforce monotonic ascent.
##
## This is the Option-2 ("block-coordinate Newton with the residue as one
## block") strategy. When non_separable is empty, the entangled block is
## empty and the iteration reduces to per-coordinate Newton on each
## per_coord. When non_separable is the whole expression, it reduces to
## ordinary multivariate Newton on the surrogate.
##
## Two safety nets defend against pathological surrogates (univariate-as-
## is convex pieces, or unknown-curvature residues):
##   * Whenever a Hessian fails to be negative definite at the current
##     iterate, the corresponding block falls back to a gradient-direction
##     step. This guarantees the proposed direction is at least uphill on
##     the surrogate.
##   * The line search runs on the actual target value, not the
##     surrogate. Even if the surrogate proposes a bad step, the iteration
##     refuses to commit to it.

#' Minorization-Maximization driver for an `mmad_expr` target
#'
#' Given a target expression `expr` and an initial parameter vector
#' `init`, this iteratively maximizes `expr` by constructing a
#' separable surrogate at the current iterate (via [minorize_at()]),
#' taking Newton steps on the separable and entangled blocks, and
#' line-searching on the target to enforce ascent.
#'
#' @param expr            An `mmad_expr` or a one-sided formula
#'   (`~ ...`). When a formula, it is lowered via [as_mmad_expr()] using
#'   `init` and `data` as the parameter vocabulary and data
#'   environment.
#' @param init            Numeric vector of initial parameter values.
#'   May be named, in which case the names define the parameter
#'   vocabulary (e.g. `init = c(alpha = 1, beta = 0.5)`) and the
#'   returned estimate carries those names.
#' @param data            Optional list / data frame of values for
#'   theta-free symbols referenced in the formula. Default: the
#'   formula's environment.
#' @param tol             Convergence tolerance for both the maximum
#'   absolute parameter change and the maximum absolute target gradient
#'   component.
#' @param max_iter        Maximum number of MM iterations (default 1000).
#' @param line_search_max Maximum number of step-halving steps per
#'   iteration's line search (default 30).
#' @param track_history   If `TRUE`, return a data frame of per-iteration
#'   diagnostics (iteration index, value, gradient norm, parameter
#'   change). Default `FALSE`.
#' @param verbose         If `TRUE`, print per-iteration diagnostics.
#'
#' @return A list with components
#'   \describe{
#'     \item{`estimate`}{the final parameter vector (named if `init` was named)}
#'     \item{`value`}{the target value at `estimate`}
#'     \item{`iterations`}{number of iterations actually run}
#'     \item{`converged`}{`TRUE` iff both convergence criteria were met}
#'     \item{`history`}{per-iteration diagnostics (data frame), or `NULL`}
#'     \item{`message`}{character: short status string}
#'   }
#' @examples
#' fit <- mmad(~ log(theta[1] + theta[2]) - theta[1] - theta[2] + 2,
#'             init = c(2, 2))
#' fit$estimate
#' \donttest{
#' # Poisson regression via sum() and X %*% theta:
#' set.seed(1); n <- 50; p <- 2
#' X <- matrix(rnorm(n * p, sd = 0.5), nrow = n)
#' y <- rpois(n, exp(X %*% c(0.3, -0.2)))
#' mmad(~ sum(y * (X %*% theta) - exp(X %*% theta)),
#'      init = rep(0, p), data = list(X = X, y = y))
#' }
#' @export
mmad <- function(expr, init, data = NULL, tol = 1e-6, max_iter = 1000,
                 line_search_max = 30, track_history = FALSE,
                 verbose = FALSE) {
  if (inherits(expr, "formula")) {
    expr <- as_mmad_expr(expr, init = init, data = data)
  }
  if (!inherits(expr, "mmad_expr")) {
    stop("mmad(): expr must be a one-sided formula or an mmad_expr.")
  }
  if (!is.numeric(init) || length(init) == 0L) {
    stop("mmad(): init must be a non-empty numeric vector.")
  }

  init_names <- names(init)
  theta      <- as.numeric(init)
  p          <- length(theta)

  ev0   <- evaluate_expr(expr, theta)
  val   <- ev0$value
  grad  <- ev0$gradient

  history     <- if (track_history) vector("list", 0L) else NULL
  converged   <- FALSE
  message_str <- ""
  n_completed <- 0L                # iterations that actually took a step

  for (iter in seq_len(max_iter)) {
    surr <- minorize_at(expr, theta)

    if (!is.null(surr$non_separable)) {
      entangled <- theta_indices(surr$non_separable)
    } else {
      entangled <- integer(0L)
    }
    separable <- setdiff(seq_len(p), entangled)

    delta <- numeric(p)

    # Per-coordinate Newton on each separable per_coord[[j]].
    for (j in separable) {
      pc <- surr$per_coord[[j]]
      if (is.null(pc)) next
      ev <- evaluate_expr(pc, theta)
      g_j <- ev$gradient[j]
      H_j <- ev$hessian[j, j]
      delta[j] <- newton_1d_step(g_j, H_j)
    }

    # Multivariate Newton on the entangled block.
    if (length(entangled) > 0L) {
      block <- surr$non_separable
      for (j in entangled) {
        if (!is.null(surr$per_coord[[j]])) {
          block <- block + surr$per_coord[[j]]
        }
      }
      ev_e <- evaluate_expr(block, theta)
      g_e  <- ev_e$gradient[entangled]
      H_e  <- ev_e$hessian[entangled, entangled, drop = FALSE]
      delta[entangled] <- newton_block_step(g_e, H_e)
    }

    # Line search on the target.
    ls <- line_search(expr, theta, delta, val, line_search_max)
    if (!ls$accepted) {
      # No improvement found. If the gradient at the current iterate is
      # already below the convergence tolerance, treat this as a clean
      # convergence (the line search has nothing to improve on at a
      # stationary point). Otherwise we are genuinely stuck -- report
      # that and exit.
      grad_here     <- evaluate_expr(expr, theta)$gradient
      grad_norm_here <- max(abs(grad_here))
      if (grad_norm_here < tol) {
        converged   <- TRUE
        message_str <- "converged (no further ascent; gradient near zero)"
      } else {
        message_str <- sprintf(
          "no ascent direction at theta = %s (grad norm = %g); ",
          paste(format(theta, digits = 6), collapse = ", "),
          grad_norm_here)
      }
      break
    }

    theta_new <- ls$theta
    val_new   <- ls$value
    n_completed <- iter

    # Diagnostics for convergence and history.
    grad_new       <- evaluate_expr(expr, theta_new)$gradient
    theta_change   <- max(abs(theta_new - theta))
    grad_norm_new  <- max(abs(grad_new))

    if (track_history) {
      history[[n_completed]] <- data.frame(
        iteration    = n_completed,
        value        = val_new,
        grad_norm    = grad_norm_new,
        theta_change = theta_change
      )
    }
    if (verbose) {
      cat(sprintf(
        "iter %4d  value=%-12.6g  grad_norm=%-10.4g  theta_change=%-10.4g\n",
        iter, val_new, grad_norm_new, theta_change))
    }

    theta <- theta_new
    val   <- val_new
    grad  <- grad_new

    if (theta_change < tol && grad_norm_new < tol) {
      converged   <- TRUE
      message_str <- "converged"
      break
    }
  }

  if (!converged && message_str == "") {
    message_str <- sprintf("max_iter (%d) reached without convergence", max_iter)
  }

  if (!is.null(init_names)) names(theta) <- init_names

  structure(list(
    estimate    = theta,
    value       = val,
    iterations  = n_completed,
    converged   = converged,
    history     = if (track_history) do.call(rbind, history) else NULL,
    message     = message_str,
    expr        = expr
  ), class = c("mmad_fit", "list"))
}

#' @export
print.mmad_fit <- function(x, ...) {
  status <- if (x$converged) sprintf("converged in %d iterations", x$iterations)
            else sprintf("not converged (%s)", x$message)
  cat(sprintf("MMAD fit: %s\n", status))
  cat("Estimate:\n")
  print(x$estimate, ...)
  cat(sprintf("Value: %s\n", format(x$value, digits = 6)))
  invisible(x)
}

#' @export
summary.mmad_fit <- function(object, ...) {
  grad      <- evaluate_expr(object$expr, as.numeric(object$estimate))$gradient
  grad_norm <- max(abs(grad))
  status    <- if (object$converged) "converged" else "not converged"
  hist_msg  <- if (is.null(object$history)) "(not tracked)"
               else sprintf("%d rows", nrow(object$history))
  out <- list(
    status     = status,
    message    = object$message,
    iterations = object$iterations,
    value      = object$value,
    grad_norm  = grad_norm,
    estimate   = object$estimate,
    history    = hist_msg
  )
  structure(out, class = c("summary.mmad_fit", "list"))
}

#' @export
print.summary.mmad_fit <- function(x, ...) {
  cat("MMAD optimization\n")
  cat(sprintf("  Status:        %s (%s)\n", x$status, x$message))
  cat(sprintf("  Iterations:    %d\n",      x$iterations))
  cat(sprintf("  Final value:   %s\n",      format(x$value,     digits = 8)))
  cat(sprintf("  Gradient norm: %s\n",      format(x$grad_norm, digits = 4)))
  cat("  Estimate:\n")
  print(x$estimate, ...)
  cat(sprintf("  History:       %s\n", x$history))
  invisible(x)
}

# ---- Internal step helpers -----------------------------------------------

# 1-d Newton step for maximization. If H is sufficiently negative (the
# 1-d surrogate is locally concave), Newton: -g/H. Otherwise (Hessian
# zero or positive -- a univariate-as-is convex piece), fall back to a
# gradient-direction step. Step size is normalized so that the line
# search has a reasonable starting magnitude even when the gradient is
# very large; the line search will refine it.
newton_1d_step <- function(g, H, eps = 1e-12) {
  if (H < -eps) return(-g / H)
  if (abs(g) < eps) return(0)
  sign(g) * min(abs(g), 1.0)
}

# Multivariate Newton step for maximization on the entangled block.
# Uses Cholesky on -H (since for maximization we need H negative
# definite). Falls back to gradient direction when -H is not positive
# definite or the Newton direction does not point uphill.
newton_block_step <- function(g, H, eps = 1e-10) {
  k <- length(g)
  if (k == 0L) return(numeric(0))

  # Try negative-definite Newton.
  delta_newton <- tryCatch({
    L <- chol(-H + eps * diag(k))
    solve_lt <- backsolve(L, forwardsolve(t(L), g), upper.tri = TRUE)
    solve_lt
  }, error = function(e) NULL)

  if (!is.null(delta_newton) && sum(g * delta_newton) > 0) {
    return(delta_newton)
  }

  # Fall back to a normalized gradient-direction step.
  gn <- sqrt(sum(g^2))
  if (gn < eps) return(numeric(k))
  g / max(gn, 1.0)
}

# Backtracking line search on the actual target. Halves alpha until
# either the target strictly improves over `val_old` or the budget is
# exhausted. Returns a list with `accepted` flag.
line_search <- function(expr, theta, delta, val_old, max_back) {
  if (all(delta == 0)) {
    return(list(accepted = FALSE, theta = theta, value = val_old))
  }
  alpha <- 1.0
  for (k in seq_len(max_back)) {
    theta_try <- theta + alpha * delta
    val_try <- tryCatch(evaluate_expr(expr, theta_try)$value,
                        error = function(e) NA_real_)
    if (is.finite(val_try) && val_try > val_old) {
      return(list(accepted = TRUE, theta = theta_try, value = val_try))
    }
    alpha <- alpha / 2
  }
  list(accepted = FALSE, theta = theta, value = val_old)
}
