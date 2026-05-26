## Phase 3: minorization rewrite engine.
##
## Builds a separable surrogate S(theta | theta_0) of the target expression
## at the current iterate theta_0. The surrogate satisfies the standard
## MM tangency conditions:
##
##   S(theta_0 | theta_0)      = expr(theta_0)
##   grad S(theta_0 | theta_0) = grad expr(theta_0)
##   S(theta | theta_0)       <= expr(theta)   for all theta in the domain
##
## Before the top-down walk, `minorize_at()` applies a bottom-up algebraic
## pre-simplification pass (`simplify_expr()` in expr.R) that rewrites
## log(exp(h)) -> h, exp(log(h)) -> h (when h > 0), and log(h^c) ->
## c*log(h) (when h > 0). This ensures that the minorization rules below
## see the structurally simplest equivalent expression.
##
## Construction proceeds top-down through the additive structure of the
## target (descending through `add`, `neg`, `scale`) and applies four
## rules at each leaf, in priority order:
##
##   (1) AFFINE: if the leaf is affine in theta, decompose it into a
##       per-coordinate slot plus a constant.
##   (2) UNIVARIATE: if the leaf depends on only one theta index, extract
##       it as-is into that coordinate's slot. This is always a valid
##       tangent lower bound (the term equals itself), and gives a
##       tighter surrogate than the supporting hyperplane when the
##       1-d shape is concave -- so we prefer it whenever possible.
##   (3) JENSEN / SUPPORTING HYPERPLANE: for a multivariate leaf, apply
##       Jensen's inequality if the effective curvature is concave and
##       the inner is an affine combination with all slot values nonneg
##       at theta_0; apply the supporting hyperplane if the effective
##       curvature is convex.
##   (4) FALLBACK: if none of the above works (mixed-sign Jensen slots,
##       concave with non-affine inner, unknown curvature), include the
##       leaf as-is in a `non_separable` bucket. The full surrogate
##       still evaluates correctly via [evaluate_expr()] -- so the user
##       can compute value/gradient/Hessian -- but it is not separable
##       across that piece. The `is_separable` flag in the return value
##       reports whether the fallback bucket is empty.

#' Build a surrogate of an `mmad_expr` at a current iterate
#'
#' Given a target `expr` and a current iterate `theta_0`, construct a
#' surrogate `S(theta | theta_0)` that lower-bounds `expr(theta)` and is
#' tangent to it at `theta_0`. The surrogate is constructed to be as
#' separable across the coordinates of theta as the structure permits;
#' anything that resists clean minorization is preserved as-is in a
#' fallback bucket so that value/gradient/Hessian computations on the
#' full surrogate remain available.
#'
#' Construction uses only Jensen's inequality and the supporting
#' hyperplane; no per-atom tight surrogates are applied.
#'
#' @param expr    An `mmad_expr`.
#' @param theta_0 A numeric vector giving the current iterate.
#'
#' @return A list with components
#'   \describe{
#'     \item{`expr`}{the full surrogate as an `mmad_expr`}
#'     \item{`constant`}{the additive constant of the surrogate}
#'     \item{`per_coord`}{a list of length `length(theta_0)` whose
#'       `j`-th entry is the 1-d surrogate term in `theta[j]`, or
#'       `NULL` if no term touches that coordinate}
#'     \item{`non_separable`}{an `mmad_expr` representing the
#'       multi-coordinate residue we could not minorize, or `NULL` if
#'       every leaf was minorized or extracted}
#'     \item{`is_separable`}{`TRUE` iff `non_separable` is `NULL`}
#'   }
#' @examples
#' expr <- log(0.5 * mmad_var(1) + 0.5 * mmad_var(2))
#' surr <- minorize_at(expr, c(2, 1))
#' surr$is_separable
#' evaluate_expr(surr$expr, c(2, 1))$value     # tangent at theta_0
#' @export
minorize_at <- function(expr, theta_0) {
  if (!inherits(expr, "mmad_expr")) {
    stop("minorize_at(): expr must be an mmad_expr.")
  }
  if (!is.numeric(theta_0) || length(theta_0) == 0L) {
    stop("minorize_at(): theta_0 must be a non-empty numeric vector.")
  }
  theta_0 <- as.numeric(theta_0)
  p       <- length(theta_0)

  # Pre-simplification pass: apply algebraic identities E3/E4/E5 bottom-up
  # so that the minorization rules work on the structurally simplest
  # equivalent expression (e.g. log(exp(h)) becomes h, log(h^c) becomes
  # c*log(h)). This does not change the value or gradient of the expression.
  expr <- simplify_expr(expr)

  state <- new.env(parent = emptyenv())
  state$constant      <- 0
  state$per_coord     <- vector("list", p)
  state$non_separable <- NULL

  process_term(expr, coef = 1, theta_0 = theta_0, p = p, state = state)

  # Reassemble the full surrogate: constant + sum of per-coord terms +
  # non-separable residue (if any).
  full <- mmad_const(state$constant)
  for (j in seq_len(p)) {
    if (!is.null(state$per_coord[[j]])) {
      full <- full + state$per_coord[[j]]
    }
  }
  if (!is.null(state$non_separable)) {
    full <- full + state$non_separable
  }

  list(expr          = full,
       constant      = state$constant,
       per_coord     = state$per_coord,
       non_separable = state$non_separable,
       is_separable  = is.null(state$non_separable))
}

# ---- Internal helpers ----------------------------------------------------

# Top-down recursion. Walks through additive structure (add/neg/scale),
# accumulates a cumulative coefficient, and dispatches each leaf to one
# of four rules: affine slot decomposition, univariate extraction,
# Jensen/hyperplane, or non-separable fallback.
process_term <- function(expr, coef, theta_0, p, state) {
  if (inherits(expr, "mmad_call")) {
    if (expr$op == "add") {
      for (a in expr$args) process_term(a, coef, theta_0, p, state)
      return(invisible(NULL))
    }
    if (expr$op == "neg") {
      process_term(expr$args[[1L]], -coef, theta_0, p, state)
      return(invisible(NULL))
    }
    if (expr$op == "scale") {
      process_term(expr$args[[1L]],
                   coef * expr$params$c, theta_0, p, state)
      return(invisible(NULL))
    }
  }

  # Rule 1: affine in theta -- decompose per-slot.
  aff <- extract_affine(expr)
  if (!is.null(aff)) {
    aff <- normalize_affine(aff)
    state$constant <- state$constant + coef * aff$constant
    for (k in seq_along(aff$indices)) {
      add_to_coord(state, aff$indices[k],
                   scale_expr(coef * aff$coefs[k],
                              mmad_var(aff$indices[k])))
    }
    return(invisible(NULL))
  }

  # Determine the set of theta indices this leaf depends on.
  idx <- theta_indices(expr)

  # Rule 2a: theta-independent expression -- evaluate once and add as
  # constant. (Rare in practice, but cheap to handle.)
  if (length(idx) == 0L) {
    val <- evaluate_expr(expr, theta_0)$value
    state$constant <- state$constant + coef * val
    return(invisible(NULL))
  }

  # Rule 2b: univariate -- extract as-is into the corresponding slot.
  # The term equals itself, so it is trivially a tangent lower bound.
  if (length(idx) == 1L) {
    add_to_coord(state, idx, scale_expr(coef, expr))
    return(invisible(NULL))
  }

  # Rule 3: multivariate. Try Jensen / supporting hyperplane / nested.
  if (inherits(expr, "mmad_call") && length(expr$args) == 1L) {
    inner          <- expr$args[[1L]]
    inner_aff      <- extract_affine(inner)
    atom_call_curv <- curvature(expr)
    eff            <- effective_curvature(atom_call_curv, coef)

    if (eff == "concave" && !is.null(inner_aff)) {
      inner_aff_n <- normalize_affine(inner_aff)
      if (jensen_check(inner_aff_n, theta_0)) {
        apply_jensen(expr$op, expr$params, inner_aff_n,
                     coef, theta_0, p, state)
        return(invisible(NULL))
      }
    }
    if (eff == "concave" && is.null(inner_aff)) {
      # Nested Jensen: inner is concave, so we can recursively minorize it
      # into a separable lower bound and apply Jensen to g(separable_LB).
      if (try_nested_jensen(expr$op, expr$params, inner,
                            coef, theta_0, p, state)) {
        return(invisible(NULL))
      }
      # Additive Jensen (extended rule E1): inner is a flat additive sum of
      # nonneg univariate/constant terms, e.g. log(1 + exp(theta_j)).
      # dcp.R already classifies such expressions as "concave" via E1, so
      # eff is correctly "concave" here; try_nested_jensen fails because
      # the inner is convex, so we fall through to this rule.
      if (try_additive_jensen(expr$op, expr$params, inner,
                              coef, theta_0, p, state)) {
        return(invisible(NULL))
      }
    }
    if (eff == "convex") {
      # For a convex composite whose inner is concave (extended rule E2,
      # e.g. exp(log(t1)+log(t2))), try_hyperplane_concave_inner gives a
      # tighter separable surrogate by linearising g and then minorizing h.
      # Fall back to the full-expression hyperplane otherwise.
      if (curvature(inner) == "concave") {
        if (try_hyperplane_concave_inner(expr$op, expr$params, inner,
                                         coef, theta_0, p, state)) {
          return(invisible(NULL))
        }
      }
      apply_hyperplane(expr, coef, theta_0, p, state)
      return(invisible(NULL))
    }
  }

  # Rule 4: fallback -- include as-is in the non-separable bucket.
  add_to_non_separable(state, scale_expr(coef, expr))
  invisible(NULL)
}

# Curvature of (coef * f) given f's curvature and the sign of coef.
# Multiplying a concave/convex function by a negative scalar swaps the
# curvature; affine and unknown pass through.
effective_curvature <- function(atom_curv, coef) {
  if (atom_curv == "affine")  return("affine")
  if (atom_curv == "unknown") return("unknown")
  if (coef >= 0) atom_curv
  else if (atom_curv == "concave") "convex"
  else if (atom_curv == "convex")  "concave"
  else "unknown"
}

# Add `term` (an mmad_expr) into per_coord[[j]], building or extending
# the coordinate's symbolic surrogate.
add_to_coord <- function(state, j, term) {
  if (is.null(state$per_coord[[j]])) {
    state$per_coord[[j]] <- term
  } else {
    state$per_coord[[j]] <- state$per_coord[[j]] + term
  }
  invisible(NULL)
}

# Add `term` into the non-separable bucket.
add_to_non_separable <- function(state, term) {
  if (is.null(state$non_separable)) {
    state$non_separable <- term
  } else {
    state$non_separable <- state$non_separable + term
  }
  invisible(NULL)
}

# Set of theta indices (as a sorted unique integer vector) referenced
# anywhere in `expr`. Used to detect univariate sub-expressions.
theta_indices <- function(expr) {
  if (inherits(expr, "mmad_var"))   return(expr$index)
  if (inherits(expr, "mmad_const")) return(integer(0))
  if (inherits(expr, "mmad_call")) {
    out <- integer(0)
    for (a in expr$args) out <- c(out, theta_indices(a))
    return(sort(unique(out)))
  }
  integer(0)
}

# Return TRUE iff the affine combo's slot values at theta_0 (one per
# theta-slot, plus the constant slot if non-zero) are all nonneg and
# their sum is strictly positive. These are the preconditions for the
# vanilla Jensen step: weights nonneg and well-defined.
jensen_check <- function(inner_aff, theta_0) {
  k <- length(inner_aff$indices)
  slot_vals <- if (k > 0L) inner_aff$coefs * theta_0[inner_aff$indices]
               else numeric(0)
  b <- inner_aff$constant
  full <- if (b != 0) c(slot_vals, b) else slot_vals
  if (length(full) == 0L) return(FALSE)
  if (any(full < 0))      return(FALSE)
  if (sum(full) <= 0)     return(FALSE)
  TRUE
}

# Test whether `expr` is affine in theta. On success returns a list
# with `coefs`, `indices`, `constant`. On failure returns NULL.
extract_affine <- function(expr) {
  if (inherits(expr, "mmad_var")) {
    return(list(coefs = 1, indices = expr$index, constant = 0))
  }
  if (inherits(expr, "mmad_const")) {
    return(list(coefs    = numeric(0),
                indices  = integer(0),
                constant = expr$value))
  }
  if (inherits(expr, "mmad_call")) {
    if (expr$op == "add") {
      out <- list(coefs = numeric(0), indices = integer(0), constant = 0)
      for (a in expr$args) {
        sub <- extract_affine(a)
        if (is.null(sub)) return(NULL)
        out$coefs    <- c(out$coefs,    sub$coefs)
        out$indices  <- c(out$indices,  sub$indices)
        out$constant <- out$constant + sub$constant
      }
      return(out)
    }
    if (expr$op == "neg") {
      sub <- extract_affine(expr$args[[1L]])
      if (is.null(sub)) return(NULL)
      return(list(coefs    = -sub$coefs,
                  indices  =  sub$indices,
                  constant = -sub$constant))
    }
    if (expr$op == "scale") {
      sub <- extract_affine(expr$args[[1L]])
      if (is.null(sub)) return(NULL)
      cv <- expr$params$c
      return(list(coefs    = cv * sub$coefs,
                  indices  =  sub$indices,
                  constant = cv * sub$constant))
    }
    if (expr$op == "pow") {
      cc <- expr$params$c
      if (cc == 0) {
        return(list(coefs = numeric(0), indices = integer(0), constant = 1))
      }
      if (cc == 1) {
        return(extract_affine(expr$args[[1L]]))
      }
    }
  }
  NULL
}

# Aggregate duplicate theta indices into a single coefficient per index.
normalize_affine <- function(aff) {
  if (length(aff$indices) == 0L) return(aff)
  agg <- tapply(aff$coefs, aff$indices, sum)
  list(coefs    = unname(as.numeric(agg)),
       indices  = as.integer(names(agg)),
       constant = aff$constant)
}

# Jensen step. Caller has already validated the preconditions via
# jensen_check(); we just emit the per-coordinate Jensen terms.
apply_jensen <- function(op, params, inner_aff, coef, theta_0, p, state) {
  k <- length(inner_aff$indices)
  slot_vals <- if (k > 0L) inner_aff$coefs * theta_0[inner_aff$indices]
               else numeric(0)
  b  <- inner_aff$constant
  S  <- sum(slot_vals) + b
  atom_value_fn <- mmad_atom(op)$value

  # Per-theta-slot contributions:
  #   (slot_val_k / S) * coef * f((S / theta_0[i_k]) * theta[i_k])
  for (kk in seq_len(k)) {
    if (slot_vals[kk] == 0) next
    i_k <- inner_aff$indices[kk]
    if (theta_0[i_k] == 0) next
    w_k        <- slot_vals[kk] / S
    factor_in  <- S / theta_0[i_k]
    inner_expr <- scale_expr(factor_in, mmad_var(i_k))
    f_call     <- mmad_call(op, list(inner_expr), params)
    add_to_coord(state, i_k, scale_expr(w_k * coef, f_call))
  }

  # Constant-slot contribution: (b / S) * coef * f(S).
  if (b != 0) {
    f_at_S <- atom_value_fn(c(S), params)
    state$constant <- state$constant + (b / S) * coef * f_at_S
  }
}

# Supporting hyperplane step. Linearizes the sub-expression at theta_0
# and writes the result into the surrogate state.
apply_hyperplane <- function(expr, coef, theta_0, p, state) {
  res <- evaluate_expr(expr, theta_0)
  V <- res$value
  G <- res$gradient
  state$constant <- state$constant + coef * (V - sum(G * theta_0))
  for (j in seq_len(p)) {
    if (G[j] != 0) {
      add_to_coord(state, j, scale_expr(coef * G[j], mmad_var(j)))
    }
  }
}

# Additive-inner Jensen for `coef * g(inner)` where g is a concave
# nondecreasing 1-arg atom and `inner` is an additive sum whose terms are
# each a constant or a univariate (single-coordinate) sub-expression.
# This handles patterns like log(1 + exp(theta_j)) that fail the DCP
# composition rule (inner is convex, not affine) yet admit a clean
# Jensen decomposition by treating each additive piece as a slot.
#
# Preconditions checked internally:
#   - inner decomposes into a flat list of additive terms
#   - every term is either theta-independent (constant slot) or
#     depends on exactly one theta index (univariate slot)
#   - every slot's value at theta_0 is strictly positive
#
# Returns TRUE on success, FALSE if any precondition fails.
try_additive_jensen <- function(op, params, inner, coef, theta_0, p, state) {
  # Collect the flat additive terms of inner.
  terms <- collect_additive_terms(inner)
  if (is.null(terms)) return(FALSE)

  atom_value_fn <- mmad_atom(op)$value

  # For each term, determine its value at theta_0 and whether it is
  # univariate or constant.
  n_terms   <- length(terms)
  slot_vals <- numeric(n_terms)
  slot_idx  <- integer(n_terms)    # 0 = constant, >0 = theta index

  for (kk in seq_len(n_terms)) {
    t_idx <- theta_indices(terms[[kk]])
    if (length(t_idx) > 1L) return(FALSE)   # multi-coordinate piece
    v <- evaluate_expr(terms[[kk]], theta_0)$value
    if (!is.finite(v) || v <= 0) return(FALSE)   # must be strictly positive
    slot_vals[kk] <- v
    slot_idx[kk]  <- if (length(t_idx) == 0L) 0L else t_idx[1L]
  }

  S <- sum(slot_vals)
  if (S <= 0) return(FALSE)

  # Emit one Jensen term per slot.
  for (kk in seq_len(n_terms)) {
    w_k   <- slot_vals[kk] / S
    if (w_k == 0) next
    j     <- slot_idx[kk]
    if (j == 0L) {
      # Constant slot: contributes a numeric constant to the surrogate.
      state$constant <- state$constant +
        w_k * coef * atom_value_fn(c(S), params)
    } else {
      # Univariate slot in theta[j].
      factor_in  <- S / slot_vals[kk]   # scale so that g(factor * s_j(theta_j))
      inner_expr <- scale_expr(factor_in, terms[[kk]])
      g_call     <- mmad_call(op, list(inner_expr), params)
      add_to_coord(state, j, scale_expr(w_k * coef, g_call))
    }
  }
  TRUE
}

# Flatten an mmad_expr into its additive terms (accounting for neg/scale
# wrappers so that, e.g., a - b returns list(a, neg(b))).  Returns NULL
# if the expression contains an additive term that is itself an `add`
# node after unwrapping -- which shouldn't happen given the smart
# constructor -- or if we detect a non-unary structure we don't handle.
collect_additive_terms <- function(expr) {
  if (inherits(expr, "mmad_call") && expr$op == "add") {
    out <- list()
    for (a in expr$args) {
      sub <- collect_additive_terms(a)
      if (is.null(sub)) return(NULL)
      out <- c(out, sub)
    }
    return(out)
  }
  # A single non-add node is itself one term.
  list(expr)
}

# Supporting-hyperplane-on-g + minorize-inner construction for
# `coef * g(inner)` where g is a convex nondecreasing 1-arg atom and
# `inner` is a concave multivariate sub-expression (DCP composition fails
# because concave-inside-convex-nondecreasing is not a recognised DCP rule).
#
# Construction (two steps):
#   1. Linearise g at u0 = h(theta_0) via its supporting hyperplane:
#        g(u) >= g(u0) + g'(u0) * (u - u0)   for all u  (g convex)
#      Substituting u = h(theta):
#        g(h(theta)) >= g(u0) + g'(u0) * (h(theta) - u0)
#   2. g'(u0) >= 0 because g is nondecreasing.  h(theta) is concave and
#      multivariate, so minorize it recursively to a separable lower bound
#      h_hat(theta) = c + sum_j s_j(theta_j) <= h(theta).  Substituting:
#        g(h(theta)) >= g(u0) - g'(u0)*u0 + g'(u0)*h_hat(theta)
#      which is separable (affine in each per-coord piece of h_hat).
#
# When coef < 0 the effective curvature of coef*g is concave and the
# effective monotonicity flips to nonincreasing; in that case we instead
# require g'(u0) <= 0, which holds because the effective sign of the
# derivative is coef * g'(u0) and coef < 0 makes it nonpositive.  The
# algebra is the same; only the sign of the coefficient of h_hat changes.
#
# Returns TRUE on success, FALSE if any precondition fails (the caller
# then falls through to the non-separable bucket).
try_hyperplane_concave_inner <- function(op, params, inner,
                                         coef, theta_0, p, state) {
  # Step 1: evaluate g and g' at u0 = h(theta_0).
  inner_res <- evaluate_expr(inner, theta_0)
  u0        <- inner_res$value
  if (!is.finite(u0)) return(FALSE)

  atom      <- mmad_atom(op)
  g_u0      <- atom$value(c(u0), params)
  g_prime   <- atom$grad(c(u0), params)[1L]   # scalar g'(u0)
  if (!is.finite(g_u0) || !is.finite(g_prime)) return(FALSE)

  # g'(u0) must be nonneg (nondecreasing atom, coef > 0) or nonpos
  # (nonincreasing effective direction, coef < 0).  In either case
  # coef * g'(u0) >= 0 is required so that the coefficient of h_hat
  # in the surrogate is nonnegative and the inequality is preserved.
  eff_deriv <- coef * g_prime
  if (eff_deriv < 0) return(FALSE)   # would flip inequality when substituting h_hat

  # Step 2: recursively minorize the concave inner h.
  inner_surr <- minorize_at(inner, theta_0)
  if (!inner_surr$is_separable) return(FALSE)

  # Constant part of the surrogate:
  #   coef * [g(u0) - g'(u0) * u0]  +  coef * g'(u0) * inner_surr$constant
  state$constant <- state$constant +
    coef * (g_u0 - g_prime * u0) +
    eff_deriv * inner_surr$constant

  # Per-coordinate parts: coef * g'(u0) * s_j(theta_j)
  for (j in seq_len(p)) {
    if (!is.null(inner_surr$per_coord[[j]])) {
      add_to_coord(state, j,
                   scale_expr(eff_deriv, inner_surr$per_coord[[j]]))
    }
  }
  TRUE
}

# Recursive ("nested") Jensen for `coef * g(inner)` where g is a 1-arg
# concave atom and `inner` is itself a non-affine but DCP-concave
# sub-expression. The construction is:
#
#   1. Recursively minorize `inner` to a tangent lower bound h_hat that
#      is *separable*: h_hat(theta) = c + sum_j s_j(theta_j), with each
#      s_j a 1-d concave function in theta_j.
#   2. Because g is nondecreasing (in the effective sense after
#      multiplying by coef), h_hat <= h implies that
#      coef * g(h_hat) <= coef * g(h). So coef * g(h_hat) is a lower
#      bound on coef * g(h) at every theta in the domain.
#   3. Apply Jensen to coef * g(h_hat) treating each s_j(theta_j_0) and
#      the constant c as Jensen "slots":
#
#        coef * g(h_hat(theta)) = coef * g(sum_k w_k * x_k(theta))
#                              >= coef * sum_k w_k * g(x_k(theta))
#
#      where w_j = s_j(theta_j_0) / S, x_j = (S / s_j(theta_j_0)) * s_j,
#      and analogously for the constant slot.
#
# Returns TRUE on success (the caller should not also fall through to
# the non-separable bucket), FALSE if any precondition fails (the
# caller falls through). Failures: inner curvature isn't concave,
# effective monotonicity isn't nondecreasing, the recursive minorize
# left a non_separable residue, or the Jensen slot validation rejects
# the slot values at theta_0.
try_nested_jensen <- function(op, params, inner, coef, theta_0, p, state) {
  if (curvature(inner) != "concave") return(FALSE)

  atom <- mmad_atom(op)
  if (is.null(atom$dcp_info)) return(FALSE)

  # Effective monotonicity of (coef * g) in inner. The chain
  # h_hat <= h => coef * g(h_hat) <= coef * g(h) requires this to be
  # "nondecreasing".
  inner_sign <- sign_of(inner)
  info       <- atom$dcp_info(c(inner_sign), params)
  atom_mono  <- info$monotonicity[1]
  eff_mono   <- if (coef >= 0) {
    atom_mono
  } else if (atom_mono == "nondecreasing") {
    "nonincreasing"
  } else if (atom_mono == "nonincreasing") {
    "nondecreasing"
  } else {
    "unknown"
  }
  if (eff_mono != "nondecreasing") return(FALSE)

  inner_surr <- minorize_at(inner, theta_0)
  if (!inner_surr$is_separable) return(FALSE)

  # Collect slot values (s_j(theta_j_0) and the constant c).
  s_vals  <- numeric(p)
  s_exprs <- vector("list", p)
  for (j in seq_len(p)) {
    if (is.null(inner_surr$per_coord[[j]])) next
    s_vals[j]  <- evaluate_expr(inner_surr$per_coord[[j]], theta_0)$value
    s_exprs[[j]] <- inner_surr$per_coord[[j]]
  }
  c_const <- inner_surr$constant

  active <- !vapply(s_exprs, is.null, logical(1))
  full_slots <- s_vals[active]
  if (c_const != 0) full_slots <- c(full_slots, c_const)
  if (length(full_slots) == 0L) return(FALSE)
  if (any(full_slots < 0))      return(FALSE)
  S <- sum(full_slots)
  if (S <= 0)                   return(FALSE)

  atom_value_fn <- atom$value

  for (j in seq_len(p)) {
    if (!active[j]) next
    if (s_vals[j] == 0) next                # zero weight slot
    w_j       <- s_vals[j] / S
    factor_in <- S / s_vals[j]
    scaled    <- scale_expr(factor_in, s_exprs[[j]])
    g_call    <- mmad_call(op, list(scaled), params)
    add_to_coord(state, j, scale_expr(w_j * coef, g_call))
  }

  if (c_const != 0) {
    state$constant <- state$constant +
      (c_const / S) * coef * atom_value_fn(c(S), params)
  }

  TRUE
}
