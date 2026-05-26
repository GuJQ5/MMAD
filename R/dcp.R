## Phase 2: DCP-style curvature, sign, and monotonicity inference for
## mmad_expr. Walks the tree bottom-up, querying each atom's dcp_info()
## hook and applying the standard composition rule:
##
##   h(g_1, ..., g_k) is convex iff h is convex AND for each j:
##     g_j is affine, OR
##     (h is nondecreasing in slot j AND g_j is convex), OR
##     (h is nonincreasing in slot j AND g_j is concave).
##
##   h is concave iff swap (convex <-> concave) above.
##
##   h is affine iff h is affine AND every g_j is affine.
##
## The composition is intentionally conservative: any expression we can't
## prove convex/concave/affine is reported as "unknown". The Phase 3
## minorization engine treats "unknown" as a hard error and points the
## user at the offending sub-expression.
##
## Extended composition rules (applied as a fallback after the standard
## rule returns "unknown"):
##
## (E3) log(exp(h)) = h  (algebraic identity):
##   Curvature and sign are those of h.  Checked before E1/E2 so the
##   exact identity takes priority over the looser curvature bounds.
##
## (E4) exp(log(h)) = h  (algebraic identity, requires h > 0):
##   Curvature and sign are those of h.  Checked before E1/E2.
##
## (E5) log(h^c) = c * log(h)  (logarithm power rule, requires h > 0):
##   The result is an affine scaling of log(h), so its curvature is
##   "concave" when c > 0, "convex" when c < 0, and "affine" when c = 0.
##   Checked before E1/E2 to give the exact curvature rather than the
##   looser bound E1 would produce.
##
## (E1) Concave+nondecreasing atom applied to a convex inner:
##   If g is concave and nondecreasing, and the inner h is a flat additive
##   sum of univariate/constant terms each with nonneg sign, then g(h) is
##   "concave".  Justification: Jensen's inequality applies per additive
##   slot (see try_additive_jensen in minorize.R).
##
## (E2) Convex+nondecreasing atom applied to a concave inner:
##   If g is convex and nondecreasing, and the inner h is concave, then
##   g(h) is "convex".  Justification: supporting hyperplane on g at
##   h(theta_0) plus recursive minorization of the concave h
##   (see try_hyperplane_concave_inner in minorize.R).

# ---- Sign lattice helpers --------------------------------------------------
# These helpers are also used inside atom dcp_info functions in atoms.R.
# Lazy binding makes the order of source-time loading irrelevant.

is_pos    <- function(s) identical(s, "positive")
is_nonneg <- function(s) s %in% c("positive", "nonneg", "zero")
is_neg    <- function(s) identical(s, "negative")
is_nonpos <- function(s) s %in% c("negative", "nonpos", "zero")

flip_sign <- function(s) {
  switch(s,
    "positive" = "negative",
    "negative" = "positive",
    "nonneg"   = "nonpos",
    "nonpos"   = "nonneg",
    "zero"     = "zero",
    "unknown")
}

# Sign of a + b under the lattice {positive, nonneg, zero, nonpos, negative,
# unknown}. Conservative: when in doubt, return "unknown".
add_signs <- function(a, b) {
  if (a == "zero") return(b)
  if (b == "zero") return(a)
  if (a == "positive" && (b == "positive" || b == "nonneg")) return("positive")
  if (b == "positive" &&  a == "nonneg")                     return("positive")
  if (a == "nonneg"   &&  b == "nonneg")                     return("nonneg")
  if (a == "negative" && (b == "negative" || b == "nonpos")) return("negative")
  if (b == "negative" &&  a == "nonpos")                     return("negative")
  if (a == "nonpos"   &&  b == "nonpos")                     return("nonpos")
  "unknown"
}

reduce_signs <- function(signs) {
  if (length(signs) == 0L) return("zero")
  Reduce(add_signs, signs)
}

# Sign of c * x where c is a numeric scalar.
scale_sign <- function(c_val, sx) {
  if (c_val == 0) return("zero")
  if (c_val > 0)  return(sx)
  flip_sign(sx)
}

# ---- Curvature composition rule -------------------------------------------

# Apply the DCP composition rule. Inputs:
#   outer_curv   - curvature of the outer atom: one of
#                  "affine"/"convex"/"concave"/"unknown"
#   mono_per_arg - character vector of length k: "nondecreasing",
#                  "nonincreasing", or "unknown"
#   inner_curv   - character vector of length k: composed curvatures of
#                  the children
#
# Returns one of "affine", "convex", "concave", "unknown".
compose_curvature <- function(outer_curv, mono_per_arg, inner_curv) {
  k <- length(inner_curv)
  if (k == 0L) return(outer_curv)

  # Test compatibility for the result being target ("convex" or "concave").
  is_compatible_for <- function(target) {
    same   <- target                                       # h nondecreasing
    flip   <- if (target == "convex") "concave" else "convex"  # h nonincreasing
    all(vapply(seq_len(k), function(j) {
      ic <- inner_curv[j]
      mo <- mono_per_arg[j]
      if (ic == "affine") return(TRUE)
      if (ic == same && mo == "nondecreasing") return(TRUE)
      if (ic == flip && mo == "nonincreasing") return(TRUE)
      FALSE
    }, logical(1)))
  }

  if (outer_curv == "affine" && all(inner_curv == "affine")) {
    return("affine")
  }

  if (outer_curv == "convex"  && is_compatible_for("convex"))  return("convex")
  if (outer_curv == "concave" && is_compatible_for("concave")) return("concave")

  if (outer_curv == "affine") {
    okv <- is_compatible_for("convex")
    okc <- is_compatible_for("concave")
    if (okv && okc) return("affine")
    if (okv)        return("convex")
    if (okc)        return("concave")
  }

  "unknown"
}

# ---- Extended composition rules --------------------------------------------

# Check whether every additive leaf of `expr` is either theta-independent
# (constant) or depends on exactly one theta index (univariate), and that
# every such leaf has a provably nonneg sign.  This is the structural
# precondition for the extended rule E1 (additive Jensen).
inner_is_nonneg_additive <- function(expr) {
  leaves <- collect_additive_leaves(expr)
  if (is.null(leaves)) return(FALSE)
  all(vapply(leaves, function(leaf) {
    s <- infer_dcp(leaf)$sign
    is_nonneg(s)
  }, logical(1)))
}

# Flatten an expression into its additive leaves (handling add/neg/scale).
# Returns NULL if any leaf is itself a multi-argument call other than add.
collect_additive_leaves <- function(expr) {
  if (inherits(expr, "mmad_call") && expr$op == "add") {
    out <- list()
    for (a in expr$args) {
      sub <- collect_additive_leaves(a)
      if (is.null(sub)) return(NULL)
      out <- c(out, sub)
    }
    return(out)
  }
  list(expr)
}

# Extended composition fallback. Called by infer_dcp when the standard
# compose_curvature returns "unknown" for a 1-argument mmad_call.
# Returns a list(curvature, sign) rather than just a curvature string,
# because E3/E4 pass through the inner's sign exactly.
extended_compose_curvature <- function(atom_op, atom_curv, atom_mono,
                                       inner_curv, inner_sign, inner_expr) {
  # E3: log(exp(h)) = h  -- exact algebraic identity.
  # Curvature and sign are exactly those of h, regardless of h's shape.
  if (atom_op == "log" &&
      inherits(inner_expr, "mmad_call") && inner_expr$op == "exp") {
    # The grandchild is the argument of exp, i.e. h in log(exp(h)).
    grandchild <- infer_dcp(inner_expr$args[[1L]])
    return(list(curvature = grandchild$curvature, sign = grandchild$sign))
  }

  # E4: exp(log(h)) = h  -- exact algebraic identity (requires h > 0).
  # We only apply this when the inner log's argument is provably positive,
  # which is the domain condition for log anyway.
  if (atom_op == "exp" &&
      inherits(inner_expr, "mmad_call") && inner_expr$op == "log") {
    grandchild <- infer_dcp(inner_expr$args[[1L]])
    # Only apply if the argument of log is provably positive (valid domain).
    if (is_pos(grandchild$sign) || is_nonneg(grandchild$sign)) {
      return(list(curvature = grandchild$curvature, sign = grandchild$sign))
    }
  }

  # E5: log(h^c) = c * log(h)  -- logarithm power rule.
  # Requires h > 0 (domain condition for log; also needed for h^c when
  # c is non-integer).  The result c*log(h) has curvature sign(c)*concave.
  if (atom_op == "log" &&
      inherits(inner_expr, "mmad_call") && inner_expr$op == "pow") {
    cc         <- inner_expr$params$c
    grandchild <- infer_dcp(inner_expr$args[[1L]])
    if (is_pos(grandchild$sign)) {
      # Curvature of c * log(h): log(h) has curvature "concave"; scaling
      # by c flips it when c < 0.
      log_h_curv <- infer_dcp(
        mmad_call("log", list(inner_expr$args[[1L]]))
      )$curvature
      result_curv <- if (cc == 0) {
        "affine"
      } else if (cc > 0) {
        log_h_curv          # same curvature as log(h)
      } else {
        # cc < 0: flip concave <-> convex
        if (log_h_curv == "concave") "convex"
        else if (log_h_curv == "convex") "concave"
        else log_h_curv     # affine or unknown pass through
      }
      return(list(curvature = result_curv, sign = "unknown"))
    }
  }

  # E1: concave+nondecreasing of a nonneg additive-convex inner => concave.
  if (atom_curv == "concave" && atom_mono == "nondecreasing" &&
      inner_curv == "convex"  && inner_is_nonneg_additive(inner_expr)) {
    return(list(curvature = "concave", sign = "unknown"))
  }

  # E2: convex+nondecreasing of a concave inner => convex.
  if (atom_curv == "convex" && atom_mono == "nondecreasing" &&
      inner_curv == "concave") {
    return(list(curvature = "convex", sign = "unknown"))
  }

  list(curvature = "unknown", sign = "unknown")
}

# ---- Public API ------------------------------------------------------------

#' Inferred curvature of an `mmad_expr`
#'
#' Returns the DCP curvature of the expression as one of `"affine"`,
#' `"convex"`, `"concave"`, or `"unknown"`. `"unknown"` means the DCP
#' rule set could not prove a definite curvature; this is a refusal,
#' not a claim of non-convexity.
#'
#' @param expr An `mmad_expr`.
#' @return A length-1 character.
#' @examples
#' curvature(log(mmad_var(1) + 1))                      # "concave"
#' curvature(exp(mmad_var(1)))                          # "convex"
#' curvature(2 * mmad_var(1) + 3 * mmad_var(2))         # "affine"
#' curvature(exp(mmad_var(1)) - exp(mmad_var(2)))       # "unknown"
#' @export
curvature <- function(expr) infer_dcp(expr)$curvature

#' Inferred sign of an `mmad_expr`
#'
#' @param expr An `mmad_expr`.
#' @return One of `"positive"`, `"nonneg"`, `"zero"`, `"nonpos"`,
#'   `"negative"`, `"unknown"`.
#' @examples
#' sign_of(exp(mmad_var(1)))                # "positive"
#' sign_of(mmad_var(1) ^ 2)                 # "nonneg"
#' @export
sign_of <- function(expr) infer_dcp(expr)$sign

#' Whether an expression's curvature is provably one of affine/convex/concave
#'
#' Convenience predicate: `TRUE` exactly when [curvature()] is not
#' `"unknown"`.
#'
#' @param expr An `mmad_expr`.
#' @return `TRUE` or `FALSE`.
#' @examples
#' is_dcp(log(mmad_var(1) + 1))                       # TRUE
#' is_dcp(exp(mmad_var(1)) - exp(mmad_var(2)))        # FALSE
#' @export
is_dcp <- function(expr) {
  curvature(expr) %in% c("affine", "convex", "concave")
}

# Internal recursion. Returns list(curvature, sign).
infer_dcp <- function(expr) {
  if (inherits(expr, "mmad_var")) {
    sg <- if (!is.null(expr$sign)) expr$sign else "unknown"
    return(list(curvature = "affine", sign = sg))
  }
  if (inherits(expr, "mmad_const")) {
    sg <- if (expr$value > 0) "positive"
          else if (expr$value < 0) "negative"
          else "zero"
    return(list(curvature = "affine", sign = sg))
  }
  if (inherits(expr, "mmad_call")) {
    children   <- lapply(expr$args, infer_dcp)
    inner_curv <- vapply(children, function(c) c$curvature, character(1))
    inner_sign <- vapply(children, function(c) c$sign,      character(1))

    atom <- mmad_atom(expr$op)
    if (is.null(atom$dcp_info)) {
      stop(sprintf("Atom '%s' has no dcp_info function registered.", expr$op))
    }
    info <- atom$dcp_info(inner_sign, expr$params)

    mono <- info$monotonicity
    # If a single monotonicity was provided for an n-ary atom, recycle.
    if (length(mono) == 1L && length(inner_curv) > 1L) {
      mono <- rep(mono, length(inner_curv))
    }

    curv_composed <- compose_curvature(info$curvature, mono, inner_curv)

    # If the standard rule cannot determine curvature and this is a
    # 1-argument atom, try the extended rules (E3, E4, E1, E2).
    if (curv_composed == "unknown" &&
        length(expr$args) == 1L &&
        length(mono) >= 1L) {
      ext <- extended_compose_curvature(
        atom_op    = expr$op,
        atom_curv  = info$curvature,
        atom_mono  = mono[1L],
        inner_curv = inner_curv[1L],
        inner_sign = inner_sign[1L],
        inner_expr = expr$args[[1L]])
      if (ext$curvature != "unknown") {
        return(list(curvature = ext$curvature, sign = ext$sign))
      }
    }

    return(list(curvature = curv_composed, sign = info$sign))
  }
  stop("infer_dcp(): unrecognised mmad_expr node of class ",
       paste(class(expr), collapse = "/"))
}
