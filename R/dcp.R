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
    return(list(curvature = curv_composed, sign = info$sign))
  }
  stop("infer_dcp(): unrecognised mmad_expr node of class ",
       paste(class(expr), collapse = "/"))
}
