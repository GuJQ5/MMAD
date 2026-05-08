------------------------------------------------------------------------
# MMAD

**Minorization-Maximization via Assembly-Decomposition Technology**

`MMAD` maximizes a target function via the minorization-maximization (MM) algorithm, using a formula-based interface, disciplined-convex inference of curvature, and a surrogate built from Jensen's inequality and the supporting hyperplane.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("GuJQ5/MMAD")
```

## Quick start

A small concave target:

``` r
library(MMAD)

fit <- mmad(~ log(theta[1] + theta[2]) - theta[1] - theta[2] + 2,
            init = c(2, 2))
fit$estimate    # ~ (0.5, 0.5)
```

The canonical multinomial-style objective from Tian et al. (2019):

``` r
fit <- mmad(~ 12 * log(0.5 * theta[1] + 0.5 * theta[2]) +
            15 * log((2/3) * theta[1] + (1/3) * theta[2]) +
             9 * log((1/3) * theta[1] + (2/3) * theta[2]) -
             6 * theta[1] - 6 * theta[2],
            init = c(4, 2))
summary(fit)
```

Poisson regression in one line:

``` r
set.seed(1)
n <- 100; p <- 3
X <- matrix(rnorm(n * p, sd = 0.5), nrow = n)
y <- rpois(n, exp(X %*% c(0.3, -0.2, 0.5)))

fit <- mmad(~ sum(y * (X %*% theta) - exp(X %*% theta)),
            init = rep(0, p),
            data = list(X = X, y = y))
fit$estimate    # close to glm(y ~ X - 1, family = poisson)$coef
```

## Pre-flight diagnostic

`Function_check()` reports the curvature, top-level summand breakdown, domain feasibility at the initial point, and whether the surrogate will be fully separable:

``` r
chk <- Function_check(~ log(theta[1] + theta[2]) - theta[1],
                      init = c(1, 1))
print(chk)
```

## Formula rules

The right-hand side of `~` is parsed into a symbolic expression tree; the layers below explain what the parser accepts.

### Parameters

Two equivalent ways to refer to the parameter vector:

- **Indexed**: `theta[1]`, `theta[2]`, ..., when `init` is unnamed (e.g. `init = c(4, 2)`).
- **Named**: bare symbols `alpha`, `beta`, ..., when `init` is named (e.g. `init = c(alpha = 4, beta = 2)`).

Both can be mixed in the same formula. Optional sign hints (`mmad_var(1, sign = "positive")`) sharpen DCP inference for atoms whose curvature depends on argument sign (notably `pow`).

### Operators

- `+`, `-`: any two operands, including two parameter-dependent expressions.
- `*`: one side must be a numeric scalar (or evaluate to one in `data`). `theta * theta` is rejected as non-DCP at parse time.
- `/`: the divisor must be a numeric scalar. `1 / theta[1]` is rejected.
- `^`: the exponent must be a numeric scalar. `theta[1] ^ theta[2]` is rejected.
- Unary `-` and `+` are supported.

### Math functions

`log`, `exp`, `sqrt`. Each applies to a single sub-expression. Domain (e.g. `log` requiring positive argument) is enforced at evaluation time, not at parse time -- the iterate's domain is the user's responsibility.

`log(1 - x)` requires no special syntax: the parser builds the natural additive structure, and the surrogate engine handles it via Jensen when the slot signs work out, or via the non-separable bucket otherwise.

### Theta-free sub-expressions

Any sub-expression that does not reference a parameter is evaluated eagerly in `data` (or the formula's environment) and frozen as a numeric constant. So `sum(X)`, `mean(Y)`, `crossprod(z)`, `2 * pi`, etc. all work as coefficients without quoting:

``` r
mmad(~ sum(X) * theta[1] - mean(X) * theta[2],
     init = c(0, 0), data = list(X = c(1, 2, 3, 4)))
```

### `sum()` over observations

When `sum()` wraps a parameter-dependent sub-expression that references vector data, the parser expands it row-by-row. For each row `i`:

- a vector `v` of length `n` in `data` is replaced with its `i`-th element `v[i]`;
- `X %*% theta` is replaced with the affine combination `X[i, 1] * theta[1] + ... + X[i, p] * theta[p]`.

The `n` row-versions are then assembled into a single sum. This makes log-likelihoods writable in one line, e.g. Poisson regression as `sum(y * (X %*% theta) - exp(X %*% theta))`.

### `X %*% theta`

The matrix-vector product is only valid **inside** `sum()`. Outside, the parser raises an error pointing the user at the `sum()` wrapper. Inside, the right-hand side must be the literal symbol `theta` (the full parameter vector, with `p = length(init)`); the left-hand side must resolve to a numeric matrix in `data` whose dimensions are `n x p`.

### What's not supported

- Products of two parameter-dependent expressions (`theta[1] * theta[2]`).
- Division by a parameter-dependent expression.
- Functions other than `log`, `exp`, `sqrt`. Custom atoms can be added via [`register_atom()`](%60?register_atom%60).
- Vector-valued returns outside `sum()`. `X %*% theta` standalone is non-scalar and currently rejected.
- Bare `theta` symbol outside `%*%` (use `theta[i]` for indexing).

## License

GPL-3
