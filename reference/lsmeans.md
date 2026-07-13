# Estimated marginal means for an `lmm` model

Marginal means use equal weights over nuisance-factor levels. Numeric
covariates not listed in `specs` are held at their observed means.

## Usage

``` r
lsmeans(
  object,
  specs,
  pairwise = FALSE,
  at = list(),
  level = 0.95,
  adjust = "none",
  ...
)
```

## Arguments

- object:

  An `lmm` object.

- specs:

  Variables defining the marginal means, supplied as a character vector
  or one-sided formula. Use `pairwise ~ factor` to request differences.

- pairwise:

  Whether to return pairwise differences as well as means.

- at:

  Named list overriding reference-grid values.

- level:

  Confidence level.

- adjust:

  Multiplicity adjustment passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- ...:

  Reserved for future extensions.

## Value

A tibble, or a list with `lsmeans` and `contrasts` tibbles when pairwise
comparisons are requested.
