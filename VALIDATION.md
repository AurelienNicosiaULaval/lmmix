# Validation against the source thesis

## Data provenance

The `multicentre` data reproduce Table 5.20, p. 144, of Mahsa Mohseni
Bonab's thesis *Programmation R et SAS pour modèles linéaires mixtes*. The
table was transcribed from the locally supplied final PDF and checked against
a rendered image of the source page.

The source PDF used for this transcription has SHA-256 checksum
`14d50797851080f64b7aa46eed230356e4085011c29ab1a9c50d24e152da8148`.

The resulting data contain 153 rows, including 28 missing responses. The model
therefore uses 125 observations, as reported in Table 5.6.

## Reference model

```r
fit <- lmm(
  multicentre,
  Y ~ Drug * Time,
  random = ~ 1 | Center,
  repeated = ~ Time | Center:Drug:Subject,
  structure = "ar1",
  method = "REML",
  ddf = "satterthwaite"
)
```

The active tests compare covariance parameters, fixed effects, type III tests,
LS-means, pairwise differences, and rowwise Satterthwaite degrees of freedom
with Tables 5.8 to 5.19.

## Annex B erratum

Annex B reports `16.25425` as the Drug 3 LS-mean. This is not the SAS LS-mean.
The data-aligned comparison in Table 5.17 gives `17.0516`, and the
automatically generated LS-mean is `17.05158`. Table 5.11 contains the same
set of estimates but presents the outer drug and time levels in reverse order.

The value `16.25425` comes from the prototype contrast on p. 109. That contrast
omits one third of the `Time3` main-effect coefficient. Its rows are equivalent
to:

```r
c(1, 0, 1, 1 / 3, 0, 0, 1 / 3, 0, 1 / 3)
```

The difference between the correct LS-mean and the prototype value is exactly
one third of the fitted `Time3` coefficient. A regression test preserves this
traceability without making the production `lsmeans()` function return an
incorrect marginal mean.

Some denominator degrees of freedom in Annex B are printed with only one or
two decimal places. Those values cannot support an absolute `1e-3` comparison.
The tests compare degrees of freedom at the precision printed in the source,
while estimates and standard errors are tested with an absolute tolerance of
`1e-3`.
