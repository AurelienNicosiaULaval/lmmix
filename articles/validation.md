# Validation against independent R implementations and PROC MIXED

## Purpose and evidentiary standard

This vignette reports the numerical comparisons used to validate
`lmmix`. Every comparison is restricted to a model class for which the
two implementations use the same fixed-effects specification, covariance
model, and estimation criterion. A small difference is meaningful only
after those definitions have been aligned.

Validation is organized in three layers:

1.  comparisons with independent R implementations for overlapping
    models;
2.  structural checks for combined random-effect and
    residual-correlation models that no single reference R function
    covers;
3.  regression comparisons with stored PROC MIXED targets for one fully
    specified combined model.

The numerical tables below are produced by executable R code when the
vignette is built. The PROC MIXED values are stored targets. The
repository does not contain a SAS execution transcript, so the SAS
section is not presented as a new independent SAS run.

## Comparison map

| lmmix model | reference function | quantities compared | tolerance |
|:---|:---|:---|:---|
| Random-intercept REML | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | fixed effects and log-likelihood: `1e-6`; variance components: `1e-4` | absolute |
| Random-slope REML | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | fixed effects and log-likelihood: `1e-5`; covariance parameters: `1e-4` | absolute |
| Random-intercept ML | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | fixed effects, log-likelihood | `1e-6` |
| Marginal CS and AR(1) REML | [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html) | fixed effects, residual variance, correlation, log-likelihood | `1e-5` |
| Marginal AR(1) REML | [`mmrm::mmrm()`](https://openpharma.github.io/mmrm/latest-tag/reference/mmrm.html) | fixed effects and log-likelihood: `1e-5`; residual covariance matrix: `1e-4` | absolute |
| Random-intercept REML | [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html) | estimates and standard errors: `1e-5`; coefficient df: `5e-3` | absolute |
| Type III Satterthwaite tests | [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html) | numerator df: `1e-6`; denominator df and F statistics: `5e-3` | absolute |
| Marginal AR(1) Kenward-Roger | [`mmrm::mmrm()`](https://openpharma.github.io/mmrm/latest-tag/reference/mmrm.html) | fixed effects: `1e-5`; standard errors: `3e-3`; denominator df: `5e-3` | absolute |
| Crossed random intercepts | [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) | log-likelihood: `1e-5`; covariance parameters: `1e-4` | absolute |
| Nested fixed-effect models | `lme4::anova()` | likelihood-ratio statistic and p-value: `1e-5` | absolute |
| Random center plus AR(1) residual | PROC MIXED stored targets | covariance parameters, fixed effects, type III tests, LS-means, differences | target-specific precision |

`nlme` provides independent likelihood comparisons for random-effects
and structured marginal models ([Pinheiro and Bates
2000](#ref-pinheiro2000); [Pinheiro et al. 2025](#ref-nlme2025)).
`lmerTest` provides an independent Satterthwaite implementation for
models fitted by `lme4` ([Kuznetsova et al. 2017](#ref-kuznetsova2017)).
`mmrm` provides a second independent comparison for the marginal
repeated-measures model ([Sabanes Bove et al. 2026](#ref-mmrm2026)).

## Shared longitudinal data

The R comparisons use the `Orthodont` data distributed with `nlme`.
`Occasion` is created explicitly so that AR(1) visits have a declared
order.

``` r

library(lmmix)

orthodont <- as.data.frame(nlme::Orthodont)
orthodont$Occasion <- factor(ave(
  seq_len(nrow(orthodont)),
  orthodont$Subject,
  FUN = seq_along
))

tibble::tibble(
  n.obs = nrow(orthodont),
  n.subjects = nlevels(orthodont$Subject),
  n.occasions = nlevels(orthodont$Occasion)
) |>
  show_table(digits = 0)
```

| n obs | n subjects | n occasions |
|------:|-----------:|------------:|
|   108 |         27 |           4 |

## Random-intercept REML comparison with `nlme::lme()`

The first comparison uses the same fixed formula, random intercept, and
REML criterion in both implementations.

``` r

fit.ri <- lmm(
  orthodont,
  distance ~ age + Sex,
  random = ~ 1 | Subject,
  method = "REML"
)

ref.ri <- nlme::lme(
  distance ~ age + Sex,
  random = ~ 1 | Subject,
  data = orthodont,
  method = "REML"
)

comparison_table(
  term = names(fixef(fit.ri)),
  lmmix = fixef(fit.ri),
  reference = nlme::fixef(ref.ri)
) |>
  show_table()
```

| term        |     lmmix | reference | abs diff |
|:------------|----------:|----------:|---------:|
| (Intercept) | 17.706713 | 17.706713 |        0 |
| age         |  0.660185 |  0.660185 |        0 |
| SexFemale   | -2.321023 | -2.321023 |        0 |

``` r

ref.ri.cov <- as.matrix(nlme::getVarCov(
  ref.ri,
  type = "random.effects"
))
ref.ri.variance <- c(
  ref.ri.cov[1L, 1L],
  ref.ri$sigma^2
)

comparison_table(
  term = c("random.intercept.var", "residual.var"),
  lmmix = VarCorr(fit.ri)$estimate,
  reference = ref.ri.variance
) |>
  show_table()
```

| term                 |    lmmix | reference | abs diff |
|:---------------------|---------:|----------:|---------:|
| random.intercept.var | 3.266790 |  3.266784 |  6.0e-06 |
| residual.var         | 2.049417 |  2.049456 |  3.9e-05 |

``` r


comparison_table(
  term = "logLik",
  lmmix = logLik(fit.ri),
  reference = logLik(ref.ri)
) |>
  show_table()
```

| term   |     lmmix | reference | abs diff |
|:-------|----------:|----------:|---------:|
| logLik | -218.7563 | -218.7563 |        0 |

## Random-slope REML comparison with `nlme::lme()`

This comparison adds a correlated random slope for age. It checks two
random variances, their correlation, the residual variance, the fixed
effects, and the REML log-likelihood.

``` r

fit.rs <- lmm(
  orthodont,
  distance ~ age + Sex,
  random = ~ 1 + age | Subject,
  method = "REML"
)

ref.rs <- nlme::lme(
  distance ~ age + Sex,
  random = ~ age | Subject,
  data = orthodont,
  method = "REML"
)

comparison_table(
  term = names(fixef(fit.rs)),
  lmmix = fixef(fit.rs),
  reference = nlme::fixef(ref.rs)
) |>
  show_table()
```

| term        |     lmmix | reference | abs diff |
|:------------|----------:|----------:|---------:|
| (Intercept) | 17.635201 | 17.635200 |    1e-06 |
| age         |  0.660185 |  0.660185 |    0e+00 |
| SexFemale   | -2.145493 | -2.145492 |    1e-06 |

``` r

ref.rs.cov <- as.matrix(nlme::getVarCov(
  ref.rs,
  type = "random.effects"
))
ref.rs.parameters <- c(
  ref.rs.cov[1L, 1L],
  ref.rs.cov[2L, 2L],
  ref.rs.cov[1L, 2L] / sqrt(ref.rs.cov[1L, 1L] * ref.rs.cov[2L, 2L]),
  ref.rs$sigma^2
)

comparison_table(
  term = c(
    "random.intercept.var",
    "random.age.var",
    "random.intercept.age.cor",
    "residual.var"
  ),
  lmmix = VarCorr(fit.rs)$estimate,
  reference = ref.rs.parameters
) |>
  show_table()
```

| term                     |     lmmix | reference | abs diff |
|:-------------------------|----------:|----------:|---------:|
| random.intercept.var     |  7.823339 |  7.823336 |    4e-06 |
| random.age.var           |  0.051269 |  0.051269 |    0e+00 |
| random.intercept.age.cor | -0.765844 | -0.765847 |    3e-06 |
| residual.var             |  1.716198 |  1.716204 |    7e-06 |

``` r


comparison_table(
  term = "logLik",
  lmmix = logLik(fit.rs),
  reference = logLik(ref.rs)
) |>
  show_table()
```

| term   |     lmmix | reference | abs diff |
|:-------|----------:|----------:|---------:|
| logLik | -217.6169 | -217.6169 |        0 |

## Maximum-likelihood comparison with `nlme::lme()`

The ML comparison verifies that agreement is not specific to the REML
criterion.

``` r

fit.ml <- lmm(
  orthodont,
  distance ~ age + Sex,
  random = ~ 1 | Subject,
  method = "ML"
)

ref.ml <- nlme::lme(
  distance ~ age + Sex,
  random = ~ 1 | Subject,
  data = orthodont,
  method = "ML"
)

rbind(
  comparison_table(
    term = names(fixef(fit.ml)),
    lmmix = fixef(fit.ml),
    reference = nlme::fixef(ref.ml)
  ),
  comparison_table(
    term = "logLik",
    lmmix = logLik(fit.ml),
    reference = logLik(ref.ml)
  )
) |>
  show_table()
```

| term        |       lmmix |   reference | abs diff |
|:------------|------------:|------------:|---------:|
| (Intercept) |   17.706713 |   17.706713 |        0 |
| age         |    0.660185 |    0.660185 |        0 |
| SexFemale   |   -2.321023 |   -2.321023 |        0 |
| logLik      | -217.428243 | -217.428243 |        0 |

## Marginal CS and AR(1) comparisons with `nlme::gls()`

These models contain no random effect. The repeated-measures covariance
is compound symmetric or AR(1), and both implementations optimize a REML
criterion.

``` r

marginal.fits <- list()
marginal.refs <- list()

for (structure in c("cs", "ar1")) {
  marginal.fits[[structure]] <- lmm(
    orthodont,
    distance ~ age + Sex,
    repeated = ~ Occasion | Subject,
    structure = structure,
    method = "REML"
  )

  correlation <- if (structure == "cs") {
    nlme::corCompSymm(form = ~ 1 | Subject)
  } else {
    nlme::corAR1(form = ~ 1 | Subject)
  }

  marginal.refs[[structure]] <- nlme::gls(
    distance ~ age + Sex,
    data = orthodont,
    correlation = correlation,
    method = "REML"
  )
}
```

``` r

marginal.fixed <- do.call(rbind, lapply(c("cs", "ar1"), function(structure) {
  fit <- marginal.fits[[structure]]
  ref <- marginal.refs[[structure]]
  table <- comparison_table(
    term = names(fixef(fit)),
    lmmix = fixef(fit),
    reference = stats::coef(ref)
  )
  tibble::add_column(table, structure = structure, .before = 1)
}))

show_table(marginal.fixed)
```

| structure | term        |     lmmix | reference | abs diff |
|:----------|:------------|----------:|----------:|---------:|
| cs        | (Intercept) | 17.706713 | 17.706713 |        0 |
| cs        | age         |  0.660185 |  0.660185 |        0 |
| cs        | SexFemale   | -2.321023 | -2.321023 |        0 |
| ar1       | (Intercept) | 17.878709 | 17.878709 |        0 |
| ar1       | age         |  0.652960 |  0.652960 |        0 |
| ar1       | SexFemale   | -2.418714 | -2.418714 |        0 |

``` r

marginal.covariance <- do.call(
  rbind,
  lapply(c("cs", "ar1"), function(structure) {
    fit <- marginal.fits[[structure]]
    ref <- marginal.refs[[structure]]
    reference.parameters <- c(
      ref$sigma^2,
      as.numeric(stats::coef(
        ref$modelStruct$corStruct,
        unconstrained = FALSE
      ))
    )
    table <- comparison_table(
      term = c("residual.var", "residual.cor"),
      lmmix = VarCorr(fit)$estimate,
      reference = reference.parameters
    )
    tibble::add_column(table, structure = structure, .before = 1)
  })
)

show_table(marginal.covariance)
```

| structure | term         |    lmmix | reference | abs diff |
|:----------|:-------------|---------:|----------:|---------:|
| cs        | residual.var | 5.316233 |  5.316240 |    7e-06 |
| cs        | residual.cor | 0.614490 |  0.614491 |    1e-06 |
| ar1       | residual.var | 5.296881 |  5.296881 |    0e+00 |
| ar1       | residual.cor | 0.625867 |  0.625867 |    0e+00 |

``` r


tibble::tibble(
  structure = c("cs", "ar1"),
  lmmix = vapply(
    marginal.fits,
    function(x) as.numeric(logLik(x)),
    numeric(1)
  ),
  reference = vapply(
    marginal.refs,
    function(x) as.numeric(logLik(x)),
    numeric(1)
  )
) |>
  transform(abs.diff = abs(lmmix - reference)) |>
  show_table()
```

| structure |     lmmix | reference | abs diff |
|:----------|----------:|----------:|---------:|
| cs        | -218.7563 | -218.7563 |        0 |
| ar1       | -222.7241 | -222.7241 |        0 |

## Marginal AR(1) comparison with `mmrm::mmrm()`

The `mmrm` comparison provides a second independent implementation of
the marginal AR(1) model. It checks the complete residual covariance
matrix rather than only the variance and correlation parameters.

``` r

ref.mmrm <- mmrm::mmrm(
  distance ~ age + Sex,
  data = orthodont,
  covariance = mmrm::cov_struct(
    type = "ar1",
    visits = "Occasion",
    subject = "Subject"
  ),
  reml = TRUE
)

fit.mmrm <- marginal.fits[["ar1"]]

comparison_table(
  term = names(fixef(fit.mmrm)),
  lmmix = fixef(fit.mmrm),
  reference = stats::coef(ref.mmrm)
) |>
  show_table()
```

| term        |     lmmix | reference | abs diff |
|:------------|----------:|----------:|---------:|
| (Intercept) | 17.878709 | 17.878710 |    1e-06 |
| age         |  0.652960 |  0.652960 |    0e+00 |
| SexFemale   | -2.418714 | -2.418715 |    1e-06 |

``` r


tibble::tibble(
  quantity = c("logLik", "residual.covariance.max.abs.diff"),
  lmmix = c(as.numeric(logLik(fit.mmrm)), NA_real_),
  reference = c(as.numeric(logLik(ref.mmrm)), NA_real_),
  abs.diff = c(
    abs(as.numeric(logLik(fit.mmrm)) - as.numeric(logLik(ref.mmrm))),
    max(abs(
      unname(fit.mmrm$covariance$residual_base) -
        unname(mmrm::component(ref.mmrm, "varcor"))
    ))
  )
) |>
  show_table()
```

| quantity                         |     lmmix | reference | abs diff |
|:---------------------------------|----------:|----------:|---------:|
| logLik                           | -222.7241 | -222.7241 |  0.0e+00 |
| residual.covariance.max.abs.diff |        NA |        NA |  4.3e-05 |

## Satterthwaite comparison with `lmerTest::lmer()`

The covariance model in this section is a random intercept with
independent residuals, which is supported by both packages. The
coefficient comparison checks estimates, standard errors, and
one-dimensional Satterthwaite degrees of freedom.

``` r

ref.lmer <- lmerTest::lmer(
  distance ~ age + Sex + (1 | Subject),
  data = orthodont,
  REML = TRUE
)

lmmix.coef <- generics::tidy(fit.ri)
lmer.coef <- summary(ref.lmer)$coefficients

tibble::tibble(
  term = lmmix.coef$term,
  lmmix.estimate = lmmix.coef$estimate,
  lmerTest.estimate = lmer.coef[, "Estimate"],
  estimate.abs.diff = abs(lmmix.estimate - lmerTest.estimate),
  lmmix.std.error = lmmix.coef$std.error,
  lmerTest.std.error = lmer.coef[, "Std. Error"],
  std.error.abs.diff = abs(lmmix.std.error - lmerTest.std.error),
  lmmix.df = lmmix.coef$df,
  lmerTest.df = lmer.coef[, "df"],
  df.abs.diff = abs(lmmix.df - lmerTest.df)
) |>
  show_table()
```

| term | lmmix estimate | lmerTest estimate | estimate abs diff | lmmix std error | lmerTest std error | std error abs diff | lmmix df | lmerTest df | df abs diff |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| (Intercept) | 17.706713 | 17.706713 | 0 | 0.833917 | 0.833922 | 5e-06 | 99.35300 | 99.35237 | 0.000639 |
| age | 0.660185 | 0.660185 | 0 | 0.061605 | 0.061606 | 1e-06 | 80.00151 | 80.00000 | 0.001514 |
| SexFemale | -2.321023 | -2.321023 | 0 | 0.761416 | 0.761417 | 0e+00 | 25.00002 | 25.00000 | 0.000021 |

The type III comparison uses the same fixed-effect terms and the same
Satterthwaite denominator-df method.

``` r

lmmix.type3 <- anova(fit.ri)
lmer.type3 <- anova(ref.lmer, type = 3)

tibble::tibble(
  term = lmmix.type3$term,
  lmmix.num.df = lmmix.type3$num.df,
  lmerTest.num.df = lmer.type3[, "NumDF"],
  num.df.abs.diff = abs(lmmix.num.df - lmerTest.num.df),
  lmmix.den.df = lmmix.type3$den.df,
  lmerTest.den.df = lmer.type3[, "DenDF"],
  den.df.abs.diff = abs(lmmix.den.df - lmerTest.den.df),
  lmmix.statistic = lmmix.type3$statistic,
  lmerTest.statistic = lmer.type3[, "F value"],
  statistic.abs.diff = abs(lmmix.statistic - lmerTest.statistic)
) |>
  show_table()
```

| term | lmmix num df | lmerTest num df | num df abs diff | lmmix den df | lmerTest den df | den df abs diff | lmmix statistic | lmerTest statistic | statistic abs diff |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| age | 1 | 1 | 0 | 80.00151 | 80 | 0.001514 | 114.840466 | 114.838287 | 0.002179 |
| Sex | 1 | 1 | 0 | 25.00002 | 25 | 0.000021 | 9.292108 | 9.292099 | 0.000009 |

## Kenward-Roger comparison with `mmrm::mmrm()`

The Kenward-Roger comparison uses a marginal AR(1) model, which both
packages fit by REML. It checks fixed effects, adjusted standard errors,
and denominator degrees of freedom.

``` r

fit.kr <- lmm(
  orthodont,
  distance ~ age + Sex,
  repeated = ~ Occasion | Subject,
  structure = "ar1",
  ddf = "kenward-roger"
)

mmrm.formula <- stats::as.formula(
  "distance ~ age + Sex + ar1(Occasion | Subject)",
  env = asNamespace("mmrm")
)
ref.kr <- mmrm::mmrm(
  mmrm.formula,
  data = orthodont,
  method = "Kenward-Roger"
)
fit.kr.table <- generics::tidy(fit.kr)
ref.kr.table <- coef(summary(ref.kr))

tibble::tibble(
  term = fit.kr.table$term,
  estimate.abs.diff = abs(fit.kr.table$estimate - coef(ref.kr)),
  std.error.abs.diff = abs(
    fit.kr.table$std.error - ref.kr.table[, "Std. Error"]
  ),
  df.abs.diff = abs(fit.kr.table$df - ref.kr.table[, "df"])
) |>
  show_table()
```

| term        | estimate abs diff | std error abs diff | df abs diff |
|:------------|------------------:|-------------------:|------------:|
| (Intercept) |             1e-06 |           0.001800 |    0.001019 |
| age         |             0e+00 |           0.000215 |    0.001574 |
| SexFemale   |             1e-06 |           0.001449 |    0.000018 |

## Crossed random effects and model comparison

The crossed-intercept comparison uses two independent grouping factors.
The likelihood-ratio comparison uses nested fixed-effect models and the
same ML covariance specification in both packages.

``` r

set.seed(42)
crossed <- expand.grid(
  site = factor(seq_len(8)),
  observer = factor(seq_len(6)),
  replicate = seq_len(3)
)
site.effect <- stats::rnorm(8, sd = 0.8)
observer.effect <- stats::rnorm(6, sd = 0.5)
crossed$response <- 2 +
  site.effect[crossed$site] +
  observer.effect[crossed$observer] +
  stats::rnorm(nrow(crossed), sd = 0.4)

fit.crossed <- lmm(
  crossed,
  response ~ 1,
  random = list(site = ~ 1 | site, observer = ~ 1 | observer)
)
ref.crossed <- lme4::lmer(
  response ~ 1 + (1 | site) + (1 | observer),
  data = crossed,
  REML = TRUE
)

tibble::tibble(
  quantity = c(
    "logLik", "site variance", "observer variance", "residual variance"
  ),
  lmmix = c(as.numeric(logLik(fit.crossed)), VarCorr(fit.crossed)$estimate),
  reference = c(
    as.numeric(logLik(ref.crossed)),
    as.data.frame(lme4::VarCorr(ref.crossed))$vcov
  ),
  abs.diff = abs(lmmix - reference)
) |>
  show_table()
```

| quantity          |      lmmix |  reference | abs diff |
|:------------------|-----------:|-----------:|---------:|
| logLik            | -96.991425 | -96.991425 |    0e+00 |
| site variance     |   0.306957 |   0.306963 |    6e-06 |
| observer variance |   0.539502 |   0.539504 |    2e-06 |
| residual variance |   0.157861 |   0.157861 |    0e+00 |

``` r


fit.reduced <- lmm(
  orthodont,
  distance ~ age,
  random = ~ 1 | Subject
)
fit.full <- lmm(
  orthodont,
  distance ~ age + Sex,
  random = ~ 1 | Subject
)
lmmix.lrt <- suppressMessages(anova(fit.reduced, fit.full))
ref.reduced <- lme4::lmer(
  distance ~ age + (1 | Subject),
  data = orthodont,
  REML = FALSE
)
ref.full <- lme4::lmer(
  distance ~ age + Sex + (1 | Subject),
  data = orthodont,
  REML = FALSE
)
reference.lrt <- anova(ref.reduced, ref.full, refit = FALSE)

tibble::tibble(
  quantity = c("likelihood-ratio statistic", "p-value"),
  lmmix = c(lmmix.lrt$Chisq[2], lmmix.lrt$p.value[2]),
  reference = c(reference.lrt$Chisq[2], reference.lrt$`Pr(>Chisq)`[2]),
  abs.diff = abs(lmmix - reference)
) |>
  show_table()
```

| quantity                   |    lmmix | reference | abs diff |
|:---------------------------|---------:|----------:|---------:|
| likelihood-ratio statistic | 8.533057 |  8.533057 |        0 |
| p-value                    | 0.003488 |  0.003488 |        0 |

## Combined covariance structures

The included `multicentre` data reproduce the multilocation
repeated-measures example in Section 28.3 of Milliken and Johnson
([2009](#ref-milliken2009)). They contain 153 rows, including 28 missing
responses. The complete-case models therefore use 125 observations.

No single external R function used above fits all of the following
components in one object: a random center effect, a structured
within-subject residual covariance, and the `lmmix` Satterthwaite
inference layer. For this combined case, all five residual structures
are checked for optimizer convergence, a positive-definite likelihood
Hessian, and a positive-definite marginal covariance matrix.

``` r

combined.fits <- lapply(
  c("id", "cs", "ar1", "toep", "un"),
  function(structure) {
    lmm(
      multicentre,
      Y ~ Drug * Time,
      random = ~ 1 | Center,
      repeated = ~ Time | Center:Drug:Subject,
      structure = structure,
      method = "REML",
      ddf = "satterthwaite"
    )
  }
)
names(combined.fits) <- c("id", "cs", "ar1", "toep", "un")

tibble::tibble(
  structure = names(combined.fits),
  convergence = vapply(
    combined.fits,
    function(x) x$convergence$code,
    integer(1)
  ),
  hessian.positive.definite = vapply(
    combined.fits,
    function(x) x$convergence$hessian_positive_definite,
    logical(1)
  ),
  min.marginal.eigenvalue = vapply(
    combined.fits,
    function(x) {
      min(eigen(
        x$covariance$v,
        symmetric = TRUE,
        only.values = TRUE
      )$values)
    },
    numeric(1)
  )
) |>
  show_table()
```

| structure | convergence | hessian positive definite | min marginal eigenvalue |
|:----------|------------:|:--------------------------|------------------------:|
| id        |           0 | TRUE                      |                9.458960 |
| cs        |           0 | TRUE                      |                0.618274 |
| ar1       |           0 | TRUE                      |                0.473359 |
| toep      |           0 | TRUE                      |                0.467054 |
| un        |           0 | TRUE                      |                0.383388 |

## Stored PROC MIXED comparison

The stored targets correspond to the following specification. The RANDOM
and REPEATED statements represent distinct covariance sources, and
`DDFM=SATTERTHWAITE` requests approximate denominator degrees of freedom
as documented by SAS Institute Inc. ([2015](#ref-sas2015)).

``` sas
proc mixed data=multicentre method=reml;
  class Center Drug Subject Time;
  model Y = Drug Time Drug*Time / ddfm=satterth solution;
  random Center;
  repeated Time / type=ar(1) subject=Subject(Center*Drug);
  lsmeans Drug Time / diff;
run;
```

The exact `lmmix` fit used in the comparison is:

``` r

sas.fit <- combined.fits[["ar1"]]
sas.fit
#> Linear mixed model fit by REML
#> Formula: `Y ~ Drug * Time`
#> Random (Center): `~1 | Center`
#> Repeated: `~Time | Center:Drug:Subject` (AR1)
#> Log-likelihood: -245.313
#> Convergence code: 0
```

### Covariance parameters and likelihood criterion

``` r

sas.covariance.target <- c(5.1737, 10.6702, 0.9351)
lmmix.covariance <- VarCorr(sas.fit)

tibble::tibble(
  parameter = paste(
    lmmix.covariance$group,
    lmmix.covariance$term,
    lmmix.covariance$component,
    sep = "."
  ),
  lmmix = lmmix.covariance$estimate,
  SAS.target = sas.covariance.target,
  abs.diff = abs(lmmix - SAS.target)
) |>
  show_table()
```

| parameter              |     lmmix | SAS target | abs diff |
|:-----------------------|----------:|-----------:|---------:|
| Center.(Intercept).var |  5.173828 |     5.1737 | 0.000128 |
| Residual.ar1.var       | 10.670170 |    10.6702 | 0.000030 |
| Residual.ar1.cor       |  0.935120 |     0.9351 | 0.000020 |

``` r


tibble::tibble(
  quantity = "deviance",
  lmmix = deviance(sas.fit),
  SAS.target = 490.62570532,
  abs.diff = abs(lmmix - SAS.target)
) |>
  show_table()
```

| quantity |    lmmix | SAS target | abs diff |
|:---------|---------:|-----------:|---------:|
| deviance | 490.6257 |   490.6257 |        0 |

### Fixed effects

``` r

sas.fixed.target <- c(
  12.0524590,
  4.7647059,
  3.9411765,
  0.9210892,
  2.3922237,
  -0.1427475,
  -0.4735052,
  0.8233152,
  0.3342713
)

sas.fixed.table <- comparison_table(
  term = names(fixef(sas.fit)),
  lmmix = fixef(sas.fit),
  reference = sas.fixed.target
)
names(sas.fixed.table) <- c("term", "lmmix", "SAS.target", "abs.diff")
show_table(sas.fixed.table)
```

| term        |     lmmix | SAS target | abs diff |
|:------------|----------:|-----------:|---------:|
| (Intercept) | 12.052381 |  12.052459 |  7.8e-05 |
| Drug2       |  4.764706 |   4.764706 |  0.0e+00 |
| Drug3       |  3.941176 |   3.941176 |  0.0e+00 |
| Time2       |  0.921092 |   0.921089 |  3.0e-06 |
| Time3       |  2.392232 |   2.392224 |  8.0e-06 |
| Drug2:Time2 | -0.142749 |  -0.142747 |  2.0e-06 |
| Drug3:Time2 | -0.473514 |  -0.473505 |  9.0e-06 |
| Drug2:Time3 |  0.823316 |   0.823315 |  1.0e-06 |
| Drug3:Time3 |  0.334252 |   0.334271 |  2.0e-05 |

### Type III tests

The stored type III targets contain F statistics rounded to two decimal
places and the interaction p-value. The other two SAS p-values were not
available in the stored target set and are therefore reported as
`not.stored` rather than reconstructed.

``` r

lmmix.sas.type3 <- anova(sas.fit)

tibble::tibble(
  term = lmmix.sas.type3$term,
  lmmix.num.df = lmmix.sas.type3$num.df,
  lmmix.den.df = lmmix.sas.type3$den.df,
  lmmix.statistic = lmmix.sas.type3$statistic,
  SAS.statistic = c(11.43, 59.27, 1.35),
  statistic.abs.diff = abs(lmmix.statistic - SAS.statistic),
  lmmix.p.value = lmmix.sas.type3$p.value,
  SAS.p.value = c("not.stored", "not.stored", "0.2597")
) |>
  show_table()
```

| term | lmmix num df | lmmix den df | lmmix statistic | SAS statistic | statistic abs diff | lmmix p value | SAS p value |
|:---|---:|---:|---:|---:|---:|---:|:---|
| Drug | 2 | 46.35956 | 11.428034 | 11.43 | 0.001966 | 0.000092 | not.stored |
| Time | 2 | 68.45549 | 59.267365 | 59.27 | 0.002635 | 0.000000 | not.stored |
| Drug:Time | 4 | 68.39047 | 1.352082 | 1.35 | 0.002082 | 0.259736 | 0.2597 |

### LS-means

The following table compares every stored Drug and Time marginal mean,
standard error, and denominator df.

``` r

drug.means <- lsmeans(sas.fit, ~Drug)
time.means <- lsmeans(sas.fit, ~Time)

lmmix.means <- rbind(
  tibble::tibble(
    factor = "Drug",
    level = as.character(drug.means$Drug),
    estimate = drug.means$estimate,
    std.error = drug.means$std.error,
    df = drug.means$df
  ),
  tibble::tibble(
    factor = "Time",
    level = as.character(time.means$Time),
    estimate = time.means$estimate,
    std.error = time.means$std.error,
    df = time.means$df
  )
)

sas.mean.estimate <- c(
  13.1568, 18.1484, 17.0516,
  14.9543, 15.6700, 17.7324
)
sas.mean.std.error <- c(
  1.5290, 1.5324, 1.5327,
  1.3961, 1.3985, 1.4032
)
sas.mean.df <- c(2.98, 3.01, 3.01, 2.08, 2.09, 2.12)

tibble::tibble(
  factor = lmmix.means$factor,
  level = lmmix.means$level,
  lmmix.estimate = lmmix.means$estimate,
  SAS.estimate = sas.mean.estimate,
  estimate.abs.diff = abs(lmmix.estimate - SAS.estimate),
  lmmix.std.error = lmmix.means$std.error,
  SAS.std.error = sas.mean.std.error,
  std.error.abs.diff = abs(lmmix.std.error - SAS.std.error),
  lmmix.df = lmmix.means$df,
  SAS.df = sas.mean.df
) |>
  show_table()
```

| factor | level | lmmix estimate | SAS estimate | estimate abs diff | lmmix std error | SAS std error | std error abs diff | lmmix df | SAS df |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| Drug | 1 | 13.15682 | 13.1568 | 2.3e-05 | 1.529053 | 1.5290 | 5.3e-05 | 2.979813 | 2.98 |
| Drug | 2 | 18.14838 | 18.1484 | 1.6e-05 | 1.532408 | 1.5324 | 8.0e-06 | 3.006723 | 3.01 |
| Drug | 3 | 17.05158 | 17.0516 | 2.2e-05 | 1.532758 | 1.5327 | 5.8e-05 | 3.008745 | 3.01 |
| Time | 1 | 14.95434 | 14.9543 | 4.2e-05 | 1.396163 | 1.3961 | 6.3e-05 | 2.077523 | 2.08 |
| Time | 2 | 15.67001 | 15.6700 | 1.3e-05 | 1.398482 | 1.3985 | 1.8e-05 | 2.091459 | 2.09 |
| Time | 3 | 17.73243 | 17.7324 | 3.0e-05 | 1.403262 | 1.4032 | 6.2e-05 | 2.120423 | 2.12 |

### Pairwise LS-mean differences

``` r

drug.pairs <- lsmeans(sas.fit, pairwise ~ Drug)$contrasts
time.pairs <- lsmeans(sas.fit, pairwise ~ Time)$contrasts

lmmix.pairs <- rbind(
  tibble::add_column(drug.pairs, factor = "Drug", .before = 1),
  tibble::add_column(time.pairs, factor = "Time", .before = 1)
)

sas.pair.estimate <- c(
  -4.9916, -3.8948, 1.0968,
  -0.7157, -2.7781, -2.0624
)
sas.pair.std.error <- c(
  1.0990, 1.0987, 1.1043,
  0.1845, 0.2714, 0.2056
)
sas.pair.df <- c(46.1, 46.1, 46.9, 68.1, 73.2, 68.6)
sas.pair.p.value <- c(
  "0.00004", "0.0009", "0.3257",
  "0.0002", "<0.0001", "<0.0001"
)

tibble::tibble(
  factor = lmmix.pairs$factor,
  contrast = lmmix.pairs$contrast,
  lmmix.estimate = lmmix.pairs$estimate,
  SAS.estimate = sas.pair.estimate,
  estimate.abs.diff = abs(lmmix.estimate - SAS.estimate),
  lmmix.std.error = lmmix.pairs$std.error,
  SAS.std.error = sas.pair.std.error,
  std.error.abs.diff = abs(lmmix.std.error - SAS.std.error),
  lmmix.df = lmmix.pairs$df,
  SAS.df = sas.pair.df,
  lmmix.p.value = lmmix.pairs$p.value,
  SAS.p.value = sas.pair.p.value
) |>
  show_table()
```

| factor | contrast | lmmix estimate | SAS estimate | estimate abs diff | lmmix std error | SAS std error | std error abs diff | lmmix df | SAS df | lmmix p value | SAS p value |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| Drug | 1 - 2 | -4.991561 | -4.9916 | 3.9e-05 | 1.098965 | 1.0990 | 3.5e-05 | 46.10872 | 46.1 | 0.000040 | 0.00004 |
| Drug | 1 - 3 | -3.894756 | -3.8948 | 4.4e-05 | 1.098726 | 1.0987 | 2.6e-05 | 46.05487 | 46.1 | 0.000914 | 0.0009 |
| Drug | 2 - 3 | 1.096806 | 1.0968 | 6.0e-06 | 1.104277 | 1.1043 | 2.3e-05 | 46.92852 | 46.9 | 0.325690 | 0.3257 |
| Time | 1 - 2 | -0.715671 | -0.7157 | 2.9e-05 | 0.184507 | 0.1845 | 7.0e-06 | 68.11869 | 68.1 | 0.000239 | 0.0002 |
| Time | 1 - 3 | -2.778088 | -2.7781 | 1.2e-05 | 0.271413 | 0.2714 | 1.3e-05 | 73.23889 | 73.2 | 0.000000 | \<0.0001 |
| Time | 2 - 3 | -2.062417 | -2.0624 | 1.7e-05 | 0.205595 | 0.2056 | 5.0e-06 | 68.60576 | 68.6 | 0.000000 | \<0.0001 |

Estimates, covariance parameters, and standard errors are tested with
absolute tolerance `1e-3`. Degrees of freedom are compared at the
precision stored in the targets. Type III F statistics are compared
after rounding to two decimal places. Inequalities such as `<0.0001` are
checked as inequalities, not as invented exact p-values.

## Automated regression tests

The numerical claims in this vignette are enforced in the package test
suite.

| test file | enforced comparison |
|:---|:---|
| `test-likelihood.R` | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html), [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html), and [`mmrm::mmrm()`](https://openpharma.github.io/mmrm/latest-tag/reference/mmrm.html) likelihood and covariance comparisons |
| `test-inference.R` | [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html) coefficient df and type III Satterthwaite comparisons |
| `test-validation-sas.R` | every stored PROC MIXED covariance, fixed-effect, LS-mean, and difference target |
| `test-covariance.R` | positive definiteness and structure-specific covariance checks |
| `test-methods.R` | S3 outputs, `ggplot2` diagnostics, and readable printed table headings |

## Reproducibility information

The versions used to build this vignette are reported rather than
assumed.

``` r

validation.packages <- c("lmmix", "nlme", "lmerTest", "mmrm")

tibble::tibble(
  package = validation.packages,
  version = vapply(
    validation.packages,
    function(package) {
      if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        "not.installed"
      }
    },
    character(1)
  )
) |>
  show_table()
```

| package  | version |
|:---------|:--------|
| lmmix    | 0.2.0   |
| nlme     | 3.1.169 |
| lmerTest | 3.2.1   |
| mmrm     | 0.3.18  |

``` r


tibble::tibble(
  R.version.string = R.version.string,
  platform = R.version$platform
) |>
  show_table()
```

| R version string             | platform            |
|:-----------------------------|:--------------------|
| R version 4.6.1 (2026-06-24) | x86_64-pc-linux-gnu |

## Boundaries of the evidence

The validation covers Kenward-Roger inference for overlapping
random-effects and marginal models, crossed random intercepts, and
fixed-effect likelihood-ratio comparisons. It does not validate
generalized responses or large-scale sparse performance because these
are outside the model scope. For the combined random-effect and
correlated-residual case, Kenward-Roger has structural tests but no
fresh independent SAS execution. An archived SAS log and output would be
required to upgrade that evidence.

Kuznetsova, Alexandra, Per B. Brockhoff, and Rune H. B. Christensen.
2017. “lmerTest Package: Tests in Linear Mixed Effects Models.” *Journal
of Statistical Software* 82 (13): 1–26.
<https://doi.org/10.18637/jss.v082.i13>.

Milliken, George A., and Dallas E. Johnson. 2009. *Analysis of Messy
Data, Volume 1: Designed Experiments*. 2nd ed. Chapman; Hall/CRC.
<https://doi.org/10.1201/EBK1584883340>.

Pinheiro, José C., and Douglas M. Bates. 2000. *Mixed-Effects Models in
s and s-PLUS*. Springer. <https://doi.org/10.1007/b98882>.

Pinheiro, José, Douglas Bates, and R Core Team. 2025. *Nlme: Linear and
Nonlinear Mixed Effects Models*.
<https://doi.org/10.32614/CRAN.package.nlme>.

Sabanes Bove, Daniel, Liming Li, Julia Dedic, et al. 2026. *Mmrm: Mixed
Models for Repeated Measures*.
<https://doi.org/10.32614/CRAN.package.mmrm>.

SAS Institute Inc. 2015. *SAS/STAT 14.1 User’s Guide: The MIXED
Procedure*. Cary, NC.
<https://support.sas.com/documentation/cdl/en/statug/68162/HTML/default/statug_mixed_syntax.htm>.
