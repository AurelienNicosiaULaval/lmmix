# Choosing lmmix and related R packages

## Purpose

Several mature R packages fit Gaussian mixed or repeated-measures
models, but their model and inference boundaries differ. This vignette
places `lmmix` within that ecosystem. It is a selection guide, not a
claim that one package is uniformly preferable.

The comparison is limited to documented capabilities that overlap with
the current `lmmix` scope. For definitive details, consult the cited
software papers and the documentation for the installed package
versions.

## Main distinctions

`lme4` is designed for mixed-effects models and uses sparse matrix
methods for random-effects computations ([Bates et al.
2015](#ref-bates2015)). `lmerTest` adds Satterthwaite and Kenward-Roger
fixed-effect inference to fitted `lmer` models ([Kuznetsova et al.
2017](#ref-kuznetsova2017)). This combination is a strong choice when
the covariance is driven by random effects and the residual errors are
conditionally independent.

`nlme` supports random effects together with within-group correlation
and variance functions ([Pinheiro and Bates 2000](#ref-pinheiro2000);
[Pinheiro et al. 2025](#ref-nlme2025)). It is a mature choice when those
model structures and its inference conventions meet the analysis
requirements.

`mmrm` targets marginal mixed models for repeated measures, with
structured within-subject covariance and small-sample inference, but
without Gaussian random effects in the fitted model ([Sabanes Bove et
al. 2026](#ref-mmrm2026)).

`lmmix` directly evaluates its profiled Gaussian ML or REML objective
for models that may combine random effects and one structured residual
covariance. Its post-processing layer supplies Satterthwaite or
Kenward-Roger inference, type III tests, estimated marginal means, and
pairwise contrasts from the same fitted object.

## Capability matrix

| Analysis requirement | `lmmix` | `lme4` and `lmerTest` | `nlme` | `mmrm` |
|:---|:---|:---|:---|:---|
| Gaussian random intercepts or slopes | Yes | Yes | Yes | No |
| Structured residual correlation | `id`, `cs`, `ar1`, `toep`, `toep(k)`, `un` | Conditionally independent residuals | Correlation structures documented by `nlme` | Repeated-measures covariance structures documented by `mmrm` |
| Random effects and correlated residuals in one model | Yes | No through the standard `lmer` interface | Yes | No |
| Satterthwaite inference | Yes | Yes through `lmerTest` | Not the package’s standard convention | Yes |
| Kenward-Roger inference | Yes for REML | Yes with supporting packages | Not the package’s standard convention | Yes |
| Type III tests and marginal means in the same package | Yes | Tests through `lmerTest`; means commonly through `emmeans` | Commonly assembled from multiple interfaces | Model inference plus `emmeans` interoperability |
| Large sparse crossed-effect focus | No | Yes | No general sparse crossed-effect claim | Subject-block computation |

The table summarizes interfaces rather than every possible extension. In
particular, `emmeans` interoperates with several of these model classes
and can provide a common marginal-means workflow.

## Choosing a package

Use `lme4` with `lmerTest` when the scientific covariance model is
adequately represented by random effects with conditionally independent
residual errors, especially for large sparse crossed designs.

Use `nlme` when you need its established combination of random effects,
correlation structures, and estimated variance functions, or when an
existing validated workflow depends on its inference conventions.

Use `mmrm` when the estimand and covariance specification are marginal,
the primary dependence is within subject, and no latent Gaussian random
effect is required.

Consider `lmmix` when the same model must contain Gaussian random
effects and a supported structured residual covariance, and the analysis
requires its integrated Satterthwaite or Kenward-Roger, type III, and
marginal-means layer. The [validation
vignette](https://aureliennicosiaulaval.github.io/lmmix/articles/validation.md)
documents the numerical comparisons that currently support this use.

## Computational boundary

`lmmix` 0.4.0 detects independent connected components induced jointly
by the random and repeated grouping structures. It factors those
covariance blocks separately during likelihood evaluation. This favors
subject-level and nested designs with many modest blocks.

The fitted object still retains a complete dense marginal covariance,
and some inference calculations remain dense. Crossed random effects can
connect all observations into one component. Therefore, `lmmix` does not
currently replace the sparse large-scale role of `lme4`, and no
universal maximum sample size is claimed. Benchmark the intended design,
covariance structure, and inference request before committing to a large
analysis.

## Current boundaries that affect selection

`lmmix` currently fits univariate Gaussian responses only. It rejects
rank-deficient fixed-effect matrices, estimates one repeated covariance
structure per model, and does not estimate residual variance strata
analogous to
[`nlme::varIdent()`](https://rdrr.io/pkg/nlme/man/varIdent.html). Known
relative precision weights are supported. Its prediction intervals
include fitted residual variance but not uncertainty in empirical random
effects. The maintained list of limitations and planned work is
available in the
[roadmap](https://aureliennicosiaulaval.github.io/lmmix/ROADMAP.md).

Bates, Douglas, Martin Mächler, Ben Bolker, and Steve Walker. 2015.
“Fitting Linear Mixed-Effects Models Using Lme4.” *Journal of
Statistical Software* 67 (1): 1–48.
<https://doi.org/10.18637/jss.v067.i01>.

Kuznetsova, Alexandra, Per B. Brockhoff, and Rune H. B. Christensen.
2017. “lmerTest Package: Tests in Linear Mixed Effects Models.” *Journal
of Statistical Software* 82 (13): 1–26.
<https://doi.org/10.18637/jss.v082.i13>.

Pinheiro, José C., and Douglas M. Bates. 2000. *Mixed-Effects Models in
s and s-PLUS*. Springer. <https://doi.org/10.1007/b98882>.

Pinheiro, José, Douglas Bates, and R Core Team. 2025. *Nlme: Linear and
Nonlinear Mixed Effects Models*.
<https://doi.org/10.32614/CRAN.package.nlme>.

Sabanes Bove, Daniel, Liming Li, Julia Dedic, et al. 2026. *Mmrm: Mixed
Models for Repeated Measures*.
<https://doi.org/10.32614/CRAN.package.mmrm>.
