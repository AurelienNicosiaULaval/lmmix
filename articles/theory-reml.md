# Estimation and inference in lmmix

## Gaussian mixed-model representation

For a response vector (Y), `lmmix` uses

``` math
Y = X\beta + Zu + \varepsilon,
\qquad
u \sim N(0, G),
\qquad
\varepsilon \sim N(0, R),
\qquad
u \perp \varepsilon.
```

The marginal model is

``` math
Y \sim N(X\beta, V),
\qquad
V(\eta) = ZG(\eta)Z^\mathsf{T} + R(\eta),
```

where (eta) is an unconstrained numerical parameter vector. This is the
standard Gaussian mixed-model decomposition described by Harville
([1977](#ref-harville1977)) and Pinheiro and Bates
([2000](#ref-pinheiro2000)).

## Covariance parameterization

The optimizer works on an unconstrained scale. The package maps those
values to admissible covariance matrices as follows.

| Component | Parameterization |
|:---|:---|
| Random-effect covariance `G` | Lower-triangular Cholesky factor with log diagonal |
| Independent residuals `id` | One variance on the log scale |
| Compound symmetry `cs` | Log variance and logistic correlation constrained to $`(-1/(q-1), 1)`$ |
| Autoregressive `ar1` | Log variance and `tanh` correlation in $`(-1, 1)`$ |
| Toeplitz `toep` | Log variance and partial autocorrelations mapped to a stationary Toeplitz correlation |
| Banded Toeplitz `toep(k)` | First `k` covariance bands estimated; longer-lag covariances fixed to zero |
| Unstructured `un` | Lower-triangular Cholesky factor with log diagonal |

Here (q) is the number of distinct repeated-measurement occasions.
Cholesky parameterizations guarantee positive definiteness for `G` and
`un`. The partial-autocorrelation construction guarantees a
positive-definite full Toeplitz correlation matrix. For `toep(k)`,
candidate banded blocks are checked by Cholesky factorization during
likelihood evaluation, and non-positive- definite candidates are
rejected.

## Profiled ML and REML criteria

For a fixed covariance matrix (V), the generalized least-squares
estimator is

``` math
\widehat\beta(\eta) =
(X^\mathsf{T}V^{-1}X)^{-1}X^\mathsf{T}V^{-1}Y.
```

Let (r = Y-X). The profiled negative ML log-likelihood is

``` math
\ell_{\mathrm{ML}}^-(\eta)
= \frac{1}{2}\left[
n\log(2\pi) + \log|V| + r^\mathsf{T}V^{-1}r
\right].
```

The profiled negative REML log-likelihood is

``` math
\ell_{\mathrm{REML}}^-(\eta)
= \frac{1}{2}\left[
(n-p)\log(2\pi)
+ \log|V|
+ \log|X^\mathsf{T}V^{-1}X|
+ r^\mathsf{T}V^{-1}r
\right].
```

REML originates with Patterson and Thompson
([1971](#ref-patterson1971)). General derivations are given by Harville
([1977](#ref-harville1977)) and LaMotte ([2007](#ref-lamotte2007)).

The implementation uses Cholesky solves for (V^{-1}X), (V^{-1}Y), log
determinants and the GLS normal equations. It does not form an explicit
matrix inverse in the likelihood calculation.

## Optimization and covariance-parameter uncertainty

[`lmm_control()`](https://aureliennicosiaulaval.github.io/lmmix/reference/lmm_control.md)
can select [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) or
[`optim()`](https://rdrr.io/r/stats/optim.html) directly. Its default
`optimizer = "auto"` starts with
[`nlminb()`](https://rdrr.io/r/stats/nlminb.html) and uses deterministic
perturbed starts and alternative optimizers only when an earlier attempt
does not converge. Every attempt is retained in
`fit$convergence$attempts`. After optimization, `numDeriv` computes the
observed likelihood Hessian on the unconstrained (eta) scale.

The inverse Hessian estimates

``` math
\widehat\Sigma_\eta =
\operatorname{Var}(\widehat\eta).
```

The package checks the Hessian eigenvalues. A non-positive-definite
Hessian triggers a warning because covariance-parameter standard errors
and Satterthwaite calculations may then be unavailable. Natural-scale
covariance standard errors are obtained by applying the delta method to
the transformation from (eta) to variances and correlations.

## Satterthwaite denominator degrees of freedom

For a one-row contrast (L), define

``` math
q(\eta) = L\operatorname{Var}(\widehat\beta)L^\mathsf{T}.
```

If

``` math
g = \frac{\partial q(\eta)}{\partial\eta},
```

the implemented Satterthwaite approximation is

``` math
\nu =
\frac{2q(\widehat\eta)^2}
{g^\mathsf{T}\widehat\Sigma_\eta g}.
```

This follows Satterthwaite ([1946](#ref-satterthwaite1946)). The
gradient of (q) is evaluated numerically on the same unconstrained
covariance scale used by the optimizer.

With `ddf = "residual"`, the package instead uses (n-p).

## Kenward-Roger adjustment

For `ddf = "kenward-roger"`, the package applies the small-sample
adjustment of Kenward and Roger ([1997](#ref-kenward1997)). This method
is available only with REML. Let

``` math
\Phi = (X^\mathsf{T}V^{-1}X)^{-1}
```

and let (W) denote the covariance matrix of the estimated unconstrained
covariance parameters. First and second numerical derivatives of (V) are
used to construct the (P_h), (Q\_{hj}), and (R\_{hj}) matrices in the
Kenward-Roger expansion. The adjusted fixed-effect covariance is

``` math
\widehat\Phi_A = \widehat\Phi +
2\widehat\Phi\left\{\sum_{h,j}W_{hj}
\left(Q_{hj} - P_h\widehat\Phi P_j - \frac{1}{4}R_{hj}\right)
\right\}\widehat\Phi.
```

Coefficient tests use this adjusted covariance. Multi-degree-of-freedom
tests also use the Kenward-Roger scale factor and denominator degrees of
freedom. `vcov(fit)` returns the adjusted covariance, while
`vcov(fit, adjusted = FALSE)` returns the unadjusted GLS covariance.

## Type III tests

Type III hypotheses are constructed with sum contrasts, then transformed
back to the fitted coefficient parameterization. For a multi-row
hypothesis matrix (L), the covariance of (L) is eigendecomposed into
orthogonal one-dimensional components. Their Satterthwaite degrees of
freedom are combined with the Fai-Cornelius moment-matching formula
([Fai and Cornelius 1996](#ref-fai1996)).

For one model, [`anova()`](https://rdrr.io/r/stats/anova.html) returns
type III tests. With two or more nested models, it performs
likelihood-ratio comparisons. Models that differ in fixed effects are
compared under ML; REML fits are refitted automatically by default.
Covariance-model comparisons must already use ML. Their chi-squared
reference distribution can be approximate when the null value lies on
the boundary of the parameter space ([Self and Liang
1987](#ref-self1987)). Setting `test = "parametric.bootstrap"` simulates
the reference distribution under each smaller fitted model, refits both
models by ML, and uses the finite- simulation corrected tail proportion
([Davison and Hinkley 1997](#ref-davison1997)).

Covariance-parameter confidence intervals are available with
`confint(fit, parm = "theta_")`. Variance intervals use a log-scale
delta method and correlation intervals use a Fisher-z delta method, so
their limits respect the natural parameter bounds. They remain local
Wald intervals based on the observed Hessian, not profile-likelihood
intervals.

## Estimated marginal means

[`lsmeans()`](https://aureliennicosiaulaval.github.io/lmmix/reference/lsmeans.md)
constructs a fixed-effects reference grid from the model terms.
Requested factors vary over their levels, nuisance factors receive equal
weights, and unrequested numeric covariates are held at their observed
means. The corresponding rows of the fixed-effects model matrix are
averaged to form the contrast matrix. This targets population marginal
means as described by Searle et al. ([1980](#ref-searle1980)).

Pairwise contrasts are differences between automatically generated
marginal mean rows. Multiplicity correction is applied to pairwise
p-values with
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). When
p-values are adjusted, `conf_adjust = "auto"` produces simultaneous
Bonferroni confidence intervals. Pointwise intervals remain available
with `conf_adjust = "none"`.

## Empirical BLUPs and prediction

For each independent random-effect term (k), covariance estimation is
followed by the empirical BLUP

``` math
\widehat u_k = G_kZ_k^\mathsf{T}V^{-1}(Y-X\widehat\beta).
```

Conditional fitted values use (X+ Zu); marginal fitted values use (X).
For new data, known grouping levels use their empirical random effects.
New levels cause an error unless `allow.new.levels = TRUE`, which
requests their population-level prediction with a zero random
contribution.

## Computational scope

Version `0.3.0` builds a dense (n n) marginal covariance matrix at each
objective evaluation. Although the random-effects design matrix is
initially sparse, covariance assembly and factorization are dense. The
package is therefore intended for small and moderate data sets.
Large-scale sparse algorithms and generalized responses are outside the
current model scope. Independent crossed and nested random-effect terms
are supported, with an unstructured covariance within each term.

Davison, A. C., and D. V. Hinkley. 1997. *Bootstrap Methods and Their
Application*. Cambridge University Press.
<https://doi.org/10.1017/CBO9780511802843>.

Fai, Alex Hrong-Tai, and Paul L. Cornelius. 1996. “Approximate f-Tests
of Multiple Degree of Freedom Hypotheses in Generalized Least Squares
Analyses of Unbalanced Split-Plot Experiments.” *Journal of Statistical
Computation and Simulation* 54 (4): 363–78.
<https://doi.org/10.1080/00949659608811740>.

Harville, David A. 1977. “Maximum Likelihood Approaches to Variance
Component Estimation and to Related Problems.” *Journal of the American
Statistical Association* 72 (358): 320–38.
<https://doi.org/10.1080/01621459.1977.10480998>.

Kenward, Michael G., and James H. Roger. 1997. “Small Sample Inference
for Fixed Effects from Restricted Maximum Likelihood.” *Biometrics* 53
(3): 983–97. <https://doi.org/10.2307/2533558>.

LaMotte, Lynn R. 2007. “A Direct Derivation of the REML Likelihood
Function.” *Statistical Papers* 48 (2): 321–27.
<https://doi.org/10.1007/s00362-006-0335-6>.

Patterson, H. D., and R. Thompson. 1971. “Recovery of Inter-Block
Information When Block Sizes Are Unequal.” *Biometrika* 58 (3): 545–54.
<https://doi.org/10.1093/biomet/58.3.545>.

Pinheiro, José C., and Douglas M. Bates. 2000. *Mixed-Effects Models in
s and s-PLUS*. Springer. <https://doi.org/10.1007/b98882>.

Satterthwaite, F. E. 1946. “An Approximate Distribution of Estimates of
Variance Components.” *Biometrics Bulletin* 2 (6): 110–14.
<https://doi.org/10.2307/3002019>.

Searle, S. R., F. M. Speed, and G. A. Milliken. 1980. “Population
Marginal Means in the Linear Model: An Alternative to Least Squares
Means.” *The American Statistician* 34 (4): 216–21.
<https://doi.org/10.1080/00031305.1980.10483031>.

Self, Steven G., and Kung-Yee Liang. 1987. “Asymptotic Properties of
Maximum Likelihood Estimators and Likelihood Ratio Tests Under
Nonstandard Conditions.” *Journal of the American Statistical
Association* 82 (398): 605–10.
<https://doi.org/10.1080/01621459.1987.10478472>.
