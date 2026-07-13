# ML, REML, and Satterthwaite inference

## Covariance model

For covariance parameters $`\eta`$, `lmmix` constructs

``` math
V(\eta) = ZG(\eta)Z^\mathsf{T} + R(\eta).
```

This decomposition and the use of structured residual covariance
matrices follow the general Gaussian mixed-model framework described by
Harville ([1977](#ref-harville1977)) and Pinheiro and Bates
([2000](#ref-pinheiro2000)).

Variances and Cholesky diagonal elements are represented on a log scale.
AR(1) correlations use the hyperbolic tangent. Compound symmetry uses a
scaled logistic transformation that respects its dimension-dependent
lower bound. Toeplitz correlations are obtained from partial
autocorrelations, and unstructured matrices use a Cholesky factor.

## Profiled likelihoods

For a given $`V`$, the generalized least-squares estimator is

``` math
\widehat\beta =
(X^\mathsf{T}V^{-1}X)^{-1}X^\mathsf{T}V^{-1}Y.
```

The profiled negative ML log-likelihood is

``` math
\frac{1}{2}\left\{
n\log(2\pi) + \log|V| +
(Y-X\widehat\beta)^\mathsf{T}V^{-1}(Y-X\widehat\beta)
\right\}.
```

The REML criterion adds $`\log|X^\mathsf{T}V^{-1}X|`$ and replaces $`n`$
by $`n-p`$ in the normalizing constant. REML originates with Patterson
and Thompson ([1971](#ref-patterson1971)); direct likelihood derivations
are given by Harville ([1977](#ref-harville1977)) and LaMotte
([2007](#ref-lamotte2007)). Cholesky solves are used instead of explicit
matrix inverses.

## Satterthwaite degrees of freedom

For a contrast row $`L`$, define

``` math
q(\eta) = L\operatorname{Var}(\widehat\beta)L^\mathsf{T}.
```

If $`g = \partial q/\partial\eta`$ and $`\widehat\Sigma_\eta`$ is the
inverse observed-information matrix, the implemented approximation is
the Satterthwaite approximation ([Satterthwaite
1946](#ref-satterthwaite1946)):

``` math
\nu = \frac{2q^2}{g^\mathsf{T}\widehat\Sigma_\eta g}.
```

The gradient and likelihood Hessian are calculated with `numDeriv`. For
a multi-row hypothesis, the covariance of $`L\widehat\beta`$ is
decomposed into orthogonal one-dimensional contrasts and combined with
the Fai-Cornelius moment-matching formula ([Fai and Cornelius
1996](#ref-fai1996)).

## Random effects

After covariance estimation, the empirical BLUP is

``` math
\widehat u = GZ^\mathsf{T}V^{-1}(Y-X\widehat\beta).
```

This empirical best linear unbiased predictor follows the mixed-model
derivation in Harville ([1977](#ref-harville1977)).

Fai, Alex Hrong-Tai, and Paul L. Cornelius. 1996. “Approximate f-Tests
of Multiple Degree of Freedom Hypotheses in Generalized Least Squares
Analyses of Unbalanced Split-Plot Experiments.” *Journal of Statistical
Computation and Simulation* 54 (4): 363–78.
<https://doi.org/10.1080/00949659608811740>.

Harville, David A. 1977. “Maximum Likelihood Approaches to Variance
Component Estimation and to Related Problems.” *Journal of the American
Statistical Association* 72 (358): 320–38.
<https://doi.org/10.1080/01621459.1977.10480998>.

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
