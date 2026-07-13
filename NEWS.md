# lmmix 0.1.0

* Added explicit profiled ML and REML estimation for Gaussian linear mixed
  models.
* Added independent, compound-symmetric, AR(1), Toeplitz, and unstructured
  residual covariance matrices.
* Added random intercepts and slopes with an unstructured covariance matrix.
* Added Satterthwaite fixed-effect inference, type III tests, estimated
  marginal means, and pairwise contrasts.
* Added standard model, broom, and `emmeans` methods.
* Added the multilocation example from Milliken and Johnson (2009, Section
  28.3) and validation tests against SAS PROC MIXED, `nlme`, `lmerTest`, and
  `mmrm`.
* Added primary theoretical references for REML, Satterthwaite inference,
  multi-degree-of-freedom tests, and estimated marginal means.
* Reworked the README, vignettes, package overview, and method documentation
  to describe the implemented API, validation evidence, and current limits.
* Added combined-model convergence tests for every residual covariance
  structure.
