# Package index

## Fit and control

- [`lmm()`](lmm.md) : Fit a Gaussian linear mixed model

- [`lmm_control()`](lmm_control.md) :

  Control numerical optimization for [`lmm()`](../reference/lmm.md)

## Inference and marginal means

- [`lsmeans()`](lsmeans.md) :

  Estimated marginal means for an `lmm` model

- [`summary(`*`<lmm>`*`)`](summary.lmm.md) :

  Summarize an `lmm` fit

## Extractors

- [`fixef()`](fixef.md) : Extract fixed effects
- [`ranef()`](ranef.md) : Extract random effects
- [`VarCorr()`](VarCorr.md) : Extract covariance components

## Tidy methods

- [`tidy(`*`<lmm>`*`)`](tidy.lmm.md) :

  Tidy an `lmm` model

- [`glance(`*`<lmm>`*`)`](glance.lmm.md) :

  One-row model summary for an `lmm` model

- [`augment(`*`<lmm>`*`)`](augment.lmm.md) : Augment data with fitted
  values and residuals

## Interoperability

- [`recover_data(`*`<lmm>`*`)`](recover_data.lmm.md)
  [`emm_basis(`*`<lmm>`*`)`](recover_data.lmm.md) :

  Support for `emmeans`

## Data

- [`multicentre`](multicentre.md) : Multi-site repeated-measures
  experiment
