# Random-coefficients example from the SAS MIXED documentation

Pharmaceutical stability measurements for three batches at six shelf
ages. SAS uses these simulated data in Example 79.5 of the SAS/STAT 14.3
documentation for the MIXED procedure.

## Usage

``` r
sas_random_coefficients
```

## Format

A data frame with 108 rows and 4 variables:

- Batch:

  Batch identifier with three levels.

- Month:

  Shelf age in months.

- Replicate:

  Replicate identifier within a batch and month.

- Y:

  Assay response, including 24 missing values.

## Source

SAS Institute Inc. (2017). *SAS/STAT 14.3 User's Guide*, Example 79.5,
Random Coefficients.
<https://support.sas.com/documentation/onlinedoc/stat/>
