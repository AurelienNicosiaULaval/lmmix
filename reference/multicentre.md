# Multi-site repeated-measures experiment

An unbalanced experiment with three drugs in three centers and three
repeated measurement times. This is the multilocation repeated-measures
example described by Milliken and Johnson (2009, Section 28.3).

## Usage

``` r
multicentre
```

## Format

A data frame with 153 rows and 5 variables:

- Center:

  Center identifier with levels R, S, and T.

- Drug:

  Drug identifier with levels 1, 2, and 3.

- Subject:

  Subject identifier within each center and drug.

- Time:

  Repeated-measurement occasion with levels 1, 2, and 3.

- Y:

  Continuous response, including 28 missing values.

## Source

Milliken, G. A., and Johnson, D. E. (2009). *Analysis of Messy Data,
Volume 1: Designed Experiments* (2nd ed.), Section 28.3. Chapman and
Hall/CRC. <https://doi.org/10.1201/EBK1584883340>
