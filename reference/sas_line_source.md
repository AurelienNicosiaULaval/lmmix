# Line-source irrigation example from the SAS MIXED documentation

Winter-wheat yields for three cultivars in three blocks. Each plot
contains twelve subplots arranged on two sides of a line-source
sprinkler. SAS uses these data in Example 79.6 of the SAS/STAT 14.3
documentation to combine three random-intercept terms with a four-band
Toeplitz residual covariance.

## Usage

``` r
sas_line_source
```

## Format

A data frame with 108 rows and 6 variables:

- Block:

  Block identifier with three levels.

- Cult:

  Winter-wheat cultivar with three levels.

- Sbplt:

  Ordered subplot position from 1 to 12.

- Irrig:

  Irrigation level with six levels.

- Dir:

  Side of the sprinkler with levels North and South.

- Y:

  Continuous yield response.

## Source

SAS Institute Inc. (2017). *SAS/STAT 14.3 User's Guide*, Example 79.6,
Line-Source Sprinkler Irrigation.
<https://support.sas.com/documentation/onlinedoc/stat/>
