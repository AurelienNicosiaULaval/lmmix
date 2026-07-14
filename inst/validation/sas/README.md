# Official SAS PROC MIXED validation targets

This directory records validation inputs and numerical targets from three
official examples in Chapter 79 of the SAS/STAT 14.3 User's Guide:

- Example 79.1, Split-Plot Design;
- Example 79.2, Repeated Measures;
- Example 79.5, Random Coefficients.

Primary source:

SAS Institute Inc. (2017). *SAS/STAT 14.3 User's Guide: The MIXED
Procedure*. Cary, NC: SAS Institute Inc.
<https://support.sas.com/documentation/onlinedoc/stat/>

The source document identifier is `statug/14.3/content/mixed.pdf`, Chapter 79.

The CSV file contains values transcribed from the published output tables.
The SAS program is a compact runnable reproduction of the model
specifications. These artifacts are documentation-based regression targets,
not output from a SAS session executed in this repository. A future archived
SAS log can be checked against the same CSV schema without changing the R
tests.

The source manual itself is not redistributed. The package includes only the
small example data, model specifications, and numerical values required for
reproducible comparison.
