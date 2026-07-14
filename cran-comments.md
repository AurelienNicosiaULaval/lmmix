## Test environments

* local macOS, R release
* GitHub Actions: macOS release
* GitHub Actions: Windows release
* GitHub Actions: Ubuntu R-devel, R-release, and R-oldrel

## R CMD check results

0 errors | 0 warnings | 1 note

The single note is the expected `New submission` note. The clean-directory
`--as-cran` check reports no package errors or warnings. The GitHub Actions
matrix also completes without package errors or warnings.

## Reverse dependencies

There are no reverse dependencies because this is a first CRAN submission.
