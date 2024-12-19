# Version 0.6.1

- Fix: error in `spatial_hetero_test.hgwrm()` when there is only one SLR effect.

# Version 0.6.0

- Feature: F-test for GLSW effects
- Feature: `print_table_md` add `booktabs` table format
- Change: variable name style
- Change: spatial heterogeneity test now use bootstrap method
- Change: a new data set `multisampling`; the old one is renamed to `mulsam.test`.
- Fix: problems caused by Bisquared kernel
- Fix: lower bound of bandwidth optimising with bisquared kernel
- Fix: effect names
- Fix: unwanted output of print.hgwr
- Fix: calculation of standard error of GLSW effects
- Fix: errors in S3 methods when slop or intercept is missing
- Fix: comapre of intercept
- Fix: update Ygf before final calibration
- Fix: deal with random effect only 1
- Fix: add a space above tables
- Fix: pass args when printing F test

# Version 0.5.0

- Feature: `spatial_hetero_test()` function for testing spatial heterogeneity
- Feature: show progress bar when testing spatial heterogeneity of GLSW effects
- Feature: support Gaussian kernel function when testing spatial heterogeneity
- Improve: efficiency of testing spatial heterogeneity
- Fix: bandwidth optimisation issues
- Fix: output message of summary of `hgwrm` objects

# Version 0.4.0

- Feature: support use `L()` to specify local fixed effects
- Feature: support `sf` objects
- Feature: enable set no intercept
- Feature: report more diagnostic information
- Update: `wuhan.hp` data

# Version 0.3.0

- Feature: bandwidth optimisation with CV criterion
- Update: `wuhan.hp` data
- Fix: fill parameter to `cat` function

# Version 0.2.3

- Feature: new data sets
- Fix: rename `multisampling$coord` to `coords`
- Fix: data document errors
- Fix: error in `wuhan.hp` data
- Fix: use group index instead of original value
- Fix: don't run examples of large data set

# Version 0.2.2

- Fix: use `gsl-config --libs` in `configure` script

# Version 0.2.1

- Fix: C++ compile errors
- Fix: first sentence of description field
- Fix: error in R function `coef.hgwrm`
- Fix: missing \value in .Rd files

# Version 0.2.0

- Feature: R interface to set kernel
- Fix: incompatible generic function definitions
- Fix: other R check warnings

# Version 0.1.0

- Feature: R function hgwr
