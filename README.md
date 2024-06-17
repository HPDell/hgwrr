# hgwrr

[![CRAN](https://www.r-pkg.org/badges/version/hgwrr)](https://cran.r-project.org/package=hgwrr)
[![Documentation](https://img.shields.io/badge/Documentation-blue)](https://hpdell.github.io/hgwrr)
[![R Package Check](https://github.com/HPDell/hgwrr/actions/workflows/R.yml/badge.svg)](https://github.com/HPDell/hgwrr/actions/workflows/R.yml)

This package is provides R interfaces to calibrate the Hierarchical and Geographically Weighted Regression (HGWR) model.

## Installation

The package now is on CRAN.

```R
install.packages("hgwrr")
```

If you want to install the latest version from GitHub, note that the packages relies on a submodule from [hpdell/hgwr](https://github.com/hpdell/hgwr). So `install_github()` from the **devtools** package is probably not working. Instead, it's better to recursively clone the package.

```bash
git clone --recursive https://github.com/hpdell/hgwrr
R CMD INSTALL hgwrr
```

## Basic Usage

Here is a quick example showing how it works.

```r
library(hgwrr)
data(multisampling)
hgwr(
  formula = y ~ L(g1 + g2) + x1 + (z1 | group),
  data = multisampling$data,
  coords = multisampling$coords,
  bw = 10
)
```

For further information, please read this [article](https://hpdell.github.io/hgwrr/articles/introduction.html).
There is a full example.

## Reference

- Hu, Yigong, Lu, Binbin, Ge, Yong, Dong, Guanpeng, 2022. Uncovering spatial heterogeneity in real estate prices via combined hierarchical linear model and geographically weighted regression. Environment and Planning B: Urban Analytics and City Science. [DOI](https://journals.sagepub.com/doi/10.1177/23998083211063885)
- Yigong Hu, Richard Harris, Richard Timmerman, and Binbin Lu. A Hierarchical and Geographically Weighted Regression Model and Its Backfitting Maximum Likelihood Estimator (Short Paper). In 12th International Conference on Geographic Information Science (GIScience 2023). Leibniz International Proceedings in Informatics (LIPIcs), Volume 277, pp. 39:1-39:6, Schloss Dagstuhl – Leibniz-Zentrum für Informatik (2023)
  [DOI](https://doi.org/10.4230/LIPIcs.GIScience.2023.39)
