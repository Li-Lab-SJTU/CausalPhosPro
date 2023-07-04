
# CausalPhosPro

<!-- badges: start -->
<!-- badges: end -->

CausalPhosPro is an R package for causal inference of phosphoproteomics data.

## Installation

You can install the released version of CausalPhosPro from github with:

``` r
install.packages("devtools")
devtools::install_github("Li-Lab-SJTU/CausalPhosPro")
```

## Directories

`extdata/` Contains built-in data sets for the package

`man/`  help files and documentation

`R/`    R functions in scripts

`simulation` includes server and ui files for `run_dia_shiny()`, a shiny-based user interface

`tests/` Includes package tests for default parameter accuracy conducted on package build

## Usage

This is a basic example which shows you how to solve a common problem:

```{r example}
library(CausalPhosPro)
## basic example code
```

## Description of the output

The following columns are available in the CausalPhosPro output:

| Column | Description |
| ------------- | ------------- |
| IV | Instrumental variables in the causal relationship |
| Exposure | The name of the exposure variable |
| Outcome | The name of the outcome variable |
| ivw_Estimate | The causal point estimate from the MR-IVW method |
| ivw_Pvalue | P-value that is obtained when testing whether this segment should be represented by one or two states. A low p-value will result in the fitting of a second copy number state |
| ivw_CILower | The lower bound of the 95% confidence interval for the estimated causal effect |
| ivw_CIUpper | The upper bound of the 95% confidence interval for the estimated causal effect |

