
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
omicdata1 <- read.csv(system.file("extdata","omicdata1.csv",package="CausalPhosPro", mustWork=TRUE),header = TRUE,fill=T,stringsAsFactors = FALSE,row.names=1)
omicdata1[1,]<-log2(as.numeric(omicdata1[1,]))
omicdata2 <- read.csv(system.file("extdata","omicdata2.csv",package="CausalPhosPro", mustWork=TRUE),header = TRUE,fill=T,stringsAsFactors = FALSE,row.names=1)
omicdata2[1,]<-log2(as.numeric(omicdata2[1,]))
 
SNP <- read.csv(system.file("extdata","SNPmatrix.csv",package="CausalPhosPro", mustWork=TRUE),header = TRUE,fill=T,stringsAsFactors = FALSE,check.names=F)
Priori <- unique(read.table(system.file("extdata","phosSNPs.txt",package="CausalPhosPro", mustWork=TRUE), sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE)[,c(1,2)])

GencodeAnnotation <- read.table(system.file("extdata","gencode.v40.annotation.gene.probeMap",package="CausalPhosPro", mustWork=TRUE),sep = "\t",header = TRUE,fill=T,stringsAsFactors = FALSE,check.names=F)

Covariates <- read.csv(system.file("extdata","Covariates.csv",package="CausalPhosPro", mustWork=TRUE),header = TRUE,fill=T,stringsAsFactors = FALSE,check.names=F,row.names=1)

CisPair<-cis_snp_detector(unlist(strsplit(rownames(omicdata1),split='p'))[1],SNP[,1:3],1e6,GencodeAnnotation)
cisPriori<-Priori[paste(Priori[,1],Priori[,2])%in%paste(CisPair[,1],CisPair[,2]),]
ratio=0.1
CausalPhosPro(omicdata1,omicdata2,SNP, cisPriori,ratio,Covariates)
```

## Description of the output

The following columns are available in the CausalPhosPro output:

| Column | Description |
| ------------- | ------------- |
| IV | Instrumental variables in the causal relationship |
| Exposure | The name of the exposure variable |
| Outcome | The name of the outcome variable |
| ivw_Estimate | The causal point estimate from the MR-IVW method |
| ivw_Pvalue | P-value associated with the causal estimate |
| ivw_CILower | The lower bound of the 95% confidence interval for the estimated causal effect |
| ivw_CIUpper | The upper bound of the 95% confidence interval for the estimated causal effect |

