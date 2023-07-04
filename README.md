
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

-  `extdata/` Contains built-in data sets for the package
-  `man/`  help files and documentation
-  `R/`    R functions in scripts
-  `simulation` includes server and ui files for `run_dia_shiny()`, a shiny-based user interface
-  `tests/` Includes package tests for default parameter accuracy conducted on package build

## Usage
 `cis_snp_detector`:  The function cis_snp_detector() is utilized to identify cis-acting SNPs from germline SNP data. The cis_snp_detector() function is called with the following parameters
 -  `gene`: The gene symbol of interest
 -  `SNP`: The germline SNP data
 -  `cisDist`: The threshold distance for defining cis-acting SNPs
 -  `GencodeAnnotation`: Human genome annotation from GENCODE
   
 `CausalPhosPro`:  The function cis_snp_detector() is utilized to identify cis-acting SNPs from germline SNP data. The cis_snp_detector() function is called with the following parameters  


## Examples
Data loading

```{r example}
# loading CausalPhosPro package
library(CausalPhosPro)
# Phosphorylation Data for exposure X
omicdata1 <- read.csv(system.file("extdata","omicdata1.csv",package="CausalPhosPro", mustWork=TRUE),row.names=1)
# Protein epression Data for outcome Y
omicdata2 <- read.csv(system.file("extdata","omicdata2.csv",package="CausalPhosPro", mustWork=TRUE),row.names=1)
# Germline SNPs as instrumental variable candidates
SNP <- read.csv(system.file("extdata","SNPmatrix.csv",package="CausalPhosPro", mustWork=TRUE))
# Covariates Data
Covariates <- read.csv(system.file("extdata","Covariates.csv",package="CausalPhosPro", mustWork=TRUE),header = TRUE,row.names=1)
# Phosphorylation-related SNPs as external prior evidence
Priori <- unique(read.table(system.file("extdata","phosSNPs.txt",package="CausalPhosPro", mustWork=TRUE),header = T)[,c(1,2)])
# human genome annotation
GencodeAnnotation <- read.table(system.file("extdata","gencode.v40.annotation.gene.probeMap",package="CausalPhosPro", mustWork=TRUE),header = T)
```
Conducting the CausalPhosPro analysis
```{r example}
# Detecting cis phosphorylation-related SNPs
CisPair<-cis_snp_detector(genelist=unlist(strsplit(rownames(omicdata1),split='p'))[1],SNP=SNP[,1:3],cisDist=1e6,GencodeAnnotation=GencodeAnnotation)
cisPriori<-Priori[paste(Priori[,1],Priori[,2])%in%paste(CisPair[,1],CisPair[,2]),]
# Running the CausalPhosPro analysis
CausalPhosPro(omicdata1=omicdata1,omicdata2=omicdata2,SNP=SNP, cisPriori=cisPriori,ratio=ratio,Covariates=Covariates)
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

