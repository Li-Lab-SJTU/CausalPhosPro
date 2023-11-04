
# Phoslink
<!-- badges: start -->
<!-- badges: end -->

Phoslink is an R package for causal inference of phosphoproteomics data.

## Installation
You can install the released version of Phoslink from github with:
``` r
devtools::install_github("Li-Lab-SJTU/Phoslink")
```

## Directories
-  `R/`    R functions in scripts
-  `Simulations/` Provides the R code for simulation analyses and the simulation results are also provided.
-  `inst/extdata/` contain built-in data sets for the package
-  `man/`  help files and documentation
   
## Usage
 `cis_snp_detector`:  The function cis_snp_detector() is utilized to identify cis-acting SNPs from germline SNP data. The cis_snp_detector() function is called with the following parameters
 -  `gene`: The gene symbol of interest
 -  `SNP`: The germline SNP data
 -  `cisDist`: The threshold distance for defining cis-acting SNPs
 -  `GencodeAnnotation`: Human genome annotation from GENCODE
   
 `Phoslink.CI`:  The Phoslink.CI() function is the main function used for conducting the Phoslink analysis using causal inference. It takes the following parameters:
-  `omicdata1`: Phosphorylation data for the exposure variable X
-  `omicdata2`: Protein expression data for the outcome variable Y
-  `SNP`: Germline SNPs serving as instrumental variable candidates
-  `Priori`: External prior evidence
-  `ratio`: SNP minor allele frequency
-  `Covariates`: Covariates data


## Examples
Data loading

```{r example}
# loading Phoslink package
library(Phoslink)
# Phosphorylation Data for exposure X
omicdata1 <- read.csv(system.file("extdata","Exposure_omicdata1.csv",package="Phoslink", mustWork=TRUE),row.names=1)
# Protein epression Data for outcome Y
omicdata2 <- read.csv(system.file("extdata","Outcome_omicdata2.csv",package="Phoslink", mustWork=TRUE),row.names=1)
# Germline SNPs as instrumental variable candidates
SNP <- read.csv(system.file("extdata","SNPmatrix.csv",package="Phoslink", mustWork=TRUE))
# Covariates Data
Covariates <- read.csv(system.file("extdata","Covariates.csv",package="Phoslink", mustWork=TRUE),header = TRUE,row.names=1)
# Phosphorylation-related SNPs as external prior evidence
Priori <- unique(read.table(system.file("extdata","phosSNPs.txt",package="Phoslink", mustWork=TRUE),header = T)[,c(1,2)])
# human genome annotation
GencodeAnnotation <- read.table(system.file("extdata","gencode.v40.annotation.gene.probeMap",package="Phoslink", mustWork=TRUE),header = T)
```
Conducting the Phoslink analysis
```{r example}
# Detecting cis phosphorylation-related SNPs
CisPair<-cis_snp_detector(genelist=unlist(strsplit(rownames(omicdata1),split='p'))[1],SNP=SNP[,1:3],cisDist=1e6,GencodeAnnotation=GencodeAnnotation)
cisPriori<-Priori[paste(Priori[,1],Priori[,2])%in%paste(CisPair[,1],CisPair[,2]),]
# Running the Phoslink analysis
Phoslink.CI(omicdata1=omicdata1,omicdata2=omicdata2,SNP=SNP, cisPriori=cisPriori,ratio=ratio,Covariates=Covariates)
```

## Description of the output

The following columns are available in the Phoslink.CI output:

| Column | Description |
| ------------- | ------------- |
| IV | Instrumental variables in the causal relationship |
| Exposure | The name of the exposure variable |
| Outcome | The name of the outcome variable |
| ivw_Estimate | The causal point estimate from the MR-IVW method |
| ivw_Pvalue | P-value associated with the causal estimate |
| ivw_CILower | The lower bound of the 95% confidence interval for the estimated causal effect |
| ivw_CIUpper | The upper bound of the 95% confidence interval for the estimated causal effect |

