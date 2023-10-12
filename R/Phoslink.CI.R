#' Return causal inference based on small sample size
#'
#' @description The function returns a data frame with the IVW estimate, p-value, and confidence intervals
#' @param omicdata1 Exposure dataframe
#' @param omicdata2 Outcome dataframe
#' @param Priori priori evidence
#' @param ratio SNP minor allele frequency
#' @param Covariates Common covariates include age, gender, and smoking status
#' @return causal inference
#' @examples
#' @export

Phoslink.CI <- function(omicdata1, omicdata2, SNP, Priori, ratio = 0.1, Covariates) {
  # intersect sample names between omicdata1, omicdata2, and SNP
  intersample <- intersect(intersect(colnames(omicdata1), colnames(omicdata2)), colnames(SNP))
  # Check if intersample is non-empty
  if (length(intersample) == 0) {
    stop("No common samples between the omicdata and SNP")
  }
  # subset omicdata1, omicdata2, and SNP based on the intersected sample names
  omicdata1 <- omicdata1[, match(intersample, colnames(omicdata1))]
  omicdata2 <- omicdata2[, match(intersample, colnames(omicdata2))]
  rownames(SNP) <- SNP[, 3]
  SNP <- SNP[, match(intersample, colnames(SNP))]
  # if Covariates is not provided, create a matrix of zeros with 3 rows (AGE, Gender, SMOKING) and the intersected sample names as columns
  if (is.null(Covariates)) {
    Covariates <- matrix(0, 3, length(intersample))
    rownames(Covariates) <- c("AGE", "Gender", "SMOKING")
    colnames(Covariates) <- intersample
  }

  # extract AGE, Gender, and SMOKING from Covariates
  AGE <- as.numeric(Covariates[2, match(intersample, colnames(Covariates))])
  Gender <- as.numeric(Covariates[1, match(intersample, colnames(Covariates))])
  SMOKING <- as.numeric(Covariates[3, match(intersample, colnames(Covariates))])
  # extract instrumental variables from SNP that are also in the Priori list
  IVs <- intersect(rownames(SNP), Priori[, 2])

  library(MendelianRandomization)
  # extract the first row of omicdata1 and omicdata2 as the exposure and outcome variables
  X <- as.numeric(omicdata1[1, ])
  Y <- as.numeric(omicdata2[1, ])
  # find the indices of missing values in X and Y
  Xisna <- which(is.na(X))
  Yisna <- which(is.na(Y))
  isna <- union(Xisna, Yisna)
  # subset SNP based on the instrumental variables and the ratio of minor allele frequency
  Z <- t(SNP[rownames(SNP) %in% IVs, , drop = F])
  if (!exists("ratio") || ratio < 0 || ratio > 1) {
    stop("Error: Ratio should be provided and be between 0 and 1.")
  }
  Z <- Z[, intersect(which(apply(Z[-Xisna, , drop = F], 2, function(x) {
    all(table(x) > dim(Z[-Xisna, ])[1] * ratio)
  })), which(apply(Z[-Yisna, , drop = F], 2, function(x) {
    all(table(x) > dim(Z[-Yisna, ])[1] * ratio)
  }))), drop = F]
  # if there is no instrumental variable left after the filtering, stop the function and return an error message
  if (ncol(Z) == 0) {
    stop("Error: No instrumental variables found for this exposure.")
  }
  # if there are instrumental variables left after the filtering, perform Mendelian randomization analysis
  if (ncol(Z) != 0) {
    J <- ncol(Z)
    N <- nrow(Z)
    betaX <- array(NA, dim = J)
    betaY <- array(NA, dim = J)
    sebetaY <- array(NA, dim = J)
    sebetaX <- array(NA, dim = J)
    PrX <- array(NA, dim = J)
    PrY <- array(NA, dim = J)
    R2X <- array(NA, dim = J)
    R2Y <- array(NA, dim = J)
    R2P <- array(NA, dim = J)
    for (iIV in 1:J) {
      if (length(Xisna) > 0) {
        regX <- lm(X[-Xisna] ~ Z[-Xisna, iIV] + AGE[-Xisna] + Gender[-Xisna] + SMOKING[-Xisna])
      } else {
        regX <- lm(X ~ Z[, iIV] + AGE + Gender + SMOKING)
      }
      if (length(Yisna) > 0) {
        regY <- lm(Y[-Yisna] ~ Z[-Yisna, iIV] + AGE[-Yisna] + Gender[-Yisna] + SMOKING[-Yisna])
      } else {
        regY <- lm(Y ~ Z[, iIV] + AGE + Gender + SMOKING)
      }

      betaX[iIV] <- summary(regX)$coefficients[2, 1]
      sebetaX[iIV] <- summary(regX)$coefficients[2, 2]
      PrX[iIV] <- summary(regX)$coefficients[2, 4]
      betaY[iIV] <- summary(regY)$coefficients[2, 1]
      sebetaY[iIV] <- summary(regY)$coefficients[2, 2]
      PrY[iIV] <- summary(regY)$coefficients[2, 4]
    }
    whichIV <- which(PrX < 0.05)
    if (length(whichIV) == 0) {
      stop("Error: No instrumental variables found for this exposure.")
    }
    if (length(whichIV) != 0) {
      MR_single <- as.data.frame(matrix(, 1, 7))
      colnames(MR_single) <- c("IV", "Exposure", "Outcome", "ivw_Estimate", "ivw_Pvalue", "ivw_CILower", "ivw_CIUpper")
      oggetto <- mr_input(
        bx = as.numeric(betaX[whichIV]),
        bxse = as.numeric(sebetaX[whichIV]),
        by = as.numeric(betaY[whichIV]),
        byse = as.numeric(sebetaY[whichIV]),
        exposure = "exposure", outcome = "outcome",
        snps = colnames(Z)[whichIV]
      )
      MRresult_ivw <- mr_ivw(oggetto, model = "default", weights = "delta", psi = cor.test(X, Y)$estimate)
      MR_single[1, 1:7] <- c(
        paste(colnames(Z)[whichIV], collapse = "; "), rownames(omicdata1), rownames(omicdata2),
        MRresult_ivw$Estimate, MRresult_ivw$Pvalue, MRresult_ivw$CILower, MRresult_ivw$CIUpper
      )
    }
  }
  return(MR_single)
}
