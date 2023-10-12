#' Identifying cis-regulatory SNPs
#'
#' @description This function identifies cis-regulatory SNPs within a certain distance from genes of interest.
#' @param gene A character vector containing the names of the genes of interest.
#' @param SNP A data frame containing SNP information, typically including chromosome, position.
#' @param cisDist An integer representing the distance threshold within which to search for Cis-regulatory SNPs, with a default value of 1e6 (1,000,000).
#' @param GencodeAnnotation  A data frame containing gene annotation information, including chromosome, gene name, and location.

#' @return cis-regulatory SNPs
#' @export

cis_snp_detector <- function(gene, SNP, cisDist = 1e6, GencodeAnnotation) {
  Totalproteinloc <- GencodeAnnotation[GencodeAnnotation[, 2] %in% gene, ]

  CisPair <- c()
  for (i in 1:nrow(Totalproteinloc)) {
    chromStartMB <- Totalproteinloc[i, 4] - cisDist
    chromEndMB <- Totalproteinloc[i, 5] + cisDist
    Samechrom <- SNP[paste0("chr", SNP[, 1]) == Totalproteinloc[i, 3], ]
    if (any(Samechrom[, 2] >= chromStartMB & Samechrom[, 2] <= chromEndMB)) {
      CisPair <- rbind(CisPair, data.frame(Gene = Totalproteinloc[i, 2], cisSNP = Samechrom[which(Samechrom[, 2] >= chromStartMB & Samechrom[, 2] <= chromEndMB), 3]))
    }
  }
  return(CisPair)
}
