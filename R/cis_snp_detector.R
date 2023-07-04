#' Identifying cis-regulatory SNPs
#'
#' @description This function identifies cis-regulatory SNPs within a certain distance from genes of interest.
#' @param genelist
#' @param SNP
#' @param cisDist
#' @param GencodeAnnotation

#' @return cis-regulatory SNPs
#' @export

cis_snp_detector <- function(genelist, SNP, cisDist, GencodeAnnotation) {
  Totalproteinloc <- GencodeAnnotation[GencodeAnnotation[, 2] %in% genelist, ]

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
