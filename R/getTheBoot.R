#' @title Random Replace Bootstrap
#'
#' @description The bootstrap of random replacement.
#'
#' @param theDataK The data matrix of factor K.
#'
#' @param thePoint The expectation matrix of factor K.
#'
#' @param theC The category of factor K.
#'
#' @param RandomReplaceRatio Decimal between 0 and 1.
#'
#' @param maxBoot The maximum steps of bootstrap.
#'
#' @return The function returns an array.
#'
#' @author zdx, \email{zhaodexuan@aliyun.com}
#'
#' @examples
#' \dontrun{
#' NA
#' }
#'
#' @importFrom stats runif
#'

getTheBoot <- function(theDataK, thePoint, theC, RandomReplaceRatio = 0.5, maxBoot = 200){
  theS <- array(NA, dim = c(length(thePoint[, 1]), 4, maxBoot))
  colnames(theS) <- c('ABS', 'RMSD', 'WABS', 'WRMSD')

  for (b in 1:maxBoot) {
    # tempDataK <- theDataK
    tempDataK <- thePoint
    for (i in 1:length(thePoint[, 1])) {
      for (j in 1:length(thePoint[1, ])) {
        # i=1
        # j=2
        if(runif(1) < RandomReplaceRatio){
          tempDataK[i, j] <- sample(0:theC[j], 1, replace = TRUE)
        }
      }
    }

    for (i in 1:length(thePoint[, 1])) {
      theI <- length(theC[!is.na(thePoint[i,])])
      theS[i, 1, b] <- sum(abs(thePoint[i, !is.na(thePoint[i,])] - tempDataK[i, !is.na(thePoint[i,])])) / theI
      theS[i, 2, b] <- sqrt(sum((thePoint[i, !is.na(thePoint[i,])] - tempDataK[i, !is.na(thePoint[i,])])^2) / theI)
      theS[i, 3, b] <- sum(abs(thePoint[i, !is.na(thePoint[i,])] - tempDataK[i, !is.na(thePoint[i,])]) / theC[!is.na(thePoint[i,])]) / theI
      theS[i, 4, b] <- sqrt(sum(((thePoint[i, !is.na(thePoint[i,])] - tempDataK[i, !is.na(thePoint[i,])]) / theC[!is.na(thePoint[i,])])^2) / theI)
    }
  }

  return(theS)
}
