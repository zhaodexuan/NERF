#' @title Random Replace Bootstrap
#'
#' @description The bootstrap of random replacement.
#'
#' @param theDataK The data matrix of factor K.
#'
#' @param thePointK The expectation list of factor K.
#'
#' @param theC The category of factor K.
#'
#' @param RandomReplaceRatio Decimal between 0 and 1.
#'
#' @param maxBoot The maximum steps of bootstrap.
#'
#' @param tempPrint = ''
#'
#' @param ifItem = FALSE, calculate the items when TRUE.
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

getTheBoot <- function(theDataK, thePointK, theC, RandomReplaceRatio = 0.5, maxBoot = 200, tempPrint = '', ifItem = FALSE){
  theWeight <- thePointK$probability
  thePoint <- thePointK$expectation

  theS <- array(NA, dim = c(length(thePoint[, 1]), 4, maxBoot))
  colnames(theS) <- c('MAD', 'RMSD', 'WMAD', 'WRMSD')
  theS_item <- array(0, dim = c(length(thePoint[, 1]), length(thePoint[1, ]), maxBoot))
  theS_item_WABS <- array(0, dim = c(length(thePoint[, 1]), length(thePoint[1, ]), maxBoot))

  for (b in 1:maxBoot) {
    tempPrintB <- paste0(tempPrint, 'boot', b, '/', maxBoot)
    print(tempPrintB)

    tempDataK <- thePoint
    for (i in 1:length(thePoint[, 1])) {
      for (j in 1:length(thePoint[1, ])) {
        if(runif(1) < RandomReplaceRatio){
          tempDataK[i, j] <- sample(0:theC[j], 1, replace = TRUE)
        }
      }
    }

    for (i in 1:length(thePoint[, 1])) {
      theI <- length(theC[!is.na(thePoint[i,])])
      tempD <- 1:length(tempDataK[i,])
      for (j in 1:length(tempDataK[i,])) {
        if(!is.na(thePoint[i,j]) && !is.na(tempDataK[i,j])){
          tempD[j] <- log(theWeight[i,j,(tempDataK[i,j]+1)] / theWeight[i,j,(thePoint[i,j]+1)]) * (tempDataK[i,j] - thePoint[i,j])
        }
      }

      theS[i, 1, b] <- sum(abs(thePoint[i, !is.na(thePoint[i,])] - tempDataK[i, !is.na(thePoint[i,])])) / theI
      theS[i, 2, b] <- sqrt(sum((thePoint[i, !is.na(thePoint[i,])] - tempDataK[i, !is.na(thePoint[i,])])^2) / theI)
      theS[i, 3, b] <- sum(abs(tempD), na.rm = TRUE) / theI
      theS[i, 4, b] <- sqrt(sum((tempD)^2, na.rm = TRUE) / theI)

      # theS[i, 3, b] <- sum(abs(thePoint[i, !is.na(thePoint[i,])] - tempDataK[i, !is.na(thePoint[i,])]) / theC[!is.na(thePoint[i,])]) / theI
      # theS[i, 4, b] <- sqrt(sum(((thePoint[i, !is.na(thePoint[i,])] - tempDataK[i, !is.na(thePoint[i,])]) / theC[!is.na(thePoint[i,])])^2) / theI)
    }

    if(ifItem){
      for (i in 1:length(thePoint[, 1])) {
        for (j in 1:length(thePoint[1, ])) {
          if(!is.na(thePoint[i,j]) && !is.na(tempDataK[i,j])){
            theS_item[i,j,b] <- thePoint[i,j]-tempDataK[i,j]
            theS_item_WABS[i,j,b] <- log(theWeight[i,j,(tempDataK[i,j]+1)] / theWeight[i,j,(thePoint[i,j]+1)]) * (tempDataK[i,j] - thePoint[i,j])
          }
        }
      }
    }

  }

  theSS <- list()
  theSS[[1]] <- theS
  theSS[[2]] <- theS_item
  theSS[[3]] <- theS_item_WABS
  return(theSS)
}
