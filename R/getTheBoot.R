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

getTheBoot <- function(theDataK, thePoint, theC, RandomReplaceRatio = 0.5, maxBoot = 200, tempPrint = '', ifItem = FALSE){
  theS <- array(NA, dim = c(length(thePoint[, 1]), 4, maxBoot))
  colnames(theS) <- c('MAD', 'RMSD', 'WMAD', 'WRMSD')
  theS_item <- array(0, dim = c(length(thePoint[, 1]), length(thePoint[1, ]), maxBoot))

  for (b in 1:maxBoot) {
    tempPrintB <- paste0(tempPrint, 'boot', b, '/', maxBoot)
    print(tempPrintB)

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

    if(ifItem){
      for (i in 1:length(thePoint[, 1])) {
        for (j in 1:length(thePoint[1, ])) {
          if(!is.na(thePoint[i,j]) && !is.na(tempDataK[i,j])){
            theS_item[i,j,b] <- thePoint[i,j]-tempDataK[i,j]

            # theS_item[i,j,1] <- theS_item[i,j,1] + abs(thePoint[i,j]-tempDataK[i,j])
            # theS_item[i,j,2] <- theS_item[i,j,2] + (thePoint[i,j]-tempDataK[i,j])^2
            # theS_item[i,j,3] <- theS_item[i,j,3] + abs((thePoint[i,j]-tempDataK[i,j])/theC[j])
            # theS_item[i,j,4] <- theS_item[i,j,4] + ((thePoint[i,j]-tempDataK[i,j])/theC[j])^2
          }
        }
      }
    }

  }

  # if(ifItem){
  #   for (i in 1:length(thePoint[, 1])) {
  #     for (j in 1:length(thePoint[1, ])) {
  #       if(!is.na(thePoint[i,j]) && !is.na(tempDataK[i,j])){
  #         theS_item[i,j,1] <- theS_item[i,j,1] / maxBoot
  #         theS_item[i,j,2] <- sqrt(theS_item[i,j,2] / maxBoot)
  #         theS_item[i,j,3] <- theS_item[i,j,3] / maxBoot
  #         theS_item[i,j,4] <- sqrt(theS_item[i,j,4] / maxBoot)
  #       }else{
  #         theS_item[i,j,1] <- NA
  #         theS_item[i,j,2] <- NA
  #         theS_item[i,j,3] <- NA
  #         theS_item[i,j,4] <- NA
  #       }
  #     }
  #   }
  # }

  theSS <- list()
  theSS[[1]] <- theS
  theSS[[2]] <- theS_item
  return(theSS)
}
