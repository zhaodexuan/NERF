#' @title Estimate by GPCM
#'
#' @description Get the estimation of the GPCM.
#'
#' @param theDataK The data matrix of factor K.
#'
#' @param theC The category of factor K.
#'
#' @param SE = FALSE
#'
#' @param method = 'EM'
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
#' @importFrom mirt mirt
#'
#' @importFrom mirt fscores
#'
#' @importFrom mirt coef
#'

getTheGPCM <- function(theDataK, theC, SE = FALSE, method = 'EM'){
  theNames <- array()
  for (n in 1:length(theC)) {
    theNames[n] <- paste0('item', n)
  }
  colnames(theDataK) <- theNames

  theMirt <- mirt(theDataK,1,itemtype = 'gpcm',SE = SE,method = method)

  theMirtTh <- fscores(theMirt,method="EAPsum",full.scores = TRUE)

  theMirtCoef <- coef(theMirt)
  theAlpha <- matrix(NA,length(theC),1)
  theAK <- matrix(NA,length(theC),(max(theC)+1))
  theDK <- matrix(NA,length(theC),(max(theC)+1))
  for (n in 1:length(theC)) {
    tempCoefN <- theMirtCoef[[n]]
    theAlpha[n] <- tempCoefN[1]
    theAK[n,1:(theC[n]+1)] <- tempCoefN[2:(theC[n]+2)]
    theDK[n,1:(theC[n]+1)] <- tempCoefN[(theC[n]+3):(2*theC[n]+3)]
  }

  theProb <- array(NA, dim = c(length(theDataK[, 1]), length(theDataK[1, ]), max(theC)+1))
  for (i in 1:length(theDataK[, 1])) {
    for (j in 1:length(theDataK[1, ])) {
      for (k in 1:(theC[j]+1)) {
        fenZi <- exp(theAK[j,k]*theAlpha[j]*theMirtTh[i]+theDK[j,k])
        fenMu <- 0
        for (c in 1:(theC[j]+1)) {
          fenMu <- fenMu + exp(theAK[j,c]*theAlpha[j]*theMirtTh[i]+theDK[j,c])
        }

        theProb[i,j,k] <- fenZi / fenMu
      }
    }
  }

  thePoint <- array(NA, dim = c(length(theDataK[, 1]), length(theDataK[1, ])))
  for (i in 1:length(theDataK[, 1])) {
    for (j in 1:length(theDataK[1, ])) {
      if(!sum(is.na(theProb[i,j,1:(theC[j]+1)]))==(theC[j]+1)){
        thePoint[i, j] <- which.max(theProb[i, j, ])
        thePoint[i, j] <- thePoint[i, j] - 1
      }
    }
  }

  return(thePoint)
}


