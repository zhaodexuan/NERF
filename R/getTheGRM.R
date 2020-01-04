#' @title Estimate by GRM
#'
#' @description Get the estimation of the GRM.
#'
#' @param theDataK The data matrix of factor K.
#'
#' @param theC The category of factor K.
#'
#' @param SE = FALSE
#'
#' @return The function returns a list.
#'
#' @author zdx, \email{zhaodexuan@aliyun.com}
#'
#' @examples
#' \dontrun{
#' NA
#' }
#'
#' @importFrom ltm grm
#'
#' @importFrom ltm factor.scores
#'

getTheGRM <- function(theDataK, theC, SE = FALSE){
  ltmGRM <- grm(theDataK)
  ltmTheta <- factor.scores(ltmGRM, resp.patterns=theDataK, method="EAP")
  ltmSummary <- summary(ltmGRM)
  theA <- array(NA, dim = c(length(theC), 1))
  theB <- array(NA, dim = c(length(theC), max(theC)))
  theT <- ltmTheta$score.dat[(length(theC)+3)]
  theT <- as.array(theT[[1]])
  for (j in 1:length(theC)) {
    # ltmCoef <- ltmGRM$coefficients[[j]]
    ltmCoef <- ltmSummary$coefficients[[j]]
    for (k in 1:theC[j]) {
      theB[j,k] <- ltmCoef[k]
    }

    theA[j] <- ltmCoef[(theC[j]+1)]
  }

  theProb <- array(0,dim = c(length(theT),length(theA),(max(theC)+1)))
  theProb[,,1] <- 1
  for (i in 1:length(theT)) {
    for (j in 1:length(theA)) {
      for (k in 1:theC[j]) {
        # theProb[i,j,k+1] <- 1/(1+exp(-1*(theA[j]*theT[i]-theB[j,k])))
        theProb[i,j,k+1] <- 1/(1+exp(-1*theA[j]*(theT[i]-theB[j,k])))
        theProb[i,j,k] <- theProb[i,j,k]-theProb[i,j,k+1]
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

  theList <- list()
  theList[[1]] <- theProb
  theList[[2]] <- thePoint
  return(theList)
}
