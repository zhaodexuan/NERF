#' @title Estimate the PCM
#'
#' @description Get the estimation of the PCM.
#'
#' @param theDataK The data matrix of factor K.
#'
#' @param theC The category of factor K.
#'
#' @param SE = FALSE
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
#' @importFrom eRm PCM
#'
#' @importFrom eRm person.parameter
#'
#' @importFrom eRm pmat
#'

getThePCM <- function(theDataK, theC, SE = FALSE){
  res <- PCM(theDataK, se = SE)  #eRm
  pres <- person.parameter(res)
  thePMatrix <- pmat(pres)
  cName <- colnames(thePMatrix)
  rName <- rownames(thePMatrix)

  theProb <- array(NA, dim = c(length(theDataK[, 1]), length(theDataK[1, ]), max(theC)+1))
  for(i in 1:length(thePMatrix[, 1])){
    for (j in 1:length(thePMatrix[1, ])) {
      theRow <- as.numeric(substr(rName[i], 2, nchar(rName[i])))
      theCD <- unlist(strsplit(cName[j], split = '.c'))
      theCol <- as.numeric(substr(theCD[1], 2, nchar(theCD[1])))
      theDep <- as.numeric(theCD[2])
      theProb[theRow,theCol,(theDep+1)] <- thePMatrix[i,j]
    }
  }

  for(i in 1:length(theDataK[, 1])){
    for (j in 1:length(theDataK[1, ])) {
      if(!sum(is.na(theProb[i,j,2:(theC[j]+1)]))==theC[j]){
        theProb[i,j,1] <- 1 - sum(theProb[i,j,2:(theC[j]+1)],na.rm = TRUE)
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
