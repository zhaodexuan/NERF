#' @title Estimate the Expectations
#'
#' @description Estimate the expectations of each person on each item.
#'
#' @param theData The prepared data matrix.
#'
#' @param theCategory The Calculated Category.
#'
#' @param theDim The Calculated Dimensions.
#'
#' @param theModel = 'GGUM' 'GPCM' 'PCM' 'GRM'
#'
#' @param pack = 'GGUM' 'mirt'
#'
#' @param SE = FALSE
#'
#' @param precision = 4
#'
#' @param N_nodes = 30
#'
#' @param max_outer = 60
#'
#' @param max_inner = 60
#'
#' @param tol = 0.001
#'
#' @param method = 'EM'
#'
#' @return The function returns a list.
#'
#' @author zdx, \email{zhaodexuan@aliyun.com}
#'
#' @examples
#' \dontrun{
#' theExpectPoint <- getTheExpectPoint(theData, theCategory, theDim)
#' }
#'
#' @export
#'

getTheExpectPoint <- function(theData, theCategory, theDim, theModel = 'GGUM', pack = 'GGUM',
                              SE = FALSE, precision = 4, N_nodes = 30, max_outer = 60,
                              max_inner = 60, tol = 0.001, method = 'EM'){
  nItem <- length(theData[1, ])
  nTestee <- length(theData[, 1])
  nFac <- max(theDim)

  theProbArray <- array(NA, dim = c(nTestee, nItem, (max(theCategory)+1)))
  theExpPoint <- array(NA, dim = c(nTestee, nItem))
  for (k in 1:nFac) {
    theDataK <- theData[, theDim==k]

    if(length(theDataK[1, ])>2){
      theC <- theCategory[theDim==k]
      theI <- length(theC)

      if(theModel=='GGUM'){
        thePoint <- getTheGGUM(theDataK, theC, pack = pack, SE = SE, precision = precision,
                               N_nodes = N_nodes, max_outer = max_outer, max_inner = max_inner,
                               tol = tol, method = method)
      }else if(theModel=='PCM'){
        thePoint <- getThePCM(theDataK, theC, SE = SE)
      }else if(theModel=='GPCM'){
        thePoint <- getTheGPCM(theDataK, theC, SE = SE)
      }else if(theModel=='GRM'){
        thePoint <- getTheGRM(theDataK, theC, SE = SE)
      }

      theProbArray[, theDim==k, ] <- thePoint[[1]]
      theExpPoint[, theDim==k] <- thePoint[[2]]
    }
  }

  theList <- list()
  theList[[1]] <- theProbArray
  theList[[2]] <- theExpPoint
  names(theList) <- c('probability', 'expectation')
  return(theList)
}
