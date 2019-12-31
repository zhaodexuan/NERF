#' @title Step Percent of Random Replacement Deviation
#'
#' @description Step by step list the proportionally random expectation deviations.
#'
#' @param theData The prepared data matrix.
#'
#' @param theExpectPoint The expectation matrix.
#'
#' @param theCategory The Calculated Category.
#'
#' @param theDim The Calculated Dimensions.
#'
#' @param theStep An array of random replacement ratio.
#'
#' @param maxBoot The maximum steps of bootstrap.
#'
#' @param theCompare = 'between', 'lower', 'upper', 'mean'
#'
#' @param theAlt = 'two.sided', 'greater', 'less'
#'
#' @param theSig = 0.05
#'
#' @param ifItem = FALSE, calculate the items when TRUE.
#'
#' @param theCompareItem = 'between', 'lower', 'upper', 'mean'
#'
#' @param theAltItem = 'two.sided', 'greater', 'less'
#'
#' @param theSigItem = 0.05
#'
#' @return The function returns a list.
#'
#' @author zdx, \email{zhaodexuan@aliyun.com}
#'
#' @examples
#' \dontrun{
#' thePercent <- getThePercent(theData, theExpectPoint, theCategory, theDim)
#' }
#'
#' @export
#'

getThePercent <- function(theData, theExpectPoint, theCategory, theDim,
                          theStep = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), maxBoot = 1000,
                          theCompare = 'upper', theAlt = 'less', theSig = 0.05, ifItem = FALSE,
                          theCompareItem = 'upper', theAltItem = 'less', theSigItem = 0.05){

  nFac <- max(theDim)
  thePR <- list()
  for (s in 1:length(theStep)) {
    tempPrintP <- paste0('percent', theStep[s], ':')
    # print(tempPrintP)
    theDeviation <- getTheDeviation(theData, theExpectPoint, theCategory, theDim,
                                    RandomReplaceRatio = theStep[s], maxBoot = maxBoot, theCompare = theCompare,
                                    theAlt = theAlt, theSig = theSig, tempPrint = tempPrintP, ifItem = ifItem,
                                    theCompareItem = theCompareItem, theAltItem = theAltItem, theSigItem = theSigItem)

    theOutMAD <- list()
    theOutRMSD <- list()
    theOutWMAD <- list()
    theOutWRMSD <- list()
    for (k in 1:nFac) {
      MAD <- unlist(theDeviation[[1]]) # nTestee, 1observe_2mean_3lower_4upper_5sign, nFac
      RMSD <- unlist(theDeviation[[2]])
      WMAD <- unlist(theDeviation[[3]])
      WRMSD <- unlist(theDeviation[[4]])

      theFacMAD <- c()
      for (i in 1:length(theData[,1])) {
        if(!is.na(MAD[i,5,k])){
          theFacMAD <- c(theFacMAD,i)
        }
      }

      theFacRMSD <- c()
      for (i in 1:length(theData[,1])) {
        if(!is.na(RMSD[i,5,k])){
          theFacRMSD <- c(theFacRMSD,i)
        }
      }

      theFacWMAD <- c()
      for (i in 1:length(theData[,1])) {
        if(!is.na(WMAD[i,5,k])){
          theFacWMAD <- c(theFacWMAD,i)
        }
      }

      theFacWRMSD <- c()
      for (i in 1:length(theData[,1])) {
        if(!is.na(WRMSD[i,5,k])){
          theFacWRMSD <- c(theFacWRMSD,i)
        }
      }

      theOutMAD[[k]] <- theFacMAD
      theOutRMSD[[k]] <- theFacRMSD
      theOutWMAD[[k]] <- theFacWMAD
      theOutWRMSD[[k]] <- theFacWRMSD
    }

    theStepPR <- list()
    theStepPR[[1]] <- theOutMAD
    theStepPR[[2]] <- theOutRMSD
    theStepPR[[3]] <- theOutWMAD
    theStepPR[[4]] <- theOutWRMSD
    theStepPR[[5]] <- theDeviation
    names(theStepPR) <- c('outMAD','outRMSD','outWMAD','outWRMSD','theDeviation')
    thePR[[s]] <- theStepPR
  }

  names(thePR) <- c(paste0('percent', theStep))
  return(thePR)
}
