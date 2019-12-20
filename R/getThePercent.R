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
                          theCompare = 'upper', theAlt = 'less', theSig = 0.05){
  nFac <- max(theDim)
  thePR <- list()
  for (s in 1:length(theStep)) {
    theDeviation <- getTheDeviation(theData, theExpectPoint, theCategory, theDim,
                                    RandomReplaceRatio = theStep[s], maxBoot = maxBoot, theCompare = theCompare,
                                    theAlt = theAlt, theSig = theSig)

    theOutABS <- list()
    theOutRMSD <- list()
    theOutWABS <- list()
    theOutWRMSD <- list()
    for (k in 1:nFac) {
      ABS <- unlist(theDeviation[[1]]) # nTestee, 1observe_2mean_3lower_4upper_5sign, nFac
      RMSD <- unlist(theDeviation[[2]])
      WABS <- unlist(theDeviation[[3]])
      WRMSD <- unlist(theDeviation[[4]])

      theFacABS <- c()
      for (i in 1:length(theData[,1])) {
        if(!is.na(ABS[i,5,k])){
          theFacABS <- c(theFacABS,i)
        }
      }

      theFacRMSD <- c()
      for (i in 1:length(theData[,1])) {
        if(!is.na(RMSD[i,5,k])){
          theFacRMSD <- c(theFacRMSD,i)
        }
      }

      theFacWABS <- c()
      for (i in 1:length(theData[,1])) {
        if(!is.na(WABS[i,5,k])){
          theFacWABS <- c(theFacWABS,i)
        }
      }

      theFacWRMSD <- c()
      for (i in 1:length(theData[,1])) {
        if(!is.na(WRMSD[i,5,k])){
          theFacWRMSD <- c(theFacWRMSD,i)
        }
      }

      theOutABS[[k]] <- theFacABS
      theOutRMSD[[k]] <- theFacRMSD
      theOutWABS[[k]] <- theFacWABS
      theOutWRMSD[[k]] <- theFacWRMSD
    }

    theStepPR <- list()
    theStepPR[[1]] <- theOutABS
    theStepPR[[2]] <- theOutRMSD
    theStepPR[[3]] <- theOutWABS
    theStepPR[[4]] <- theOutWRMSD
    theStepPR[[5]] <- theDeviation
    names(theStepPR) <- c('outABS','outRMSD','outWABS','outWRMSD','theDeviation')
    thePR[[s]] <- theStepPR
  }

  names(thePR) <- c(paste0('percent', theStep))
  return(thePR)
}
