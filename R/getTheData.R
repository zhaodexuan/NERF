#' @title Prepare Data for Analysis
#'
#' @description Transform the data matrix into the specific format for analysis.
#'
#' @param theTempData The data matrix to analysis.
#'
#' @param ifInf = 'NA' 'edge'
#'
#' @param ifNA = 'omit' 'knn' 'central'
#'
#' @param k = 10
#'
#' @return The function returns an array.
#'
#' @author zdx, \email{zhaodexuan@aliyun.com}
#'
#' @examples
#' \dontrun{
#' theData <- getTheData(theTempData)
#' }
#'
#' @export
#'
#' @importFrom stats na.omit
#'
#' @importFrom DMwR knnImputation
#'
#' @importFrom DMwR centralImputation
#'

getTheData <- function(theTempData, ifInf = 'edge', ifNA = 'central', k = 10){

  # NaN <- NA
  theTempData[is.na(theTempData)] <- NA

  if(ifInf=='edge'){
    # -Inf <- min
    # Inf <- max
    for (j in 1:length(theTempData[1, ])) {
      theTempData[theTempData[!is.na(theTempData[, j]), j]==-Inf, j] <- min(theTempData[theTempData[, j]!=-Inf, j], na.rm = TRUE)
      theTempData[theTempData[!is.na(theTempData[, j]), j]==Inf, j] <- max(theTempData[theTempData[, j]!=Inf, j], na.rm = TRUE)
    }
  }else if(ifInf=='NA'){
    theTempData[theTempData==-Inf] <- NA
    theTempData[theTempData==Inf] <- NA
  }

  if(sum(is.na(theTempData)) > 0){
    if(ifNA=='omit'){
      theTempData <- na.omit(theTempData)
    }else if(ifNA=='knn'){
      theTempData <- knnImputation(theTempData, k = k)
    }else if(ifNA=='central'){
      theTempData <- centralImputation(theTempData)
    }
  }

  theData <- array(0, dim = c(length(theTempData[, 1]),length(theTempData[1, ])))
  for (j in 1:length(theTempData[1, ])) {
    theData[, j] <- factor(rank(theTempData[, j]))
  }

  theData <- theData - 1
  return(theData)
}
