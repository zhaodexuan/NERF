#' @title Non-Effortful Responses Filter
#'
#' @description Detect non-effortful responses in questionnaires or examinations.
#'
#' @param theTempData The data matrix to analysis.
#'
#' @param ifInf = 'NA' 'edge'
#'
#' @param ifNA = 'omit' 'knn' 'central'
#'
#' @param k = 10
#'
#' @param center = TRUE
#'
#' @param scale = TRUE
#'
#' @param fm = 'mle'
#'
#' @param rotate = 'varimax'
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
#' theNERF <- NERF(theTempData)
#' }
#'
#' @export
#'

NERF <- function(theTempData, ifInf = 'edge', ifNA = 'central', k=10,
                 center = TRUE, scale = TRUE, fm = 'mle', rotate = 'varimax',
                 theModel = 'GPCM', pack ='mirt', SE = FALSE, precision = 4,
                 N_nodes = 30, max_outer = 60, max_inner = 60, tol = 0.001,
                 method = 'EM', theStep = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                 maxBoot = 1000, theCompare = 'upper', theAlt = 'less', theSig = 0.05,
                 ifItem = FALSE, theCompareItem = 'upper', theAltItem = 'less', theSigItem = 0.05){

  theData <- getTheData(theTempData, ifInf = ifInf, ifNA = ifNA, k = k)

  theDim <- getTheDim(theData, center = center, scale = scale, fm = fm, rotate = rotate)

  theCategory <- getTheCategory(theData)

  theExpectPoint <- getTheExpectPoint(theData, theCategory, theDim, theModel = theModel,
                                      pack = pack, SE = SE, precision = precision, N_nodes = N_nodes,
                                      max_outer = max_outer, max_inner = max_inner, tol = tol,
                                      method = method)

  # theDeviation <- getTheDeviation(theData, theExpectPoint, theCategory, theDim,
  #                                 RandomReplaceRatio = RandomReplaceRatio, maxBoot = maxBoot,
  #                                 theCompare = theCompare, theAlt = theAlt, theSig = theSig)

  thePercent <- getThePercent(theData, theExpectPoint, theCategory, theDim,
                              theStep = theStep, maxBoot = maxBoot,
                              theCompare = theCompare, theAlt = theAlt, theSig = theSig, ifItem = ifItem,
                              theCompareItem = theCompareItem, theAltItem = theAltItem, theSigItem = theSigItem)

  theNERF <- list()
  theNERF[[1]] <- theData
  theNERF[[2]] <- theDim
  theNERF[[3]] <- theCategory
  theNERF[[4]] <- theExpectPoint
  theNERF[[5]] <- thePercent
  names(theNERF) <- c('theData','theDim','theCategory','theExpectPoint','thePercent')
  return(theNERF)
}
