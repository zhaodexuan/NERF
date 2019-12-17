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
#' @param theModel = 'GGUM' 'GPCM' 'PCM'
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
#' @param RandomReplaceRatio Decimal between 0 and 1.
#'
#' @param maxBoot The maximum steps of bootstrap.
#'
#' @param theSig = 0.05
#'
#' @return The function returns a list.
#'
#' @author zdx, \email{zhaodexuan@aliyun.com}
#'
#' @examples
#' \dontrun{
#' theDeviation <- NERF(theTempData)
#' }
#'
#' @export
#'

NERF <- function(theTempData, ifInf = 'edge', ifNA = 'central', k=10,
                 center = TRUE, scale = TRUE, fm = 'mle', rotate = 'varimax',
                 theModel = 'GPCM', pack ='mirt', SE = FALSE, precision = 4,
                 N_nodes = 30, max_outer = 60, max_inner = 60, tol = 0.001,
                 method = 'EM', RandomReplaceRatio = 0.5, maxBoot = 200, theSig = 0.05){

  theData <- getTheData(theTempData, ifInf = ifInf, ifNA = ifNA, k = k)

  theDim <- getTheDim(theData, center = center, scale = scale, fm = fm, rotate = rotate)

  theCategory <- getTheCategory(theData)

  theExpectPoint <- getTheExpectPoint(theData, theCategory, theDim, theModel = theModel,
                                      pack = pack, SE = SE, precision = precision, N_nodes = N_nodes,
                                      max_outer = max_outer, max_inner = max_inner, tol = tol,
                                      method = method)

  theDeviation <- getTheDeviation(theData, theExpectPoint, theCategory, theDim,
                                  RandomReplaceRatio = RandomReplaceRatio,
                                  maxBoot = maxBoot, theSig = theSig)

  return(theDeviation)
}
