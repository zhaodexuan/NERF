#' @title Calculate the Dimensions
#'
#' @description Calculate the dimensions of the data by factor analysis.
#'
#' @param theData The prepared data matrix.
#'
#' @param center = TRUE
#'
#' @param scale = TRUE
#'
#' @param fm = 'mle'
#'
#' @param rotate = 'varimax'
#'
#' @return The function returns an array.
#'
#' @author zdx, \email{zhaodexuan@aliyun.com}
#'
#' @examples
#' \dontrun{
#' theDim <- getTheDim(theData)
#' }
#'
#' @export
#'
#' @importFrom psych fa
#'
#' @importFrom psych fa.parallel
#'

getTheDim <- function(theData, center = TRUE, scale = TRUE, fm = 'mle', rotate = 'varimax'){
  nItem <- length(theData[1, ])
  nTestee <- length(theData[, 1])
  theX <- scale(theData, center = center, scale = scale)
  theM <- fa.parallel(theX)
  nFac <- theM$nfact
  theFA <- fa(theX, nfactors = nFac, fm=fm, rotate = rotate)
  theDim <- max.col(theFA$loadings)

  return(theDim)
}
