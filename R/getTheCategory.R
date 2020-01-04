#' @title Calculate the Category
#'
#' @description Calculate the category of the data.
#'
#' @param theData The prepared data matrix.
#'
#' @return The function returns an array.
#'
#' @author zdx, \email{zhaodexuan@aliyun.com}
#'
#' @examples
#' \dontrun{
#' theCategory <- getTheCategory(theData)
#' }
#'
#' @export
#'

getTheCategory <- function(theData){
  theC <- c(1:length(theData[1, ]))

  for (j in 1:length(theData[1, ])) {
    theC[j] <- max(theData[, j])
  }

  return(theC)
}
