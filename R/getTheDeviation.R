#' @title Calculate the Deviations
#'
#' @description Calculate deviations between observations and expectations.
#'
#' @param theData The prepared data matrix.
#'
#' @param theExpectPoint The expectation matrix.
#'
#' @param theCategory The Calculated Category.
#'
#' @param theDim The Calculated Dimensions.
#'
#' @param RandomReplaceRatio Decimal between 0 and 1.
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
#' theDeviation <- getTheDeviation(theData, theExpectPoint, theCategory, theDim)
#' }
#'
#' @export
#'
#' @importFrom stats quantile
#'
#' @importFrom stats t.test
#'

getTheDeviation <- function(theData, theExpectPoint, theCategory, theDim,
                            RandomReplaceRatio = 0.5, maxBoot = 200, theCompare = 'upper',
                            theAlt = 'less', theSig = 0.05){
  
  
  nTestee <- length(theData[, 1])
  nFac <- max(theDim)
  theABS <- array(NA, dim = c(nTestee, nFac))
  theRMSD <- array(NA, dim = c(nTestee, nFac))
  theWABS <- array(NA, dim = c(nTestee, nFac))
  theWRMSD <- array(NA, dim = c(nTestee, nFac))
  theRandomABS <- array(NA, dim = c(nTestee, 5, nFac)) #1observe 2mean 3lower 4upper 5sign
  theRandomRMSD <- array(NA, dim = c(nTestee, 5, nFac))
  theRandomWABS <- array(NA, dim = c(nTestee, 5, nFac))
  theRandomWRMSD <- array(NA, dim = c(nTestee, 5, nFac))
  for (k in 1:nFac) {
    theC <- theCategory[theDim==k]
    theDataK <- theData[, theDim==k]
    thePoint <- theExpectPoint[, theDim==k]

    for (i in 1:nTestee){
      theI <- length(theC[!is.na(thePoint[i,])])
      theABS[i, k] <- sum(abs(thePoint[i, !is.na(thePoint[i,])] - theDataK[i, !is.na(thePoint[i,])])) / theI
      theRMSD[i, k] <- sqrt(sum((thePoint[i, !is.na(thePoint[i,])] - theDataK[i, !is.na(thePoint[i,])])^2) / theI)
      theWABS[i, k] <- sum(abs(thePoint[i, !is.na(thePoint[i,])] - theDataK[i, !is.na(thePoint[i,])]) / theC[!is.na(thePoint[i,])]) / theI
      theWRMSD[i, k] <- sqrt(sum(((thePoint[i, !is.na(thePoint[i,])] - theDataK[i, !is.na(thePoint[i,])]) / theC[!is.na(thePoint[i,])])^2) / theI)
    }

    theS <- getTheBoot(theDataK, thePoint, theC, RandomReplaceRatio = RandomReplaceRatio, maxBoot = maxBoot)

    for (i in 1:nTestee){
      # theRandomABS[i, 2, k] <- mean(theS[i, 1, ], na.rm = TRUE)
      # theRandomRMSD[i, 2, k] <- mean(theS[i, 2, ], na.rm = TRUE)
      # theRandomWABS[i, 2, k] <- mean(theS[i, 3, ], na.rm = TRUE)
      # theRandomWRMSD[i, 2, k] <- mean(theS[i, 4, ], na.rm = TRUE)
      #
      # theRandomABS[i, 3, k] <- quantile(theS[i, 1, ],(theSig/2), na.rm = TRUE)
      # theRandomRMSD[i, 3, k] <- quantile(theS[i, 2, ],(theSig/2), na.rm = TRUE)
      # theRandomWABS[i, 3, k] <- quantile(theS[i, 3, ],(theSig/2), na.rm = TRUE)
      # theRandomWRMSD[i, 3, k] <- quantile(theS[i, 4, ],(theSig/2), na.rm = TRUE)
      #
      # theRandomABS[i, 4, k] <- quantile(theS[i, 1, ],(1-theSig/2), na.rm = TRUE)
      # theRandomRMSD[i, 4, k] <- quantile(theS[i, 2, ],(1-theSig/2), na.rm = TRUE)
      # theRandomWABS[i, 4, k] <- quantile(theS[i, 3, ],(1-theSig/2), na.rm = TRUE)
      # theRandomWRMSD[i, 4, k] <- quantile(theS[i, 4, ],(1-theSig/2), na.rm = TRUE)
      
      T_ABS <- t.test(theS[i, 1, ], alternative = theAlt, conf.level = (1-theSig))
      T_RMSD <- t.test(theS[i, 2, ], alternative = theAlt, conf.level = (1-theSig))
      T_WABS <- t.test(theS[i, 3, ], alternative = theAlt, conf.level = (1-theSig))
      T_WRMSD <- t.test(theS[i, 4, ], alternative = theAlt, conf.level = (1-theSig))
      
      theRandomABS[i, 2, k] <- T_ABS$estimate[[1]]
      theRandomRMSD[i, 2, k] <- T_RMSD$estimate[[1]]
      theRandomWABS[i, 2, k] <- T_WABS$estimate[[1]]
      theRandomWRMSD[i, 2, k] <- T_WRMSD$estimate[[1]]

      theRandomABS[i, 3, k] <- T_ABS$conf.int[[1]]
      theRandomRMSD[i, 3, k] <- T_RMSD$conf.int[[1]]
      theRandomWABS[i, 3, k] <- T_WABS$conf.int[[1]]
      theRandomWRMSD[i, 3, k] <- T_WRMSD$conf.int[[1]]

      theRandomABS[i, 4, k] <- T_ABS$conf.int[[2]]
      theRandomRMSD[i, 4, k] <- T_RMSD$conf.int[[2]]
      theRandomWABS[i, 4, k] <- T_WABS$conf.int[[2]]
      theRandomWRMSD[i, 4, k] <- T_WRMSD$conf.int[[2]]
    }
  }

  theRandomABS[, 1, ] <- theABS
  theRandomRMSD[, 1, ] <- theRMSD
  theRandomWABS[, 1, ] <- theWABS
  theRandomWRMSD[, 1, ] <- theWRMSD

  for (k in 1:nFac) {
    if(theCompare == 'mean'){
      theRandomABS[theRandomABS[, 1, k] > theRandomABS[, 2, k], 5, k] <- 1
      theRandomRMSD[theRandomRMSD[, 1, k] > theRandomRMSD[, 2, k], 5, k] <- 1
      theRandomWABS[theRandomWABS[, 1, k] > theRandomWABS[, 2, k], 5, k] <- 1
      theRandomWRMSD[theRandomWRMSD[, 1, k] > theRandomWRMSD[, 2, k], 5, k] <- 1
    }else if(theCompare == 'between'){
      for (i in 1:nTestee){
        if(theRandomABS[i, 1, k] > theRandomABS[i, 3, k] && theRandomABS[i, 1, k] < theRandomABS[i, 4, k]){
          theRandomABS[i, 5, k] <- 1
        }

        if(theRandomRMSD[i, 1, k] > theRandomRMSD[i, 3, k] && theRandomRMSD[i, 1, k] < theRandomRMSD[i, 4, k]){
          theRandomRMSD[i, 5, k] <- 1
        }

        if(theRandomWABS[i, 1, k] > theRandomWABS[i, 3, k] && theRandomWABS[i, 1, k] < theRandomWABS[i, 4, k]){
          theRandomWABS[i, 5, k] <- 1
        }

        if(theRandomWRMSD[i, 1, k] > theRandomWRMSD[i, 3, k] && theRandomWRMSD[i, 1, k] < theRandomWRMSD[i, 4, k]){
          theRandomWRMSD[i, 5, k] <- 1
        }
      }
    }else if(theCompare == 'lower'){
      theRandomABS[theRandomABS[, 1, k] < theRandomABS[, 3, k], 5, k] <- 1
      theRandomRMSD[theRandomRMSD[, 1, k] < theRandomRMSD[, 3, k], 5, k] <- 1
      theRandomWABS[theRandomWABS[, 1, k] < theRandomWABS[, 3, k], 5, k] <- 1
      theRandomWRMSD[theRandomWRMSD[, 1, k] < theRandomWRMSD[, 3, k], 5, k] <- 1
    }else if(theCompare == 'upper'){
      theRandomABS[theRandomABS[, 1, k] > theRandomABS[, 4, k], 5, k] <- 1
      theRandomRMSD[theRandomRMSD[, 1, k] > theRandomRMSD[, 4, k], 5, k] <- 1
      theRandomWABS[theRandomWABS[, 1, k] > theRandomWABS[, 4, k], 5, k] <- 1
      theRandomWRMSD[theRandomWRMSD[, 1, k] > theRandomWRMSD[, 4, k], 5, k] <- 1
    }
    
  }

  colnames(theRandomABS) <- c('observe', 'mean', 'lower', 'upper', 'sign')
  colnames(theRandomRMSD) <- c('observe', 'mean', 'lower', 'upper', 'sign')
  colnames(theRandomWABS) <- c('observe', 'mean', 'lower', 'upper', 'sign')
  colnames(theRandomWRMSD) <- c('observe', 'mean', 'lower', 'upper', 'sign')

  theRandomS <- list()
  theRandomS[[1]] <- theRandomABS
  theRandomS[[2]] <- theRandomRMSD
  theRandomS[[3]] <- theRandomWABS
  theRandomS[[4]] <- theRandomWRMSD

  names(theRandomS) <- c('ABS', 'RMSD', 'WABS', 'WRMSD')

  return(theRandomS)
}
