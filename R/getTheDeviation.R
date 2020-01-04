#' @title Calculate the Deviations
#'
#' @description Calculate deviations between observations and expectations.
#'
#' @param theData The prepared data matrix.
#'
#' @param theExpectPoint The expectation list.
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
#' @param tempPrint = ''
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
                            theAlt = 'less', theSig = 0.05, tempPrint = '', ifItem = FALSE,
                            theCompareItem = 'upper', theAltItem = 'less', theSigItem = 0.05){

  nTestee <- length(theData[, 1])
  nFac <- max(theDim)
  theMAD <- array(NA, dim = c(nTestee, nFac))
  theRMSD <- array(NA, dim = c(nTestee, nFac))
  theWMAD <- array(NA, dim = c(nTestee, nFac))
  theWRMSD <- array(NA, dim = c(nTestee, nFac))
  theRandomMAD <- array(NA, dim = c(nTestee, 5, nFac)) #1observe 2mean 3lower 4upper 5sign
  theRandomRMSD <- array(NA, dim = c(nTestee, 5, nFac))
  theRandomWMAD <- array(NA, dim = c(nTestee, 5, nFac))
  theRandomWRMSD <- array(NA, dim = c(nTestee, 5, nFac))
  theSign_item <- array(NA, dim = c(nTestee, length(theCategory)))
  theSign_item_WABS <- array(NA, dim = c(nTestee, length(theCategory)))
  for (k in 1:nFac) {
    # k <- 1
    theC <- theCategory[theDim==k]
    theDataK <- theData[, theDim==k]
    theWeight <- theExpectPoint$probability[, theDim==k, ]
    thePoint <- theExpectPoint$expectation[, theDim==k]
    theSign_item_K <- theSign_item[, theDim==k]
    theSign_item_K_WABS <- theSign_item_WABS[, theDim==k]

    for (i in 1:nTestee){
      # i <- 1
      theI <- length(theC[!is.na(thePoint[i,])])
      tempD <- 1:length(theDataK[i,])
      for (j in 1:length(theDataK[i,])) {
        # j <- 1
        if(!is.na(thePoint[i,j]) && !is.na(theDataK[i,j])){
          tempD[j] <- log(theWeight[i,j,(theDataK[i,j]+1)] / theWeight[i,j,(thePoint[i,j]+1)]) * (theDataK[i,j] - thePoint[i,j])
        }
      }

      theMAD[i, k] <- sum(abs(thePoint[i, !is.na(thePoint[i,])] - theDataK[i, !is.na(thePoint[i,])])) / theI
      theRMSD[i, k] <- sqrt(sum((thePoint[i, !is.na(thePoint[i,])] - theDataK[i, !is.na(thePoint[i,])])^2) / theI)
      theWMAD[i, k] <- sum(abs(tempD), na.rm = TRUE) / theI
      theWRMSD[i, k] <- sqrt(sum((tempD)^2, na.rm = TRUE) / theI)

      # theWMAD[i, k] <- sum(abs(thePoint[i, !is.na(thePoint[i,])] - theDataK[i, !is.na(thePoint[i,])]) / theC[!is.na(thePoint[i,])]) / theI
      # theWRMSD[i, k] <- sqrt(sum(((thePoint[i, !is.na(thePoint[i,])] - theDataK[i, !is.na(thePoint[i,])]) / theC[!is.na(thePoint[i,])])^2) / theI)
    }

    thePointK <- list()
    thePointK[[1]] <- theWeight
    thePointK[[2]] <- thePoint
    names(thePointK) <- c('probability','expectation')
    tempPrintD <- paste0(tempPrint, 'factor', k, '/', nFac, ':')
    # print(tempPrintD)
    theSS <- getTheBoot(theDataK, thePointK, theC, RandomReplaceRatio = RandomReplaceRatio, maxBoot = maxBoot, tempPrint = tempPrintD, ifItem = ifItem)
    theS <- unlist(theSS[[1]])
    theS_item <- unlist(theSS[[2]])
    theS_item_WABS <- unlist(theSS[[3]])

    if(ifItem){
      for (i in 1:nTestee){
        for (j in 1:length(theC)){
          # i<-1
          # j<-1
          if (!is.na(theDataK[i,j])){
            d_item <- abs(thePoint[i,j] - theDataK[i,j])
            t_item <- t.test(abs(theS_item[i,j,]),  alternative = theAltItem, conf.level = (1-theSigItem))

            if(!is.na(t_item$statistic[[1]])){
              if(theCompareItem == 'mean'){
                if(d_item > t_item$estimate[[1]]){
                  theSign_item_K[i,j] <- 1
                }
              }else if(theCompareItem == 'between'){
                if(d_item < t_item$conf.int[[1]] || d_item > t_item$conf.int[[2]]){
                  theSign_item_K[i,j] <- 1
                }
              }else if(theCompareItem == 'lower'){
                if(d_item < t_item$conf.int[[1]]){
                  theSign_item_K[i,j] <- 1
                }
              }else if(theCompareItem == 'upper'){
                if(d_item > t_item$conf.int[[2]]){
                  theSign_item_K[i,j] <- 1
                }
              }
            }

            d_item <- abs(log(theWeight[i,j,(theDataK[i,j]+1)] / theWeight[i,j,(thePoint[i,j]+1)]) * (theDataK[i,j] - thePoint[i,j]))
            t_item <- t.test(abs(theS_item_WABS[i,j,]),  alternative = theAltItem, conf.level = (1-theSigItem))

            if(!is.na(t_item$statistic[[1]])){
              if(theCompareItem == 'mean'){
                if(d_item > t_item$estimate[[1]]){
                  theSign_item_K_WABS[i,j] <- 1
                }
              }else if(theCompareItem == 'between'){
                if(d_item < t_item$conf.int[[1]] || d_item > t_item$conf.int[[2]]){
                  theSign_item_K_WABS[i,j] <- 1
                }
              }else if(theCompareItem == 'lower'){
                if(d_item < t_item$conf.int[[1]]){
                  theSign_item_K_WABS[i,j] <- 1
                }
              }else if(theCompareItem == 'upper'){
                if(d_item > t_item$conf.int[[2]]){
                  theSign_item_K_WABS[i,j] <- 1
                }
              }
            }

          }
        }
      }

      theSign_item[, theDim==k] <- theSign_item_K
      theSign_item_WABS[, theDim==k] <- theSign_item_K_WABS
    }

    for (i in 1:nTestee){
      # theRandomMAD[i, 2, k] <- mean(theS[i, 1, ], na.rm = TRUE)
      # theRandomRMSD[i, 2, k] <- mean(theS[i, 2, ], na.rm = TRUE)
      # theRandomWMAD[i, 2, k] <- mean(theS[i, 3, ], na.rm = TRUE)
      # theRandomWRMSD[i, 2, k] <- mean(theS[i, 4, ], na.rm = TRUE)
      #
      # theRandomMAD[i, 3, k] <- quantile(theS[i, 1, ],(theSig/2), na.rm = TRUE)
      # theRandomRMSD[i, 3, k] <- quantile(theS[i, 2, ],(theSig/2), na.rm = TRUE)
      # theRandomWMAD[i, 3, k] <- quantile(theS[i, 3, ],(theSig/2), na.rm = TRUE)
      # theRandomWRMSD[i, 3, k] <- quantile(theS[i, 4, ],(theSig/2), na.rm = TRUE)
      #
      # theRandomMAD[i, 4, k] <- quantile(theS[i, 1, ],(1-theSig/2), na.rm = TRUE)
      # theRandomRMSD[i, 4, k] <- quantile(theS[i, 2, ],(1-theSig/2), na.rm = TRUE)
      # theRandomWMAD[i, 4, k] <- quantile(theS[i, 3, ],(1-theSig/2), na.rm = TRUE)
      # theRandomWRMSD[i, 4, k] <- quantile(theS[i, 4, ],(1-theSig/2), na.rm = TRUE)

      T_MAD <- t.test(theS[i, 1, ], alternative = theAlt, conf.level = (1-theSig))
      T_RMSD <- t.test(theS[i, 2, ], alternative = theAlt, conf.level = (1-theSig))
      T_WMAD <- t.test(theS[i, 3, ], alternative = theAlt, conf.level = (1-theSig))
      T_WRMSD <- t.test(theS[i, 4, ], alternative = theAlt, conf.level = (1-theSig))

      theRandomMAD[i, 2, k] <- T_MAD$estimate[[1]]
      theRandomRMSD[i, 2, k] <- T_RMSD$estimate[[1]]
      theRandomWMAD[i, 2, k] <- T_WMAD$estimate[[1]]
      theRandomWRMSD[i, 2, k] <- T_WRMSD$estimate[[1]]

      theRandomMAD[i, 3, k] <- T_MAD$conf.int[[1]]
      theRandomRMSD[i, 3, k] <- T_RMSD$conf.int[[1]]
      theRandomWMAD[i, 3, k] <- T_WMAD$conf.int[[1]]
      theRandomWRMSD[i, 3, k] <- T_WRMSD$conf.int[[1]]

      theRandomMAD[i, 4, k] <- T_MAD$conf.int[[2]]
      theRandomRMSD[i, 4, k] <- T_RMSD$conf.int[[2]]
      theRandomWMAD[i, 4, k] <- T_WMAD$conf.int[[2]]
      theRandomWRMSD[i, 4, k] <- T_WRMSD$conf.int[[2]]
    }
  }

  theRandomMAD[, 1, ] <- theMAD
  theRandomRMSD[, 1, ] <- theRMSD
  theRandomWMAD[, 1, ] <- theWMAD
  theRandomWRMSD[, 1, ] <- theWRMSD

  for (k in 1:nFac) {
    if(theCompare == 'mean'){
      theRandomMAD[theRandomMAD[, 1, k] > theRandomMAD[, 2, k], 5, k] <- 1
      theRandomRMSD[theRandomRMSD[, 1, k] > theRandomRMSD[, 2, k], 5, k] <- 1
      theRandomWMAD[theRandomWMAD[, 1, k] > theRandomWMAD[, 2, k], 5, k] <- 1
      theRandomWRMSD[theRandomWRMSD[, 1, k] > theRandomWRMSD[, 2, k], 5, k] <- 1
    }else if(theCompare == 'between'){
      for (i in 1:nTestee){
        if(theRandomMAD[i, 1, k] < theRandomMAD[i, 3, k] || theRandomMAD[i, 1, k] > theRandomMAD[i, 4, k]){
          theRandomMAD[i, 5, k] <- 1
        }

        if(theRandomRMSD[i, 1, k] < theRandomRMSD[i, 3, k] || theRandomRMSD[i, 1, k] > theRandomRMSD[i, 4, k]){
          theRandomRMSD[i, 5, k] <- 1
        }

        if(theRandomWMAD[i, 1, k] < theRandomWMAD[i, 3, k] || theRandomWMAD[i, 1, k] > theRandomWMAD[i, 4, k]){
          theRandomWMAD[i, 5, k] <- 1
        }

        if(theRandomWRMSD[i, 1, k] < theRandomWRMSD[i, 3, k] || theRandomWRMSD[i, 1, k] > theRandomWRMSD[i, 4, k]){
          theRandomWRMSD[i, 5, k] <- 1
        }
      }
    }else if(theCompare == 'lower'){
      theRandomMAD[theRandomMAD[, 1, k] < theRandomMAD[, 3, k], 5, k] <- 1
      theRandomRMSD[theRandomRMSD[, 1, k] < theRandomRMSD[, 3, k], 5, k] <- 1
      theRandomWMAD[theRandomWMAD[, 1, k] < theRandomWMAD[, 3, k], 5, k] <- 1
      theRandomWRMSD[theRandomWRMSD[, 1, k] < theRandomWRMSD[, 3, k], 5, k] <- 1
    }else if(theCompare == 'upper'){
      theRandomMAD[theRandomMAD[, 1, k] > theRandomMAD[, 4, k], 5, k] <- 1
      theRandomRMSD[theRandomRMSD[, 1, k] > theRandomRMSD[, 4, k], 5, k] <- 1
      theRandomWMAD[theRandomWMAD[, 1, k] > theRandomWMAD[, 4, k], 5, k] <- 1
      theRandomWRMSD[theRandomWRMSD[, 1, k] > theRandomWRMSD[, 4, k], 5, k] <- 1
    }

  }

  colnames(theRandomMAD) <- c('observe', 'mean', 'lower', 'upper', 'sign')
  colnames(theRandomRMSD) <- c('observe', 'mean', 'lower', 'upper', 'sign')
  colnames(theRandomWMAD) <- c('observe', 'mean', 'lower', 'upper', 'sign')
  colnames(theRandomWRMSD) <- c('observe', 'mean', 'lower', 'upper', 'sign')

  theRandomS <- list()
  theRandomS[[1]] <- theRandomMAD
  theRandomS[[2]] <- theRandomRMSD
  theRandomS[[3]] <- theRandomWMAD
  theRandomS[[4]] <- theRandomWRMSD
  theRandomS[[5]] <- theSign_item
  theRandomS[[6]] <- theSign_item_WABS
  names(theRandomS) <- c('MAD', 'RMSD', 'WMAD', 'WRMSD', 'sign_itemABS', 'sign_itemWABS')

  return(theRandomS)
}
