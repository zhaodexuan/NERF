#' @title Estimate the GGUM
#'
#' @description Get the estimation of the GGUM.
#'
#' @param theDataK The data matrix of factor K.
#'
#' @param theC The category of factor K.
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
#' @return The function returns an array.
#'
#' @author zdx, \email{zhaodexuan@aliyun.com}
#'
#' @examples
#' \dontrun{
#' NA
#' }
#'
#' @importFrom GGUM GGUM
#'
#' @importFrom GGUM Theta.EAP
#'
#' @importFrom GGUM probs.GGUM
#'
#' @importFrom mirt mirt
#'
#' @importFrom mirt fscores
#'
#' @importFrom mirt coef
#'

getTheGGUM <- function(theDataK, theC, pack ='GGUM',
                       SE = FALSE, precision = 4, N_nodes = 30, max_outer = 60, max_inner = 60, tol = 0.001,
                       method = 'EM'){
  thePoint <- array(NA, dim = c(length(theDataK[, 1]), length(theDataK[1, ])))
  if(pack=='GGUM'){
    theFit <- GGUM(theDataK, theC, SE = SE, precision = precision, N.nodes = N_nodes, max.outer = max_outer, max.inner = max_inner, tol = tol)
    theTh <- Theta.EAP(theFit)
    theProb <- probs.GGUM(theFit$alpha, theFit$delta, theFit$taus, theTh[, 2], theC)

    for (i in 1:length(theDataK[, 1])) {
      for (j in 1:length(theDataK[1, ])) {
        if(!sum(is.na(theProb[i,j,1:(theC[j]+1)]))==(theC[j]+1)){
          thePoint[i, j] <- which.max(theProb[i, j, ])
          thePoint[i, j] <- thePoint[i, j] - 1
        }
      }
    }
  }else if(pack=='mirt'){
    theNames <- array()
    for (n in 1:length(theC)) {
      theNames[n] <- paste0('item', n)
    }

    colnames(theDataK) <- theNames

    theMirt <- mirt(theDataK,1,itemtype = 'ggum',SE = SE,method = method)

    theMirtTh <- fscores(theMirt,method="EAPsum",full.scores = TRUE)

    theMirtCoef <- coef(theMirt)
    theAlpha <- matrix(NA,length(theC),1)
    theBeta <- matrix(NA,length(theC),1)
    theTaus <- matrix(NA,length(theC),(max(theC)+1))
    for (n in 1:length(theC)) {
      tempCoefN <- theMirtCoef[[n]]
      theAlpha[n] <- tempCoefN[1]
      theBeta[n] <- tempCoefN[2]
      theTaus[n,2:(theC[n]+1)] <- tempCoefN[3:(2+theC[n])]
      theTaus[n,1] <- 0
    }

    theProb <- array(NA,dim = c(length(theDataK[,1]),length(theDataK[1,]),(max(theC)+1)))
    for (i in 1:length(theDataK[,1])) {
      for (j in 1:length(theDataK[1,])) {
        for (x in 1:(theC[j]+1)) {
          z <- x-1 # 0:z = 1:x
          fenZi1_1 <- z*sqrt(theAlpha[j]^2 * (theMirtTh[i]-theBeta[j])^2)
          fenZi1_2 <- sum(theAlpha[j] * theTaus[j,1:x])
          fenZi1 <- exp(fenZi1_1 + fenZi1_2)

          M <- 2 * theC[j] + 1
          fenZi2_1 <- (M-z)*sqrt(theAlpha[j]^2 * (theMirtTh[i]-theBeta[j])^2)
          fenZi2_2 <- sum(theAlpha[j] * theTaus[j,1:x])
          fenZi2 <- exp(fenZi2_1 + fenZi2_2)

          fenMu <- 0
          for (y in 1:(theC[j]+1)) {
            w <- y-1 # 0:w = 1:y
            fenMu1_1 <- w*sqrt(theAlpha[j]^2 * (theMirtTh[i]-theBeta[j])^2)
            fenMu1_2 <- sum(theAlpha[j] * theTaus[j,1:y])
            fenMu1 <- exp(fenMu1_1 + fenMu1_2)

            fenMu2_1 <- (M-w)*sqrt(theAlpha[j]^2 * (theMirtTh[i]-theBeta[j])^2)
            fenMu2_2 <- sum(theAlpha[j] * theTaus[j,1:y])
            fenMu2 <- exp(fenMu2_1 + fenMu2_2)

            fenMu <- fenMu + t(fenMu1 +fenMu2)
          }

          theProb[i,j,x] <- (fenZi1 + fenZi2) / fenMu
        }
      }
    }

    for (i in 1:length(theDataK[, 1])) {
      for (j in 1:length(theDataK[1, ])) {
        if(!sum(is.na(theProb[i,j,1:(theC[j]+1)]))==(theC[j]+1)){
          thePoint[i, j] <- which.max(theProb[i, j, ])
          thePoint[i, j] <- thePoint[i, j] - 1
        }
      }
    }
  }

  return(thePoint)
}
