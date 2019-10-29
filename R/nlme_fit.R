################################################################################
# Copyright (c) 2015, 2016, 2017 Genome Research Ltd. 
# Copyright (c) 2015, 2016, 2017 The Netherlands Cancer Institute (NKI)
#  
# Author: Howard Lightfoot <cancerrxgene@sanger.ac.uk> 
# Author: Dieudonne van der Meer
# Author: Daniel J Vis
# 
# This file is part of gdscIC50. 
# 
# gdscIC50 is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 3 of the License, or (at your option) any later 
# version. 
#  
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details. 
#  
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see <http://www.gnu.org/licenses/>. 
################################################################################
#' @import nlme
NULL

groupNlmeData <- function(nlme_data) {
  gDat <- nlme::groupedData(y ~ x | drug/CL, data = nlme_data, FUN = mean,
                      labels = list(x = "Concentration", y = "Viability"),
                      units = list(x = "uM/l", y = "percentage killed"))
  return(gDat)
}

#' Fit an nlme dose reponse model to GDSC data.
#' 
#' @param nlme_data a data frame as specified by \code{\link{prepNlmeData}}
#' @param isLargeData logical. Default \code{(isLargeData = T)}. Set to \code{F}
#'  if model fitting fails to converge, e.g., 
#'  \code{step halving factor reduced below minimum in PNLS step}
#'  
#'  With this \code{isLargeData = TRUE}, the covariance between the position and
#'  scale parameter on the cell line level are assumed to be correlated. This 
#'  further stabilizes the fit. In small bespoke screens this is set to FALSE as
#'  the model otherwise struggles to converge.
#' 
#' @return a fitted nlme model
#' 
#' @examples 
#' data("gdsc_example")
#' gdsc_example <- removeFailedDrugs(gdsc_example)
#' gdsc_example <- removeMissingDrugs(gdsc_example)
#' gdsc_example <- normalizeData(gdsc_example)
#' gdsc_example <- setConcsForNlme(gdsc_example)
#' nlme_data <- prepNlmeData(gdsc_example, "COSMIC_ID")
#' \dontrun{
#' nlme_model <- fitModelNlmeData(nlme_data, isLargeData = F)
#' }
#' 
#' @seealso \code{\link{prepNlmeData}}, \code{\link{getModelCoef}}, 
#'  \code{\link{calcNlmeStats}}.
#' 
#' @export
fitModelNlmeData <- function (nlme_data, isLargeData = T) {
  gDat <- groupNlmeData(nlme_data)
  nlme_model <- fitModel(gDat, bLargeScale = isLargeData)
  return(nlme_model)
}

#' fitModel
#' 
#' Used internally by fitModelNlmeData
#' 
#' @param gDat grouped data frame input
#' @param vStart starting parameter values for xmid and scal
#' @param bLargeScale equivalent to isLargeData in fitModelNlmeData
#' @param bSilent verbose output or not
#'
#' @return a fitted nlme model
#'
fitModel <- function(gDat, vStart = c(8.886464, 1.495953 ), bLargeScale=TRUE, bSilent=TRUE){
  # start with sanity checks; catch most common mistakes early
  # first, do we have all requires columns?
  if(length(which(colnames(gDat)=='CL'))==0){
    stop('gDat is required to contain the cell lines names in a column with name CL')  
  }
  if(length(which(colnames(gDat)=='x'))==0){
    stop('gDat is required to contain the X-from concentration step names in a column with name x')  
  }
  if(length(which(colnames(gDat)=='y'))==0){
    stop('gDat is required to contain the relative kill (1-viability) in a column with name y')  
  }
  if(length(which(colnames(gDat)=='drug'))==0){
    stop('gDat is required to contain the drug in a column with name drug')  
  }
  if(length(which(colnames(gDat)=='maxc'))==0){
    stop('gDat is required to contain the maxc (maximum concentration) in a column with name maxc')  
  }
  # coding of relative kill (1-viability) correct; mean of gDat$y at gDat$x==9 should be higher than at min(gDat$x)
  tmp.minx <- min(gDat$x)
  tmp.whichxmin <- which(gDat$x == tmp.minx)
  tmp.whichxmax <- which(gDat$x == 9)
  tmp.ymin <- mean(stats::na.omit(gDat$y[tmp.whichxmin]))
  tmp.ymax <- mean(stats::na.omit(gDat$y[tmp.whichxmax]))
  if(tmp.ymin > tmp.ymax){
    stop('Coding of relative viabilities seems incorrect; note that y is defined as 1-viability.')
  }
  if(bLargeScale){
    fmv5 <- nlme::nlme(y ~ gdscIC50::logist3(x, xmid, scal),
                       fixed= xmid+scal~1,
                       random=list(CL = nlme::pdSymm(xmid+scal~1), 
                                   drug = nlme::pdDiag(xmid~1)),
                       data = gDat, start=vStart, method='REML')
    if(!bSilent){
      summary(fmv5)
    }
    return(fmv5)
  }else{
    fmv5 <- nlme::nlme(y ~ gdscIC50::logist3(x, xmid, scal),
                       fixed = xmid + scal ~  1, 
                       random = list(CL = nlme::pdDiag(xmid + scal ~ 1),
                                     drug = nlme::pdDiag(xmid ~ 1)), 
                       data = gDat, start = vStart, method = "REML",
                       control = nlme::nlmeControl(pnlsTol = 0.2,
                                             msVerbose = FALSE,
                                             tolerance=1e-4,
                                             returnObject = T)
                       )
    if(!bSilent){
      summary(fmv5)
    }
    return(fmv5)
  }
}

#' logist3
#' 
#' The 2 parameter logistic function used to model GDSC dose response data. This
#' model is described further in published in Vis, D.J. et al. Pharmacogenomics
#'  2016, 17(7):691-700)
#'
#' @param x concentration of drug treatment usually on transformed scale. 
#'   See setXFromConc().
#' @param xmid location parameter of logistic curve
#' @param scal scale parameter of logistic curve
#'
#' @return a function object of class selfStart
#' @export
logist3 <- stats::selfStart( ~ 1/(1 + exp(-(x - xmid)/scal)),
                             initial = function(mCall, LHS, data){  
                               xy <- stats::sortedXyData(mCall[["x"]], LHS, data)
                               if(nrow(xy) < 3) {
                                 stop("Too few distinct input values to fit a logistic")
                                 }
                               xmid <- stats::NLSstClosestX(xy, 0.5 ) 
                               scal <- stats::NLSstClosestX(xy, 0.75 ) - xmid
                               value <- c(xmid, scal)
                               names(value) <- mCall[c("xmid", "scal")]
                               value
                               },
                             parameters = c("xmid", "scal")
                             )


logistInit4 <- function(mCall, LHS, data){
    xy <- sortedXyData(mCall[["x"]], LHS, data)
    if(nrow(xy) < 3) {
        stop("Too few distinct input values to fit a logistic")
    }
    xmid <- NLSstClosestX(xy, 0.5 )
    scal <- NLSstClosestX(xy, 0.75 ) - xmid
    ANCHOR_VIAB <- 0.9
    value <- c(ANCHOR_VIAB,xmid, scal)
    names(value) <- mCall[c("xmid")]
    value
}

#' logist4
#' 
#' The 2 parameter logistic function used to model the combination expected curve. 
#' Takes the anchor viability into account
#'
#' @param x concentration of drug treatment usually on transformed scale. 
#'   See setXFromConc().
#' @param xmid location parameter of logistic curve
#' @param scal scale parameter of logistic curve
#' @param ANCHOR_VIAB the anchor viability
#'
#' @return a function object of class selfStart
#' @export
logist4 <- stats::selfStart( ~ (1-ANCHOR_VIAB)+(ANCHOR_VIAB/(1 + exp(-(x - xmid)/(scal)))), 
                             initial = logistInit4, parameters = c("xmid"))


logistInit5 <- function(mCall, LHS, data)
{
    xy <- sortedXyData(mCall[["x"]], LHS, data)
    if(nrow(xy) < 3) {
        stop("Too few distinct input values to fit a logistic")
    }
    xmid <- NLSstClosestX(xy, 0.5 )
    scal <- NLSstClosestX(xy, 0.75 ) - xmid
    # ANCHOR_VIAB <- 0.99
    value <- c(ANCHOR_VIAB,scal,xmid)
    names(value) <- mCall[c("xmid","scal")]
    value
}

#' logist5
#' 
#' The 2 parameter logistic function used to model combination observed curve. 
#'
#' @param x concentration of drug treatment usually on transformed scale. 
#'   See setXFromConc().
#' @param SYNERGY_XMID location parameter of logistic curve of the combination observed
#' @param LIBRARY_SACL scale parameter of logistic curve for the combination observed comes
#'   from the library alone.
#' @param ANCHOR_VIAB the anchor viability
#'
#' @return a function object of class selfStart
#' @export
logist5 <- stats::selfStart( ~ (1-ANCHOR_VIAB)+(ANCHOR_VIAB/(1 + exp(-(x - xmid)/(scal)))), 
                             initial = logistInit5, parameters = c("xmid","scal"))



#' getIntegral
#'
#' @param x concentration of drug treatment usually on transformed scale. 
#'   See setXFromConc()
#' @param xmid location parameter of logistic curve
#' @param scal scale parameter of logistic curve
#'
#' @return integral value
#' @export
getIntegral <- function(x,xmid,scal){
  a <- xmid
  b <- scal
  return(b*log(exp(a/b)+exp(x/b))-a)
}

