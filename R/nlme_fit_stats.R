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

#' @import ggplot2
NULL

#' Get the statistics from a fitted nlme dose reponse model.
#' 
#' @param nlme_model the results of \code{\link{fitModelNlmeData}}.
#' @param nlme_data the input data frame to the model fit.
#' 
#' @return a data frame of fitted nlme model statistics - xmid (IC50)
#'  and scal (slope).
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
#' nlme_stats <- calcNlmeStats(nlme_model, nlme_data)
#' }
#' 
#' @seealso \code{\link{fitModelNlmeData}}, \code{\link{calcNlmeStats}}.
#' 
#' @export
getModelCoef <- function(nlme_model, nlme_data) {
  drug_specifier <- nlme_data %>% select_(~drug_spec) %>% distinct()
  stopifnot(nrow(drug_specifier) == 1)
  drug_specifier <- strsplit(drug_specifier$drug_spec, "\\+")[[1]]
  model_coef <- stats::coef(nlme_model)
  model_coef <- model_coef %>% mutate_(combo = ~rownames(.)) %>%
    tidyr::separate_("combo", into = c("CL", "drug"), sep = "/", convert = T)
  model_coef <- model_coef %>% mutate_(drug = ~ as.character(drug))
  model_coef <- left_join(model_coef, nlme_data, by = c("CL", "drug") )
  return(model_coef)
}

calcIC50 <- function(model_coef) {
  model_coef <- model_coef %>%
    mutate_(
      IC50 = ~log(getConcFromX(xmid, maxc)))
  return(model_coef)
}

calcAuc <- function(model_coef) {
  model_coef <- model_coef %>% 
    group_by_(~CL, ~drug) %>% 
    mutate_(xmin = ~min(x)) %>%
    mutate_(xmax = ~max(x)) %>%
    mutate_(auc = 
              ~1 - (getIntegral(xmax, xmid, scal) - 
                      getIntegral(xmin, xmid, scal)) / (xmax - xmin))
  model_coef <- model_coef %>% 
    ungroup() %>%
    select_(~-xmax, ~-xmin)
  return(model_coef)
}

calcAucTrap <- function(model_stats) {
  model_stats <- model_stats %>% 
    group_by_(~CL, ~drug) %>% arrange_(~x) %>% 
    mutate_(area = ~ (max(x) - min(x)) * max(c(1,yhat))) %>% 
    mutate_(AUCtrap = ~((caTools::trapz(x, 1 - yhat)) / area))
  model_stats <- model_stats %>% ungroup() %>% select_(~-area)
  return(model_stats)
}

calcNlmeFit <- function(model_coef){
  model_coef <- model_coef %>%
    group_by_(~CL, ~drug) %>%
    mutate_(yhat = ~ logist3(x, xmid, scal))
  attributes(model_coef$yhat) <- NULL
  model_coef <- model_coef %>% mutate_(yres = ~y - yhat)
  model_coef <- model_coef %>% mutate_(RMSE = ~sqrt(mean(yres ^ 2)))
  model_coef <- model_coef %>% ungroup()
  return(model_coef)
}

#' Calculate all statistics from dose reponse fit
#' 
#' @param nlme_model the results of \code{\link{fitModelNlmeData}}.
#' @param nlme_data the nlme data frame used for fitting \code{\link{prepNlmeData}}.
#' 
#' @return a data frame inluding IC50, AUC and RMSE.
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
#' nlme_stats <- calcNlmeStats(nlme_model, nlme_data)
#' }
#' 
#' @seealso \code{\link{fitModelNlmeData}}, \code{\link{getModelCoef}}, 
#'  \code{\link{getIC50Matrix}}
#' 
#' @export
calcNlmeStats <- function (nlme_model, nlme_data) {
#   gDat <- groupNlmeData(nlme_data)
  model_coef <- getModelCoef(nlme_model = nlme_model, nlme_data = nlme_data)
  model_stats <- calcNlmeFit(model_coef)
  model_stats <- model_stats %>%
    arrange_(~desc(x)) %>%
    mutate_(x_micromol = ~getConcFromX(x, maxc))
  model_stats <- calcIC50(model_stats)
  model_stats <- calcAuc(model_stats)
  model_stats <- calcAucTrap(model_stats)
  return(model_stats)
}

#' Extract a matrix of IC50 values from a model statistics data frame.
#' 
#' @param model_stats the results of \code{\link{calcNlmeStats}}.
#' @param drug_identifier a character vector indicating which column of 
#'   nlme_stats to use as the "drug" in the IC50matrix. Normally this is the
#'    \code{drug} column as used in the model fit formula but you may want to 
#'    change to another identifier for user friendliness e.g. DRUG_ID or 
#'    DRUG_NAME. The function checks whether there is a one to one mapping 
#'    between the chosen identifier and the \code{drug} column before proceeding.
#'    
#'    The cell line column is named according to the \code{cl_spec} contained in 
#'    the \code{model_stats}.
#' @param re_name logical. Rename the drug columns to GDSCtools standard 
#'   (Drug_..._IC50). Default is TRUE.
#' @param measure character. Default to "IC50" but can also be set to "auc".
#'    
#' @return a matrix of IC50s.
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
#' nlme_stats <- calcNlmeStats(nlme_model, nlme_data)
#' ic50s <- getIC50Matrix(nlme_stats)
#' aucs <- getIC50Matrix(nlme_stats, measure = "auc")
#' }
#' 
#' @seealso \code{\link{getModelCoef}}
#' 
#' @export
getIC50Matrix <- function(model_stats, drug_identifier = "drug",
                          re_name = T, measure = "IC50") {
  
  # Make sure there is a 1-to-1 correspondence between chosen drug identifier and the drug column
  possible_duplicates <- model_stats %>% 
    select_(drug_identifier, ~drug) %>% 
    distinct() %>% 
    count_(drug_identifier) %>% 
    filter(n > 1) %>% 
    nrow()
  stopifnot(possible_duplicates == 0)
  
  
  
  cl_spec <- unique(model_stats$CL_SPEC)
  
  if(re_name){
    IC50_matrix <- model_stats %>%  
      select_(~CL, drug_identifier, measure) %>% 
      distinct() %>%
      mutate_(drug = lazyeval::interp(~ paste("Drug", foo, measure, sep = "_"),
                                      foo = as.name(drug_identifier)))
    
    if (identical(drug_identifier, "drug")){
      IC50_matrix <- IC50_matrix %>% tidyr::spread_("drug", measure)
    } 
    else {
      IC50_matrix <- IC50_matrix %>% 
        select_(.dots = list(paste("-", drug_identifier))) %>%
        tidyr::spread_("drug", measure)
    }
  }
  else{
    IC50_matrix <- model_stats %>%  
      select_(~CL, drug_identifier, measure) %>% 
      distinct() %>%
      tidyr::spread_(drug_identifier, measure)
  }
  
  IC50_matrix <- IC50_matrix %>% rename_(.dots = stats::setNames('CL', cl_spec))
  
  return(IC50_matrix)
}


l3_model <- function(x, xmid, scal){
  #   yhat <- 1/(1 + exp(1) ^ ((x - xmid) / scal))
  yhat <- 1 - logist3(x, xmid, scal)
  return(yhat)
}


l3_model2 <- function(lx, maxc, xmid, scal){
  x <- getXfromConc(exp(lx), maxc)
#     yhat <- 1/(1 + exp(1) ^ ((x - xmid) / scal))
  yhat <- 1 - logist3(x, xmid, scal)
  return(yhat)
}

getConcFromX <- function(x, maxc) {
  xc <- maxc * 2 ^ (x - 9)
  return(xc)
}

getXfromConc <- function(xc, maxc) {
  x <- (log(xc / maxc)/log(2))+ 9
  return(x)
}

# Is the following correct? Are we using y  or 1-y???
# IC50 is fine because symmetric
# getX <- function(y, xmid, scal){
#   x <- xmid - scal*(log((1 - y) / y))
#   return(x)
#   }
getX <- function(y, xmid, scal){
  x <- xmid - scal*(log((1 - y) / y))
  return(x)
}

#' plotResponse
#'
#' Plot dose reponse curve.
#' @param model_stats dataframe of fitted values as produced by \code{calcNlmeStats}
#' @param cell_line as identified from the \code{CL} column of the model stats data frame.
#' @param drug_identifier character string as identified by the \code{drug} column of the 
#' model stats data frame.
#'
#' @return ggplot of dose response curve.
#' @export
#'
#' @examples
#' \dontrun{
#' nlme_model <- fitModelNlmeData(nlme_data, isLargeData = F)
#' nlme_stats <- calcNlmeStats(nlme_model, nlme_data)
#' plotResponse(model_stats = nlme_stats, 
#'              cell_line = 1503364,
#'              drug_identifier = "1032")
#' }

plotResponse <- function(model_stats, cell_line, drug_identifier) {
  plot_data <- model_stats %>%
    filter_(~CL == cell_line, ~drug == drug_identifier)
  stopifnot(nrow(plot_data) > 0)
  
  IC50 <- unique(plot_data$IC50)
  stopifnot(length(IC50) == 1)
  
  auc <- unique(plot_data$auc)
  stopifnot(length(auc) == 1)
  
  rmse <- unique(plot_data$RMSE)
  stopifnot(length(rmse) == 1)
  
  drug_id <- unique(plot_data$DRUG_ID_lib)
  stopifnot(length(drug_id) == 1)
  
  cell_line_name <- unique(plot_data$CELL_LINE_NAME)
  cell_line_name <- ifelse(is.null(cell_line_name),  "", cell_line_name)
  stopifnot(length(cell_line_name) == 1)
  
  max_conc <- unique(plot_data$maxc)
  stopifnot(length(max_conc) == 1)

  plot_data <- plot_data %>%
    mutate_(lx = ~log(getConcFromX(x, maxc)),
            lxmid = ~log(getConcFromX(xmid, maxc))
            )
  
  plot_xmid <- plot_data %>% select_(~xmid) %>% distinct()
  plot_scal <- plot_data %>% select_(~scal) %>% distinct()
  plot_maxc <- plot_data %>% select_(~maxc) %>% distinct()
  
  plot_low_x <- 1 - plot_scal$scal * log((1 - 1e-3) / 1e-3) + plot_xmid$xmid
  plot_low_x <- log(getConcFromX(plot_low_x, +plot_maxc$maxc))
  plot_low_x <- min(c(plot_data$lx, plot_low_x))
  
  plot_high_x <- 1 - plot_scal$scal * log(1e-3 / (1 - 1e-3)) + plot_xmid$xmid
  plot_high_x <- log(getConcFromX(plot_high_x, plot_maxc$maxc))
  plot_high_x <- max(c(plot_data$lx, plot_high_x))
  
  p <- ggplot(plot_data) + aes_(x = ~lx, y = ~ 1 - yhat) + geom_point(shape = 3) +
  scale_x_continuous(limits = c(plot_low_x, plot_high_x))

  p <- p + stat_function(aes_(x = ~lx),
                         fun = l3_model2,
                         args = list(maxc = plot_maxc$maxc,
                                     xmid = plot_xmid$xmid,
                                     scal = plot_scal$scal))
  p <- p + annotate("rect",
                     xmin = min(plot_data$lx),
                     xmax = max(plot_data$lx),
                     ymin = 0,
                     ymax =  max(c(1, 1 - plot_data$yhat)), 
                     alpha = 0.2
  )
  
  p <- p + aes_(x = ~lx, y = ~ 1 - yhat) + geom_point(shape = 4)

  p <- p + geom_point(aes_(x = ~lx, y = ~1 - y), shape = 1)
  p <- p + geom_point(aes_(x = ~lxmid, y = 0.5), colour = "red") + theme(legend.position = "none")
  p <- p + geom_point(aes_(x = ~lxmid, y = 0.5), colour = "green", shape = 5, size = 3) + theme(legend.position = "none")
  p <- p + geom_linerange(aes_(x = ~IC50, ymax = 0.5, ymin = 0, colour = "red"), linetype = "dashed")
  p <- p + annotate("label",
                    x = IC50 + 1,
                    y = 0.5,
                    hjust = "left",
                    label = sprintf("IC50==%.3f~log[e]~mu*M", IC50), parse = T)
  p <- p + annotate("label",
                    x = min(plot_data$lx) + 0.1,
                    y = 0.25,
                    hjust = "left",
                    label = sprintf("auc = %.3f", auc))
  p  <- p + ylab("Response: normalized intensity") + 
    theme(axis.title.y = element_text(size=14))
  p <- p + xlab(expression(Dose/log[e]~mu*M)) + 
    theme(axis.title.x = element_text(size=14))
  
  p <- p + annotate("label", 
                     # plot label at 2/3 x width
                     x = plot_low_x + ((plot_high_x - plot_low_x) * 2 / 3), 
                     y =  max(c(1, 1 - plot_data$yhat)),
                     hjust = "left", 
                     vjust = "inward", 
                     # label = sprintf("rmse = %.3f", rmse)
                     label = sprintf("Cell-line: %s\nDrug id: %d\nMax dose = %.3f uM\nrmse = %.3f",
                                    cell_line_name, drug_id, max_conc, rmse))
  # label = paste(sprintf("Cell-line:~%s", cell_line_name), "\n", sprintf("Drug~id:~%d", drug_id)), parse = T)

  p <- p + ggtitle(paste("Dose response: cell line ", cell_line, "; drug ", drug_identifier, "."))
  
  return(p)
}
