################################################################################
# Copyright (c) 2015, 2016, 2017, 2018, 2019, 2020 Genome Research Ltd. 
# Copyright (c) 2015, 2016, 2017, 2018, 2019, 2020 The Netherlands Cancer Institute (NKI)
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
#' @import dplyr
NULL

#' Removes well positions where the drug is now considered unsuitable for
#'  screening from GDSC raw data, i.e., TAG = 'FAIL'
#' 
#' \code{removeFailedDrugs} removes rows from GDSC raw data where the
#'  \code{TAG == 'FAIL'}.
#' 
#' @param myDat a GDSC raw data data frame.
#' 
#' @seealso \code{\link{removeMissingDrugs}}, \code{\link{normalizeData}},
#'  \code{\link{setConcsForNlme}}, \code{\link{prepNlmeData}}
#'
#' @examples
#' data("gdsc_example")
#' gdsc_example <- removeFailedDrugs(gdsc_example)
#' gdsc_example <- removeMissingDrugs(gdsc_example)
#' gdsc_example <- normalizeData(gdsc_example)
#' gdsc_example <- setConcsForNlme(gdsc_example)
#' nlme_data <- prepNlmeData(gdsc_example, "COSMIC_ID")
#'  
#' @export
removeFailedDrugs <- function(myDat){
  # Remove the entire position because the failed drug might be part of a combination.
  failed_positions <- myDat %>% 
    filter(TAG == 'FAIL') %>%
    select(SCAN_ID, POSITION)
  myDat <-  anti_join(myDat, failed_positions, by = c("SCAN_ID", "POSITION"))
  return(myDat)
}

#' Removes library drugs listed with a drug id of NA from GDSC raw data. 
#' 
#' In contrast to drugs with 'FAIL' tags, the drugset contained NA for 
#'  the DRUG_ID, i.e., no drug id was assigned for that tag in the experimental
#'  design.
#' 
#' \code{removeMissingDrugs} removes rows from GDSC raw data where the
#'  \code{DRUG_ID} is NA.
#' 
#' @param myDat a GDSC raw data data frame.
#' 
#' @seealso  \code{\link{removeFailedDrugs}},  \code{\link{normalizeData}},
#'   \code{\link{setConcsForNlme}},  \code{\link{prepNlmeData}}
#'  
#' @examples
#' data("gdsc_example")
#' gdsc_example <- removeFailedDrugs(gdsc_example)
#' gdsc_example <- removeMissingDrugs(gdsc_example)
#' gdsc_example <- normalizeData(gdsc_example)
#' gdsc_example <- setConcsForNlme(gdsc_example)
#' nlme_data <- prepNlmeData(gdsc_example, "COSMIC_ID")
#' 
#' @export
removeMissingDrugs <- function(myDat){
  na_libs <- myDat %>%
    filter(grepl("^(L|R|A)\\d+", TAG)) %>% 
    filter(is.na(DRUG_ID))
  myDat <- anti_join(myDat, na_libs, by = c("SCAN_ID", "POSITION"))
  return(myDat)
}

#' Remove reference drugs from screen data
#'
#' @param myDat a GDSC raw data data frame.
#'
#' @return a GDSC raw data data frame with TAG s, such as R1-D1-S, removed.
#' @export
removeReferenceDrugs <- function(myDat){
  ref_libs <- myDat %>%
    filter(grepl("^(R)\\d+", TAG))
  myDat <- anti_join(myDat, ref_libs, by = c("SCAN_ID", "POSITION"))
  return(myDat)
}


#' findMultiLibs
#' 
#' Print out positions of multi-library treatments e.g. L1-Dx-S + L2-Dx-S
#' 
#' In some later GDSC combination screens single treatments were sometimes
#' composed of multiple library drugs. These cannot be processed as single 
#' drug treatments for nlme dose response fitting.  These positions can be 
#' filtered out using \code{\link{removeMultiLibs}}
#' from the raw screen data.
#'
#' @param myDat data frame of GDSC screen data 
#'
#' @return dataframe of replicate positions
#' @export
#'
#' @examples
#' \dontrun{
#' findMultiLibs(screen_data)
#' }
findMultiLibs <- function(myDat){
  repl_pos <- myDat %>% 
    filter(grepl("^L\\d+-D\\d+-S$", TAG)) %>% 
    distinct(DRUGSET_ID, POSITION, TAG, DRUG_ID) %>% 
    group_by(DRUGSET_ID, POSITION) %>% 
    summarise(ntags = n()) %>% 
    filter(ntags > 1)
  
  print("Multiple libraries in -S well. ")
  repl_pos %>% 
    purrr::pmap(function(DRUGSET_ID, POSITION, TAG, ...) print(
      paste0("Drugset: ", DRUGSET_ID, ", position: ", POSITION)
      )
    )
  return(repl_pos)
}

#' Removes -S treatments with more than one library drug 
#' 
#' In some later GDSC combination screens single treatments were sometimes
#' composed of multiple library drugs. These cannot be processed as single 
#' drug treatments for nlme dose response fitting and need to be filtered out
#' from the raw screen data.
#' 
#' \code{removeMultiLibs} removes rows from GDSC raw data where the
#'  \code{DRUG_ID} is NA.
#' 
#' @param myDat a GDSC raw data data frame.
#' 
#' @seealso  \code{\link{removeFailedDrugs}},  \code{\link{normalizeData}},
#'   \code{\link{setConcsForNlme}},  \code{\link{prepNlmeData}}
#'  
#' @examples
#' \dontrun{
#' screen_data <- screen_data %>% removeMissingDrugs() %>% removeFailedDrugs()
#' screen_data <- removeMultiLibs(screen_data)
#' norm_data <- normalizeData(screen_data)
#' nlme_data <- prepNlmeData(norm_data, "COSMIC_ID")
#' }
#' @export
removeMultiLibs <- function(myDat){
  repl_pos <- findMultiLibs(myDat)
  print("\nRemoved multilibs fom screening data\n")
  myDat <- myDat %>% 
    anti_join(repl_pos %>% select(-ntags),
              by = c("DRUGSET_ID", "POSITION"))
  return(myDat)
}

#' Normalizes GDSC raw data intensities with respect to controls.
#' 
#' \code{normalizeData} returns normalized intensities for the drug treated
#' wells - column \code{normalized_intensity}. 
#' 
#' Replaces \code{TAG} column with columns: \code{lib_drug}, \code{dose} , and 
#'  \code{treatment} (single or combination - S or C respectively).
#'   
#' If anchor drug treatments (anchor + library combinations) exist in the data 
#'   then additional columns are added: \code{DRUG_ID_anch}, \code{CONC_anch}
#'   and \code{anchor}. 
#'   
#' Combination treatments (C) are treated as if the library alone is to be 
#'   fitted.
#'
#' @param myDat a GDSC raw data data frame.
#' @param trim logical indicating whether to trim normalized values to the range
#' 0 to 1. default \code{(trim = T)}
#' @param neg_control The tag used to recognise a negative control well - the 
#'   upper end of the dynamic range.
#' @param pos_control The tag used to recognise a positive control well - the 
#'   lower end of the dynamic range.
#' 
#' @seealso  \code{\link{removeFailedDrugs}},  \code{\link{removeMissingDrugs}},
#'   \code{\link{setConcsForNlme}},  \code{\link{prepNlmeData}}
#' 
#' @examples
#' data("gdsc_example")
#' gdsc_example <- removeFailedDrugs(gdsc_example)
#' gdsc_example <- removeMissingDrugs(gdsc_example)
#' gdsc_example <- normalizeData(gdsc_example)
#' gdsc_example <- setConcsForNlme(gdsc_example)
#' nlme_data <- prepNlmeData(gdsc_example, "COSMIC_ID")
#'
#' @export
normalizeData <- function(myDat, trim = T, neg_control = 'NC-1',
                          pos_control = 'B'){
  
  nc1 <- calcTagMean(myDat, tag_name = neg_control, mean_col_name = "NC")
  pc1 <- calcTagMean(myDat, tag_name = pos_control, mean_col_name = "PC")
  
  # Take account of historic data with no MASTER_CELL_ID column - use ends_with
  normalized_data <- myDat %>% 
    filter(grepl("(A|L|R|PC)\\d+(-D\\d+)?-(S|C)", TAG)) %>%
    select(SCAN_ID, BARCODE,  RESEARCH_PROJECT, DATE_CREATED, DRUGSET_ID,
            CELL_LINE_NAME, ends_with("CELL_ID"),  COSMIC_ID, POSITION,
            TAG, DRUG_ID, CONC, INTENSITY) %>%
    mutate(lib_drug = sub("((L|R|PC)\\d+)(-D\\d+)?-(S|C)", "\\1", TAG),
            lib_drug = ifelse(grepl("^A.+", lib_drug), yes = NA, no = lib_drug),
            anchor = sub("(A\\d+)-(S|C)", "\\1", TAG),
            anchor = ifelse(grepl("^(L|R|PC).+", anchor), yes = NA, no = anchor),
            dose = sub("(A|L|R|PC)\\d+-?(D\\d+)?-(S|C)", "\\2", TAG),
            treatment = sub("((A|L|R|PC)\\d+)(-D\\d+)?-(S|C)", "\\4", TAG)
    ) %>%
    select(-TAG) 
  
  libraries <- normalized_data %>% 
    filter(!is.na(lib_drug)) %>% 
    select(-anchor, DRUG_ID_lib = DRUG_ID, CONC = CONC)
  
  anchors <- normalized_data %>% 
    filter(!is.na(anchor)) %>% 
    select(-lib_drug, -dose, DRUG_ID_anch = DRUG_ID, CONC_anch = CONC)
  if (nrow(anchors) > 0){
    suppressMessages(normalized_data <- full_join(libraries, anchors))
  }
  else {
    normalized_data <- libraries
  }

  normalized_data <- left_join(normalized_data, nc1, by=("SCAN_ID"))
  normalized_data <- left_join(normalized_data, pc1, by=("SCAN_ID"))
  normalized_data <- normalized_data %>%
    mutate(normalized_intensity = ((INTENSITY - PC) / (NC - PC)))

  if(trim){
    normalized_data <- normalized_data %>%
      mutate(normalized_intensity = 
                (ifelse(normalized_intensity > 1, 1, normalized_intensity))) %>%
      mutate(normalized_intensity =
                (ifelse(normalized_intensity < 0, 0, normalized_intensity)))
  }
  
  normalized_data <- normalized_data %>% 
    mutate(norm_neg_pos = paste(neg_control, pos_control, sep = "+"))

  normalized_data <- normalized_data %>% 
    mutate(time_stamp = Sys.time()) %>%
    mutate(sw_version = set_package_used())

  return(normalized_data)
}


calcTagMean <- function(myDat, tag_name, mean_col_name = "tag_mean") {
  suppressMessages(check_for_tag <- left_join(myDat %>% 
                               select(SCAN_ID) %>% 
                               distinct(),
                             myDat %>% 
                               group_by(SCAN_ID) %>% 
                               filter(TAG == tag_name) %>% 
                               count()
                             )
  )

  e1 <- simpleError(paste("calcTagMean:", 
                          tag_name, 
                          "is not present for some or all of the SCAN_IDs in your data.", 
                          sep = " "))
  if (any(is.na(check_for_tag$n))){
    stop(e1)
  }
  
  tag_means <- myDat %>% 
    group_by(SCAN_ID) %>% 
    filter(TAG == tag_name) %>% 
    summarise(tag_mean = mean(INTENSITY, na.rm = T))
  
  e2 <- simpleError(paste("calcTagMean:", 
                          tag_name, 
                          "has a mean of NaN for some or all of the SCAN_IDs in your data.", 
                          sep = " "))
  if (any(is.nan(tag_means$tag_mean))){
    stop(e2)
  }
  
  tag_means <- tag_means %>% rename(!!mean_col_name := tag_mean)
  
  return(tag_means)
}


#' normalizeComboData
#' 
#' Now deprecated - using \code{normalizeData} will return identical results.
#' 
#' \code{normalizeComboData} returns normalized intensities for the drug treated
#'   wells with respect to controls - column \code{normalized_intensity}.
#'   
#' Like \code{normalizeData}, the \code{TAG} column is replaced with 
#'   \code{lib_drug} and \code{dose} columns, but also \code{treatment} (S or C)
#'   , \code{DRUG_ID_anch}, \code{CONC_anch} and \code{anchor} columns. 
#'   Combination treatments (C) are treated as if the library alone is to be 
#'   fitted.
#'
#' @param myDat a GDSC raw data data frame.
#' @param trim logical indicating whether to trim normalized values to the range
#' 0 to 1. default \code{(trim = T)}
#' @param neg_control The tag used to recognise a negative control well - the 
#'   upper end of the dynamic range.
#' @param pos_control The tag used to recognise a positive control well - the 
#'   lower end of the dynamic range.
#' 
#' @seealso  \code{\link{removeFailedDrugs}},  \code{\link{removeMissingDrugs}},
#'   \code{\link{setConcsForNlme}},  \code{\link{prepNlmeData}}, 
#'   \code{\link{normalizeData}}
#' 
#' @examples
#' data("gdsc_example") # Need a combo example here
#' gdsc_example <- removeFailedDrugs(gdsc_example)
#' gdsc_example <- removeMissingDrugs(gdsc_example)
#' gdsc_example <- normalizeComboData(gdsc_example)
#' gdsc_example <- setConcsForNlme(gdsc_example)
#' nlme_data <- prepNlmeData(gdsc_example, "COSMIC_ID")
#'
#' @export
normalizeComboData <- function(myDat, trim = T, neg_control = 'NC-1',
                               pos_control = 'B'){
  normalized_data <- normalizeData(myDat, trim, neg_control, pos_control)
  return(normalized_data)
}

averageControlData <- function(screen_data, pos_control, neg_control){
  nc1 <- screen_data %>% group_by(SCAN_ID) %>% filter(TAG == neg_control) %>% 
    summarise(NC = mean(INTENSITY))
  pc1 <- screen_data %>% group_by(SCAN_ID) %>% filter(TAG == pos_control) %>%
    summarise(PC = mean(INTENSITY))
  average_controls <- suppressMessages(inner_join(nc1, pc1))
  return(average_controls)
}


#' normalizeMultiComboData
#' 
#' Normalize combination data and condense data for >2 treatments
#' 
#' Some combination treatments have more than 2 drugs. This function will 
#' condense the treatments into an anchor + a library, e.g., A10, A12, L3, L4
#' becomes A10_A12 with L3_L4.
#' 
#' The concentration for the analysis is at the moment chosen with the lowest 
#' anchor/library number taking precendence, so in the above example the
#' concentration of A10 is used as the concentration of the combined A10_A12 anchor.
#' 
#' This code will run slower than normalizeComboData if used for simple 2 drug 
#' combinations.
#' 
#' @param screen_data a GDSC raw data data frame.
#' @param trim logical indicating whether to trim normalized values to the range
#' 0 to 1. default \code{(trim = T)}
#' @param neg_control The tag used to recognise a negative control well - the 
#'   upper end of the dynamic range.
#' @param pos_control The tag used to recognise a positive control well - the 
#'   lower end of the dynamic range.
#' 
#' @seealso  \code{\link{normalizeComboData}}
#'
#' @export
normalizeMultiComboData <- function(screen_data, trim = T, neg_control = 'NC-1', pos_control = 'B'){
  condensed_screen_data <- condenseScreenData(screen_data = screen_data, 
                     neg_control = neg_control, 
                     pos_control = pos_control)
  
  normalized_data <- condensed_screen_data %>%
    mutate(normalized_intensity = ((INTENSITY - PC) / (NC - PC)))
  
  if(trim){
    normalized_data <- normalized_data %>%
      mutate(normalized_intensity = 
                (ifelse(normalized_intensity > 1, 1, normalized_intensity))) %>%
      mutate(normalized_intensity =
                (ifelse(normalized_intensity < 0, 0, normalized_intensity)))
  }

  normalized_data <- normalized_data %>% 
    mutate(norm_neg_pos = paste(neg_control, pos_control, sep = "+")) %>%
    mutate(time_stamp =  Sys.time()) %>%
    mutate(sw_version =  set_package_used())
  
  return(normalized_data)
}


condenseScreenData <- function(screen_data, neg_control, pos_control){
  
  average_controls <- averageControlData(screen_data, 
                                         neg_control = neg_control,
                                         pos_control = pos_control)
  
  drugset_layouts <- screen_data %>% 
    select(DRUGSET_ID, POSITION, TAG, DRUG_ID, CONC) %>% 
    distinct()
  
  # drugset_layouts <- drugset_layouts %>% 
  #   group_by_(~DRUGSET_ID) %>% 
  #   do(condensed_layouts = condenseDruggedLayout(.)) %>%
  #   tidyr::unnest(condensed_layouts)
  drugset_layouts <- drugset_layouts %>% 
    split(.$DRUGSET_ID) %>% 
    purrr::map_df(condenseDruggedLayout)
  
  screen_data <- screen_data %>% select(RESEARCH_PROJECT,
                         BARCODE, 
                         SCAN_ID,
                         DATE_CREATED,
                         SCAN_DATE,
                         CELL_ID,
                         COSMIC_ID,
                         MASTER_CELL_ID,
                         CELL_LINE_NAME,
                         SEEDING_DENSITY,
                         DRUGSET_ID,
                         ASSAY,
                         DURATION,
                         POSITION,
                         INTENSITY) %>%
    distinct()
  
  condensed_screen_data <- 
    inner_join(screen_data, drugset_layouts,
               by = c("DRUGSET_ID", "POSITION")) %>%
    inner_join(average_controls, by = c("SCAN_ID"))
  
  return(condensed_screen_data)
}

condenseDruggedLayout <- function(drugset_layout){
  drugset_id <- unique(drugset_layout$DRUGSET_ID)
  drugset_layout <- drugset_layout %>%
    filter(grepl("^(A|L)", TAG))
  if (nrow(drugset_layout) == 0){
    stop(paste("No drugged wells in drugset ", drugset_id, sep = ""))
  }
  
  condensed_drugged_layout <- drugset_layout %>%
    split(.$POSITION) %>%
    purrr::map_df(condenseWellPosition) 
  return(condensed_drugged_layout)
}

condenseWellPosition <- function(position_data){
  position_data <- position_data %>% filter(grepl("^(A|L)", TAG))
  if (nrow(position_data) == 0){
    stop("Attempting condenseWellPosition for non-drugged wells (A or L)")
  }
  drugset_id <- position_data %>% select(DRUGSET_ID) %>% distinct()
  well_pos <- position_data %>% select(POSITION) %>% distinct()
  
  if (nrow(drugset_id) != 1) {
    stop("Attempting to condense positions from different drugsets.")
  }
  if (nrow(well_pos) != 1) {
    stop("Attempting to condense different well positions.")
  }
  
  position_data <- position_data %>% mutate(lib_drug = sub("((L|R)\\d+)(-D\\d+)?-(S|C)", "\\1", TAG),
                                             lib_drug = ifelse(grepl("^A.+", lib_drug), yes = NA, no = lib_drug),
                                             anchor = sub("(A\\d+)-(S|C)", "\\1", TAG),
                                             anchor = ifelse(grepl("^(L|R).+", anchor), yes = NA, no = anchor),
                                             dose = sub("(A|L|R)\\d+-?(D\\d+)?-(S|C)", "\\2", TAG),
                                             dose = ifelse(dose == "", yes = NA, no = dose), 
                                             treatment = sub("((A|L|R)\\d+)(-D\\d+)?-(S|C)", "\\4", TAG), 
                                             treatment = ifelse(treatment == "", yes = NA, no = treatment))
  
  
  libraries <- position_data %>% 
      filter( !is.na(lib_drug)) %>%
      select(DRUG_ID, CONC, lib_drug, dose, treatment) %>%
      # This arrange is important for the rare case where the same DRUG_ID is 
      # used multiple times in as an S treatment (e.g. L21-S and L22-S combined
      # into a single well as 'S'-treatment, where L21 and L22 are the same DRUG_ID)
      arrange( lib_drug)  
  
  if (nrow(libraries) == 0){
    libraries <- data_frame(lib_drug = as.character(NA),
                            dose = as.character(NA),
                            treatment_lib = as.character(NA),
                            DRUG_ID_lib = as.character(NA),
                            CONC_lib = as.character(NA), 
                            CONC_lib_analysis = as.numeric(NA))
  } else {
    # # Check the dose level is the same for all libraries
    # if (length(unique(libraries$dose)) > 1){
    #   stop(paste("Non matching dose levels for library drugs: drugset ",
    #              libraries$DRUGSET_ID, ", position ", libraries$POSITION,
    #              sep = ""))
    # }
    # Check all treatments are the same
    if (length(unique(libraries$treatment)) > 1){
      stop(paste("Non matching treatments for library drugs: drugset ",
                 libraries$DRUGSET_ID, ", position ", libraries$POSITION, 
                 sep = ""))
    }
    
    # Select the arbitrary analysis concentration by the library number
    # L1 before L2 before L3 etc.
    lib_conc_analysis <- libraries %>% 
      mutate(lib_number = as.numeric(sub("L(\\d+)", "\\1", lib_drug))) %>% 
        select(lib_number, CONC) %>%
        arrange(lib_number) %>% 
        slice(1) %>% 
        select(CONC)
    
    libraries <- libraries %>% 
      mutate(lib_drug = paste(.$lib_drug, collapse = "|"), 
             dose = paste(.$dose, collapse = "|"),
             treatment_lib = treatment,
             DRUG_ID_lib = paste(.$DRUG_ID, collapse = "|"), 
             CONC_lib = paste(.$CONC, collapse = "|"),
             CONC_lib_analysis = as.numeric(lib_conc_analysis$CONC)
      ) %>% 
      select(-CONC, -DRUG_ID, -treatment) %>%
      distinct()
    
    if(nrow(libraries) > 1){
      stop(paste("Failed to condense library drugs into one row:",
                 libraries$DRUGSET_ID, ", position ", libraries$POSITION, 
                 sep = ""))
    }
  }
  
  anchors <- position_data %>% 
      filter( !is.na(anchor)) %>%
      select(DRUG_ID, CONC, anchor, treatment) %>%
      
      # This arrange is important for the rare case where the same DRUG_ID is 
      # used multiple times in as an S treatment (e.g. A21-S and A22-S combined
      # into a single well as 'S'-treatment, where A21 and A22 are the same DRUG_ID)
      arrange(anchor) 
  
  if (nrow(anchors) == 0){
    anchors <- data_frame(anchor = as.character(NA),
                          treatment_anch = as.character(NA),
                          DRUG_ID_anch = as.character(NA),
                          CONC_anch = as.character(NA), 
                          CONC_anch_analysis = as.numeric(NA))
  } else{
    # Check all treatments are the same
    if (length(unique(anchors$treatment)) > 1){
      stop(paste("Non matching treatments for library drugs: drugset ",
                 anchors$DRUGSET_ID, ", position ", anchors$POSITION, 
                 sep = ""))
    }
    
    anch_conc_analysis <- anchors %>% 
      mutate(anch_number = as.numeric(sub("A(\\d+)", "\\1", anchor))) %>% 
      select(anch_number, CONC) %>%
      arrange(anch_number) %>% 
      slice(1) %>% 
      select(CONC)

    anchors <- anchors %>% 
      mutate(anchor = paste(.$anchor, collapse = "|"), 
             treatment_anch = treatment,
             DRUG_ID_anch = paste(.$DRUG_ID, collapse = "|"),
             CONC_anch = paste(.$CONC, collapse = "|"),
             CONC_anch_analysis = as.numeric(anch_conc_analysis$CONC)
      ) %>% 
      select(-CONC, -DRUG_ID, -treatment) %>%
      distinct()
    
    if(nrow(anchors) > 1){
      stop(paste("Failed to condense anchor drugs into one row:",
                 anchors$DRUGSET_ID, ", position ", anchors$POSITION, 
                 sep = ""))
    }
  }
  
  condensed_position <- cbind(drugset_id, well_pos, libraries, anchors)
  condensed_position <- condensed_position %>% 
    mutate(treatment_lib =  ifelse(is.na(treatment_lib) && !is.na(treatment_anch), 
                                  yes = treatment_anch, 
                                  no = treatment_lib)) %>%
    mutate(treatment_anch =  ifelse(is.na(treatment_anch) && !is.na(treatment_lib), 
                                   yes = treatment_lib, 
                                   no = treatment_anch)) %>%
    mutate(treatment =  ifelse(treatment_lib == treatment_anch,
                              yes = treatment_lib, 
                              no = 'mismatch')
    ) %>%
    select(-treatment_lib, -treatment_anch) 
    if(nrow(condensed_position) != 1){
    stop(paste("Failed to condense library and anchor drugs into one row: drugset",
               condensed_position$DRUGSET_ID,
               ", position ", condensed_position$POSITION, 
               sep = ""))
  }
  
  if (condensed_position$treatment == 'mismatch'){
    stop(paste("Treatments (-C or -S) mismatch: drugset ", 
               condensed_position$DRUGSET_ID, ", position ", 
               condensed_position$POSITION, sep = "")
    )
  }
  return(condensed_position)
}


# condenseNormalizedPosition <- function(position_data){
#   
#   position_annotation <- position_data %>%
#     select(SCAN_ID,
#            BARCODE,
#            DATE_CREATED,
#            DRUGSET_ID,
#            CELL_LINE_NAME,
#            CELL_ID,
#            COSMIC_ID,
#            POSITION,
#            INTENSITY) %>%
#     distinct()
#   
#   if (nrow(position_annotation) != 1){
#     stop(paste("Cannot condense position data where annotation is not unique: scan_id ",
#                paste(foo$SCAN_ID, collapse = " "),
#                "; position ",
#                paste(foo$POSITION, collapse = " "),
#                sep = ""
#     )
#     )
#   }
#   
#   libraries <- position_data %>% filter(!is.na(lib_drug)) %>%
#     select(DRUG_ID, CONC, lib_drug, dose, treatment)
#   
#   if (nrow(libraries) == 0){
#     libraries <- position_annotation %>% 
#       mutate(lib_drug = NA,
#              dose = NA,
#              treatment_lib = NA,
#              DRUG_ID_lib = NA,
#              CONC_lib = NA, 
#              CONC_lib_analysis = NA
#       )
#   } else {
#     
#     lib_conc_analysis <- libraries %>% select(lib_drug, CONC) %>%
#       arrange(lib_drug) %>% slice(1) %>% select(CONC)
#     
#     libraries <- libraries %>% 
#       mutate(lib_drug = paste(.$lib_drug, collapse = "|"), 
#              dose = ifelse(length(unique(dose)) == 1,
#                            yes = unique(dose),
#                            no ='mismatch'),
#              treatment_lib = ifelse(length(unique(treatment)) == 1,
#                                     yes = unique(treatment),
#                                     no = 'mismatch'),
#              DRUG_ID_lib = paste(.$DRUG_ID, collapse = "|"), 
#              CONC_lib = paste(.$CONC, collapse = "|"),
#              CONC_lib_analysis = as.numeric(lib_conc_analysis$CONC)
#       ) %>% 
#       select(-CONC, -DRUG_ID, -treatment) %>%
#       distinct()
#     
#     stopifnot(nrow(libraries) < 2)
#     
#     if (nrow(libraries) == 1 && libraries$dose == 'mismatch'){
#       stop(paste("Library doses mismatch: scan id ", libraries$SCAN_ID, 
#                  ", drugset ", libraries$DRUGSET_ID,
#                  ", position ", libraries$POSITION, sep = "")
#       )
#     }
#     if (nrow(libraries) == 1 && libraries$treatment_lib == 'mismatch'){
#       stop(paste("Library treatments (-C or -S) mismatch: scan id ", libraries$SCAN_ID, 
#                  ", drugset ", libraries$DRUGSET_ID,
#                  ", position ", libraries$POSITION, sep = "")
#       )
#     }
#     libraries <- cbind(position_annotation, libraries)
#   }
#   
#   anchors <- position_data %>% filter(!is.na(anchor)) %>%
#     select(DRUG_ID, CONC, anchor, treatment)
#   
#   if (nrow(anchors) == 0){
#     anchors <- position_annotation %>% 
#       mutate(anchor = NA,
#              treatment_anch = NA,
#              DRUG_ID_anch = NA,
#              CONC_anch = NA, 
#              CONC_anch_analysis = NA
#       )
#   }
#   else{
#     anch_conc_analysis <- anchors %>% select(anchor, CONC) %>%
#       arrange(anchor) %>% slice(1) %>% select(CONC)
#     
#     anchors <- anchors %>% 
#       mutate(anchor = paste(.$anchor, collapse = "|"), 
#              treatment_anch = ifelse(length(unique(treatment)) == 1,
#                                      yes = unique(treatment),
#                                      no = 'mismatch'),
#              DRUG_ID_anch = paste(.$DRUG_ID, collapse = "|"),
#              CONC_anch = paste(.$CONC, collapse = "|"),
#              CONC_anch_analysis = as.numeric(anch_conc_analysis$CONC)
#       ) %>% 
#       select(-CONC, -DRUG_ID, -treatment) %>%
#       distinct()
#     
#     stopifnot(nrow(anchors) < 2)
#     if (nrow(anchors) == 1 && anchors$treatment_anch == 'mismatch'){
#       stop(paste("Anchor treatments (-C or -S) mismatch: scan id ", anchors$SCAN_ID, 
#                  ", drugset ", anchors$DRUGSET_ID,
#                  ", position ", anchors$POSITION, sep = "")
#       )
#     }
#     
#     anchors <- cbind(position_annotation, anchors)
#   }
#   
#   normalized_position_data <- suppressMessages(inner_join(
#     libraries, 
#     anchors)) %>%
#     mutate(treatment_lib = ifelse(is.na(treatment_lib) && !is.na(treatment_anch), 
#                                   yes = treatment_anch, 
#                                   no = treatment_lib)) %>%
#     mutate(treatment_anch = ifelse(is.na(treatment_anch) && !is.na(treatment_lib), 
#                                    yes = treatment_lib, 
#                                    no = treatment_anch)) %>%
#     mutate(treatment = ifelse(treatment_lib == treatment_anch,
#                               yes = treatment_lib, 
#                               no = 'mismatch')
#     ) %>%
#     select(-treatment_lib, -treatment_anch) 
#   
#   stopifnot(nrow(normalized_position_data) == 1)
#   
#   if (normalized_position_data$treatment == 'mismatch'){
#     stop(paste("Treatments (-C or -S) mismatch: scan id ", normalized_position_data$SCAN_ID, 
#                ", drugset ", normalized_position_data$DRUGSET_ID,
#                ", position ", normalized_position_data$POSITION, sep = "")
#     )
#   }
#   return(normalized_position_data)
# }

setMaxConc <- function(normalized_data, drug_specifiers){
  normalized_data %>% 
    group_by_at(drug_specifiers) %>%
    mutate(maxc =  max(CONC)) %>%
    ungroup()
}

setXFromConc <- function(normalized_data){
  # Instead of using djvMixedIC50::getXfromConcSeries
  normalized_data <- normalized_data %>%
    mutate(x =  (log(CONC / maxc) / log(2)) + 9)
  return(normalized_data)
}

#' Adds columns \code{maxc} and \code{x} to GDSC normalized-data data frame.
#' 
#' This is run prior to running prepNlmeData. The functionality is separate 
#' because the resulting `maxc` column might be used as a `drug_specifier`
#'  in prepNlmeData.
#' 
#' \code{maxc} is the maximum micromolar concentration of the treatment drug.
#' \code{x} is the conversion of the micromolar screening concentration to 
#' a range where the maximum dose = 9
#' 
#' @param normalized_data a GDSC normalized-data data frame.
#' @param group_conc_ranges logical. If TRUE then instances where the same drug
#'  has been used in different concentration ranges will be grouped together
#'   with a single maxc. Default is FALSE.
#' @param conc_col a string to identify the column used for the concentration 
#'   values - useful for combination drug treatments.
#' 
#' @seealso  \code{\link{removeFailedDrugs}},  \code{\link{removeMissingDrugs}},
#'   \code{\link{normalizeData}},  \code{\link{prepNlmeData}}
#' 
#' @examples
#' data("gdsc_example")
#' gdsc_example <- removeFailedDrugs(gdsc_example)
#' gdsc_example <- removeMissingDrugs(gdsc_example)
#' gdsc_example <- normalizeData(gdsc_example)
#' gdsc_example <- setConcsForNlme(gdsc_example)
#' nlme_data <- prepNlmeData(gdsc_example, "COSMIC_ID")
#' 
#' @export
setConcsForNlme <- function(normalized_data, 
                            group_conc_ranges = F,
                            conc_col = "CONC"
                           ) {
  if (group_conc_ranges){
    drug_specifiers <- "DRUG_ID_lib"
    message("Grouping all dilution series per DRUG_ID_lib to get maximum
            concentrations and to set x values.")
  }
  else {
    drug_specifiers = c("DRUGSET_ID", "lib_drug")
    message("Different dilution series per DRUG_ID_lib are being treated
            separately to get maximum concentration and to set x values.")
  }
  
  if(conc_col != "CONC" && !("CONC" %in% names(normalized_data)) ){
      normalized_data["CONC"] <- normalized_data[conc_col]
  }
    
  normalized_data <- normalized_data %>% 
    setMaxConc(drug_specifiers = drug_specifiers) %>% 
    setXFromConc()
  
  return(normalized_data)
}

#' Converts GDSC normalized-data data frame to input format for nlme curve fit.
#' 
#' @param normalized_data a GDSC normalized-data data frame.
#' @param cl_id spcifies which cell line identifier to be used for nlme input - 
#'   currently "COSMIC_ID" or "CELL_ID"
#' @param drug_specifiers a character vector containing the names of the columns 
#'   to be taken from the \code{normalized_data} and combined to make the new 
#'   \code{drug} column. The drug column will be used in the nlme model.
#' @param include_combos logical indicating whether or not to include -C tags in
#'   the output - you will need to be careful with the drug_specifiers if T, e.g., 
#'   \code{drug_specifiers = c("SCAN_ID", "lib_drug", "anchor")}
#' 
#' @return data frame with columns \code{DRUG_ID, CELL_LINE_NAME, CL, maxc, x,
#'  y, drug, SCAN_ID, norm_neg_pos, drug_spec}
#' \code{CL} is the cell line identifier.
#' \code{x} is the conversion of the micromolar screening concentration to 
#' a range where the maximum dose = 9
#' \code{y} is the normalized intensity from the experiment (sometimes called
#'  viability)
#' \code{drug} a concatenation of the drug_specifier columns 
#'   (default DRUG_ID and maxc).
#' \code{SCAN_ID} The id of the plate scan used to provide the raw data.
#' \code{norm_neg_pos} The identifiers of the control tags used to normalise the
#'  data separated by '+'.
#' \code{drug_spec} the column names used to make the drug column separated by
#'  '+'. The original columns will be also included as separate columns.
#' 
#' @examples
#' data("gdsc_example")
#' gdsc_example <- removeFailedDrugs(gdsc_example)
#' gdsc_example <- removeMissingDrugs(gdsc_example)
#' gdsc_example <- normalizeData(gdsc_example)
#' gdsc_example <- setConcsForNlme(gdsc_example)
#' nlme_data <- prepNlmeData(gdsc_example, 
#'                           cl_id = "COSMIC_ID",
#'                            drug_specifiers = c("DRUGSET_ID", "lib_drug"))
#' 
#' 
#' @seealso  \code{\link{removeFailedDrugs}},  \code{\link{removeMissingDrugs}},
#'   \code{\link{normalizeData}},  \code{\link{setConcsForNlme}}, 
#'   \code{\link{fitModelNlmeData}}
#'  
#' @export
prepNlmeData <- function(normalized_data, cl_id = "",
                         drug_specifiers = c("DRUG_ID_lib", "maxc"),
                         include_combos = F){
  
  # Check that setConcsForNlme has been run first. 
  e1 <- simpleError("No maxc or x columns in the normalized_data. Run setConcsForNlme() before prepNlmeData")
  if (is.null(normalized_data$maxc) || is.null(normalized_data$x)){
    stop(e1)
  }
  
  # Check normalized_data has the required columns
  e2 <- simpleError("Your normalized data does not contain the columns specified to make the drug column.")
  if (!all({{drug_specifiers}} %in% names(normalized_data))){
    stop(e2)
  }
  
  if(include_combos){
    normalized_data <- normalized_data %>% 
      filter(!is.na(lib_drug))
  }
  else { # Don't fit single Anchors, e.g., A1-S
    normalized_data <- normalized_data %>% 
      filter(treatment == 'S', !is.na(lib_drug))
  }
  
  if(!(cl_id) %in% c("COSMIC_ID", "CELL_ID", "MASTER_CELL_ID")){
      stop('choose a suitable cl_id: COSMIC_ID, MASTER_CELL_ID or CELL_ID')

  }
  nlme_data <- normalized_data %>% 
      select(CELL_LINE_NAME, CL = one_of({{cl_id}}), maxc, x, y = normalized_intensity, 
          one_of({{drug_specifiers}}), BARCODE, SCAN_ID, POSITION, DRUGSET_ID, norm_neg_pos) %>% 
    mutate(CL_SPEC = {{cl_id}}) %>% 
    tidyr::unite(col = drug, {{drug_specifiers}}, sep = "_", remove = F) %>% 
    mutate(drug_spec = paste({{drug_specifiers}}, collapse = "+"),
           y =  1 - y,
           time_stamp = normalized_data$time_stamp[1],
           sw_version = normalized_data$sw_version[1]) %>% 
    arrange(CL, drug, x)
  
  # Check that there is a 1-to-1 correspondence between the drug and maxc for that drug
  maxc_check <- nlme_data %>% distinct(drug, maxc) %>% count(drug) %>% filter(n > 1)
  
  if (nrow(maxc_check) > 0){
      warning(paste("There is more than one maximum concentration for drug ",
                    maxc_check$drug,
                    "",
                    sep = "")
              )
    }
  
  
  # Add extra annotation
  if (!is.null(normalized_data$RESEARCH_PROJECT)){
      nlme_data <- nlme_data %>%
          left_join(normalized_data %>% distinct(BARCODE, RESEARCH_PROJECT),
                    by = "BARCODE")
  }
  
  return(nlme_data)
}

#' Utility function to get a time and data stamp for use in file names.
#' 
#' @return character string in date format "%d%b%y_%H%M"
#' 
#' @examples 
#' time_stamp = getTimeStamp()
#' 
#' @export
getTimeStamp <- function(){
  time_stamp <- format(Sys.time(), "%d%b%y_%H%M")
  return(time_stamp)
}

#' Get the package name and version information
#'
#' @return package name and version used to prepare the data as a string.
#' 
set_package_used <- function(){
  package_used <- paste(
    utils::packageDescription("gdscIC50")$Package,
    utils::packageDescription("gdscIC50")$Version,
    sep = "_")
  return(package_used)
}


