context("Test setConcsForNlme")

test_data<- readRDS("norm_single_agents_test.rds")
norm_data <- normalizeData(test_data)
norm_data <- setConcsForNlme(norm_data)

test_that("group_conc_ranges=F returns separate maxc for drug id 9002",{
  set_concs_data <- setConcsForNlme(norm_data, group_conc_ranges = F)
  number_of_maxc <- unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_L48 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 202 &
                                    set_concs_data$lib_drug == "L48" & 
                                    set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                      set_concs_data$lib_drug == "L50" & 
                                      set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  testthat::expect_true(nrow(number_of_maxc) == 2)
  testthat::expect_false(maxc_L48$maxc == maxc_L50$maxc)
})

test_that("group_conc_ranges=T returns a single maxc for drug id 9002", {
  set_concs_data <- setConcsForNlme(norm_data, group_conc_ranges = T)
  number_of_maxc <- unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_L48 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 202 &
                                      set_concs_data$lib_drug == "L48" & 
                                      set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                      set_concs_data$lib_drug == "L50" & 
                                      set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  expect_true(nrow(number_of_maxc) == 1)
  expect_true(maxc_L48$maxc == maxc_L50$maxc)
  
})

test_that("group_conc_ranges=F returns the same set of x values for L48 and L50",{
  set_concs_data <- setConcsForNlme(norm_data, group_conc_ranges = F)
  x_L48 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 202 &
                                             set_concs_data$lib_drug == "L48" & 
                                             set_concs_data$DRUG_ID_lib == 9002, "x"])
  x_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                             set_concs_data$lib_drug == "L50" & 
                                             set_concs_data$DRUG_ID_lib == 9002, "x"])
  testthat::expect_true(isTRUE(all.equal(sort(x_L48$x), sort(x_L50$x))))
})

test_that("group_conc_ranges=T returns a different set of x values for L48 and L50",{
  set_concs_data <- setConcsForNlme(norm_data, group_conc_ranges = T)
  x_L48 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 202 &
                                   set_concs_data$lib_drug == "L48" & 
                                   set_concs_data$DRUG_ID_lib == 9002, "x"])
  x_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                   set_concs_data$lib_drug == "L50" & 
                                   set_concs_data$DRUG_ID_lib == 9002, "x"])
  testthat::expect_false(isTRUE(all.equal(sort(x_L48$x), sort(x_L50$x))))
})

test_that("group_conc_ranges=T returns an informative message",{
  testthat::expect_message(
    setConcsForNlme(norm_data, group_conc_ranges = T),
    "^Grouping all dilution series per DRUG_ID_lib to get maximum\\s+concentrations and to set x values", perl = T)
})

test_that("group_conc_ranges=F returns an informative message",{
  testthat::expect_message(
    setConcsForNlme(norm_data, group_conc_ranges = F),
    "^Different dilution series per DRUG_ID_lib are being treated\\s+separately to get maximum concentration and to set x values", perl = T
  )
})



