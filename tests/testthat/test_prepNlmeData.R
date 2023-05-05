context("Test prepNlmeData")
# library(dplyr)

test_data<- readRDS("norm_single_agents_test.rds")
norm_data <- normalizeData(test_data)

context("Test prepNlmeData")

test_that("prepNlmeData gives an error if setGroupsForNlme has not been run first", {
  testthat::expect_error(
    prepNlmeData(norm_data))
})

test_that("y is equal to 1 - viabilty", {
  norm_data <- setGroupsForNlme(norm_data)
  nlme_data <- prepNlmeData(norm_data)
  test_data <- merge(norm_data[c("BARCODE", "POSITION", "normalized_intensity")],
                          nlme_data[c("BARCODE", "POSITION", "y")])
  testthat::expect_identical(1 - test_data$normalized_intensity, test_data$y)
})

# test_that("prepNlmeData gives an error if no cl_id is chosen", {
#   norm_data <- setDrugsForNlme(norm_data)
#   testthat::expect_error(
#     prepNlmeData(norm_data),
#     "choose a suitable cl_id: COSMIC_ID, MASTER_CELL_ID or CELL_ID"
#     )
# })

# test_that("prepNlmeData gives an error if a non-column drug_specifier is chosen", {
#   norm_data <- setConcsForNlme(norm_data)
#   testthat::expect_error(
#     prepNlmeData(norm_data,
#                  cl_id = "COSMIC_ID",
#                  drug_specifiers = ("drug_identifier")),
#     "Your normalized data does not contain the columns specified to make the drug column."
#     )
# })
# 
# 
# test_that("prepNlmeData gives a warning if setConcsForNlme has not grouped 
#           concentrations and drug specifiers doesn't use maxc when there is a
#           drug at two ranges (9002).", {
#   norm_data <- setConcsForNlme(norm_data, group_conc_ranges = F)
#   testthat::expect_warning(
#     prepNlmeData(norm_data, cl_id = "COSMIC_ID", drug_specifiers = c("DRUG_ID_lib")),
#     "There is more than one maximum concentration for drug "
#     )
# })
