context("Test setConcsForNlme")

test_data<- readRDS("norm_single_agents_test.rds")

# Add extra data for testing group conc ranges changes
barcode1235 <- test_data %>% 
  filter(BARCODE == 1234) %>% 
  mutate(BARCODE = 1235, SCAN_ID = 1200, DRUGSET_ID = 102, CELL_ID = 211) %>% 
  mutate(CONC = ifelse(DRUG_ID == 9002, CONC / 10, CONC))

barcode1235$INTENSITY <- sample(barcode1235$INTENSITY, length(barcode1235$INTENSITY), replace = FALSE)

test_data <- test_data %>% rbind(barcode1235)

norm_data <- normalizeData(test_data)
# norm_data <- setConcsForNlme(norm_data)

# DRUG_ID 9002 with maxc of:  20 (L50, ds101); 2 (L50, ds102); 10 (L48, ds202).

test_that("group_conc_ranges=F returns separate maxc for drug id 9002",{
  set_concs_data <- setConcsForNlme(norm_data, group_conc_ranges = F)
  number_of_maxc <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, "maxc"]))
  number_of_cell_lines <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, c("MASTER_CELL_ID")]))
  
  maxc_L48 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 202 &
                                    set_concs_data$lib_drug == "L48" & 
                                    set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_L50 <- unique(set_concs_data[
    # set_concs_data$DRUGSET_ID == 101 &
                                      set_concs_data$lib_drug == "L50" & 
                                      set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  testthat::expect_true(number_of_maxc == 3)
  expect_true(number_of_maxc > number_of_cell_lines)
  testthat::expect_true(nrow(maxc_L50) == 2)
  testthat::expect_true(nrow(maxc_L48) == 1)
  # testthat::expect_true(maxc_L48$maxc != maxc_L50$maxc)
})

test_that("group_conc_ranges=T returns a single maxc for drug id 9002", {
  set_concs_data <- setConcsForNlme(norm_data, group_conc_ranges = T)
  number_of_maxc <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, "maxc"]))
  maxc_L48 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 202 &
                                      set_concs_data$lib_drug == "L48" & 
                                      set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                      set_concs_data$lib_drug == "L50" & 
                                      set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  number_of_cell_lines <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, "MASTER_CELL_ID"]))
  # Expect there to be different maxc per cell line
  expect_true(number_of_maxc == 2)
  expect_true(number_of_maxc == number_of_cell_lines)
  testthat::expect_true(nrow(maxc_L50) == 1)
  testthat::expect_true(nrow(maxc_L48) == 1)
  expect_true(maxc_L48$maxc != maxc_L50$maxc)
  
})

test_that("group_conc_ranges=F returns the same set of x values for L50 ds101 and L50 ds102",{
  set_concs_data <- setConcsForNlme(norm_data, group_conc_ranges = F)
  x_L50_101 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                             set_concs_data$lib_drug == "L50" & 
                                             set_concs_data$DRUG_ID_lib == 9002, "x"])
  x_L50_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                             set_concs_data$lib_drug == "L50" & 
                                             set_concs_data$DRUG_ID_lib == 9002, "x"])
  testthat::expect_true(isTRUE(all.equal(sort(x_L50_101$x), sort(x_L50_102$x))))
})

test_that("group_conc_ranges=T returns a different set of x values for L50 ds101 and L50 ds102",{
  set_concs_data <- setConcsForNlme(norm_data, group_conc_ranges = T)
  x_L50_101 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  x_L50_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  testthat::expect_false(isTRUE(all.equal(sort(x_L50_101$x), sort(x_L50_102$x))))
})

test_that("group_conc_ranges=T returns different maxc if dividing cell lines by expansion_id",{
  set_concs_data <- setConcsForNlme(norm_data, group_conc_ranges = T, cell_line_spec = "MASTER_CELL_ID")
  number_of_maxc_master_id <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002 &
                                                          set_concs_data$MASTER_CELL_ID == 111, "maxc"]))
  testthat::expect_true(number_of_maxc_master_id == 1)
  
  set_concs_data <- setConcsForNlme(norm_data, group_conc_ranges = T, cell_line_spec = "CELL_ID")
  number_of_maxc_cell_id<- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002 &
                                                        set_concs_data$MASTER_CELL_ID == 111, "maxc"]))
  testthat::expect_true(number_of_maxc_cell_id == 2)
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



