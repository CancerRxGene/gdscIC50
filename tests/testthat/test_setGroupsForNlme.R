context("Test setGroupsForNlme")
# --- Function prototype: ----
# gdscIC50::setGroupsForNlme(norm_data, 
#                           drug_spec = "DRUG_ID_lib", 
#                           cell_line_spec = "MASTER_CELL_ID", 
#                           conc_col = "CONC")

test_data<- readRDS("norm_single_agents_test.rds")

# Add extra data for testing group conc ranges changes
barcode1235 <- test_data %>% 
  filter(BARCODE == 1234) %>% 
  mutate(BARCODE = 1235, SCAN_ID = 1200, DRUGSET_ID = 102, CELL_ID = 211) %>% 
  mutate(CONC = ifelse(DRUG_ID == 9002, CONC / 10, CONC)) %>% 
  mutate(DRUG_ID = ifelse(grepl("^L52-D\\d+-S", TAG), 9002, DRUG_ID))


barcode1235$INTENSITY <- sample(barcode1235$INTENSITY, length(barcode1235$INTENSITY), replace = FALSE)

test_data <- test_data %>% rbind(barcode1235)

norm_data <- normalizeData(test_data)
# norm_data <- setConcsForNlme(norm_data)

# DRUG_ID 9002 with maxc of:  20 (L50, ds101); 2 & 5 (L50 & L52, ds102); 10 (L48, ds202).
# Bill sees ds101 and 102; Ben sees ds202

test_that('drug_spec = "DRUG_ID_lib" returns 2 maxc for drug id 9002',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = "DRUG_ID_lib", 
                                    cell_line_spec = "MASTER_CELL_ID",
                                    conc_col = "CONC")
  number_of_maxc <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, "maxc"]))
  maxc_ds202_L48 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 202 &
                                      set_concs_data$lib_drug == "L48" & 
                                      set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_ds101_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                      set_concs_data$lib_drug == "L50" & 
                                      set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_ds102_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                            set_concs_data$lib_drug == "L50" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_ds102_L52 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                            set_concs_data$lib_drug == "L52" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_Bill_9002 <- unique(set_concs_data[set_concs_data$CELL_LINE_NAME == "Bill" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_Ben_9002 <- unique(set_concs_data[set_concs_data$CELL_LINE_NAME == "Ben" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  number_of_cell_lines <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, "MASTER_CELL_ID"]))
  # Expect there to be different maxc per cell line
  testthat::expect_true(number_of_maxc == 2)
  testthat::expect_true(number_of_maxc == number_of_cell_lines)
  
  testthat::expect_true(maxc_ds101_L50$maxc == 20)
  testthat::expect_true(maxc_ds102_L50$maxc == 20)
  testthat::expect_true(maxc_ds102_L52$maxc == 20)
  testthat::expect_true(maxc_ds202_L48$maxc == 10)
  
  testthat::expect_true(maxc_Bill_9002 == 20)
  testthat::expect_true(maxc_Ben_9002 == 10)
  
  testthat::expect_identical(maxc_Bill_9002$maxc, max(maxc_ds101_L50$maxc, maxc_ds102_L50$maxc, maxc_ds102_L52$maxc))
  testthat::expect_identical(maxc_Ben_9002$maxc, maxc_ds202_L48$maxc)
})


test_that('drug_spec = "DRUG_ID_lib" returns a different set of x values for L50 ds101 and L50 ds102',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = "DRUG_ID_lib", 
                                    cell_line_spec = "MASTER_CELL_ID",
                                    conc_col = "CONC") 
  x_L50_101 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  x_L50_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  testthat::expect_false(isTRUE(all.equal(sort(x_L50_101$x), sort(x_L50_102$x))))
})


test_that('drug_spec = "DRUG_ID_lib" returns a different set of x values for L50 ds101 and L52 ds102',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = "DRUG_ID_lib", 
                                    cell_line_spec = "MASTER_CELL_ID",
                                    conc_col = "CONC") 
  x_L50_101 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  x_L52_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L52" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  testthat::expect_false(isTRUE(all.equal(sort(x_L50_101$x), sort(x_L52_102$x))))
})


test_that('drug_spec = "DRUG_ID_lib" returns a different set of x values for L50 ds102 and L52 ds102',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = "DRUG_ID_lib", 
                                    cell_line_spec = "MASTER_CELL_ID",
                                    conc_col = "CONC") 
  x_L50_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  x_L52_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L52" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  testthat::expect_false(isTRUE(all.equal(sort(x_L50_102$x), sort(x_L52_102$x))))
})



test_that('drug_spec = c("DRUGSET_ID", "lib_drug") returns 4 separate maxc for drug id 9002',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = c("DRUGSET_ID", "lib_drug"), 
                                    cell_line_spec = "MASTER_CELL_ID",
                                    conc_col = "CONC"
  )
  number_of_maxc <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, "maxc"]))
  maxc_ds202_L48 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 202 &
                                            set_concs_data$lib_drug == "L48" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_ds101_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                            set_concs_data$lib_drug == "L50" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_ds102_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                            set_concs_data$lib_drug == "L50" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_ds102_L52 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                            set_concs_data$lib_drug == "L52" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_Bill_9002 <- unique(set_concs_data[set_concs_data$CELL_LINE_NAME == "Bill" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_Ben_9002 <- unique(set_concs_data[set_concs_data$CELL_LINE_NAME == "Ben" & 
                                           set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  
  unique_specifiers <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, 
                                                     c("MASTER_CELL_ID", "DRUGSET_ID", "lib_drug")]))
  
  # Expect there to be different maxc per cell line
  expect_true(number_of_maxc == 4)
  expect_true(number_of_maxc == unique_specifiers)
  
  testthat::expect_true(nrow(maxc_Bill_9002) == 3)
  testthat::expect_true(nrow(maxc_Ben_9002) == 1)
  
  testthat::expect_true(maxc_ds101_L50$maxc == 20)
  testthat::expect_true(maxc_ds102_L50$maxc == 2)
  testthat::expect_true(maxc_ds102_L52$maxc == 5)
  testthat::expect_true(maxc_ds202_L48$maxc == 10)
  
  testthat::expect_equal(length(setdiff(maxc_Bill_9002$maxc, c(maxc_ds101_L50$maxc,
                                                       maxc_ds102_L50$maxc,
                                                       maxc_ds102_L52$maxc)
                                )),
                         0)
  testthat::expect_identical(maxc_Ben_9002$maxc, maxc_ds202_L48$maxc)
})


test_that('drug_spec = c("DRUGSET_ID", "lib_drug") returns the same set of x values for L50 ds101 and L50 ds102',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = c("DRUGSET_ID", "lib_drug"), 
                                    cell_line_spec = "MASTER_CELL_ID",
                                    conc_col = "CONC") 
  x_L50_101 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  x_L50_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  # testthat::expect_true(isTRUE(all.equal(sort(x_L50_101$x), sort(x_L50_102$x))))
  expect_equal(sort(x_L50_101$x), sort(x_L50_102$x))
  
})


test_that('drug_spec = c("DRUGSET_ID", "lib_drug") returns the same set of x values for L50 ds101 and L52 ds102',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = c("DRUGSET_ID", "lib_drug"), 
                                    cell_line_spec = "MASTER_CELL_ID",
                                    conc_col = "CONC") 
  x_L50_101 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  x_L52_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L52" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  # testthat::expect_true(isTRUE(all.equal(sort(x_L50_101$x), sort(x_L52_102$x))))
  expect_equal(sort(x_L50_101$x), sort(x_L52_102$x))
})


test_that('drug_spec = c("DRUGSET_ID", "lib_drug") returns the same set of x values for L50 ds102 and L52 ds102',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = c("DRUGSET_ID", "lib_drug"), 
                                    cell_line_spec = "MASTER_CELL_ID",
                                    conc_col = "CONC") 
  x_L50_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  x_L52_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L52" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  # testthat::expect_true(isTRUE(all.equal(sort(x_L50_102$x), sort(x_L52_102$x))))
  expect_equal(sort(x_L50_102$x), sort(x_L52_102$x))
})

# DRUG_ID 9002 with maxc of:  20 (L50, ds101); 2 & 5 (L50 & L52, ds102); 10 (L48, ds202).
# Bill sees ds101 and 102; Ben sees ds202

test_that('drug_spec = c("BARCODE", "DRUG_ID_lib") returns 3 separate maxc for drug id 9002',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = c("BARCODE", "DRUG_ID_lib"), 
                                    cell_line_spec = "MASTER_CELL_ID",
                                    conc_col = "CONC"
  )
  number_of_maxc <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, "maxc"]))
  maxc_ds202_L48 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 202 &
                                            set_concs_data$lib_drug == "L48" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_ds101_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                            set_concs_data$lib_drug == "L50" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_ds102_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                            set_concs_data$lib_drug == "L50" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_ds102_L52 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                            set_concs_data$lib_drug == "L52" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_Bill_9002 <- unique(set_concs_data[set_concs_data$CELL_LINE_NAME == "Bill" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_Ben_9002 <- unique(set_concs_data[set_concs_data$CELL_LINE_NAME == "Ben" & 
                                           set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  
  unique_specifiers <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, 
                                                  c("MASTER_CELL_ID", "BARCODE", "DRUG_ID_lib")]))
  
  # Expect there to be different maxc per cell line
  expect_true(number_of_maxc == 3)
  expect_true(number_of_maxc == unique_specifiers)
  
  testthat::expect_true(nrow(maxc_Bill_9002) == 2)
  testthat::expect_true(nrow(maxc_Ben_9002) == 1)
  
  testthat::expect_true(maxc_ds101_L50$maxc == 20)
  testthat::expect_true(maxc_ds102_L50$maxc == 5)
  testthat::expect_true(maxc_ds102_L52$maxc == 5)
  testthat::expect_true(maxc_ds202_L48$maxc == 10)
  
  testthat::expect_equal(length(setdiff(maxc_Bill_9002$maxc, 
                                        c(maxc_ds101_L50$maxc,
                                          max(maxc_ds102_L50$maxc,
                                              maxc_ds102_L52$maxc))
                                )),
                         0)
  
  testthat::expect_identical(maxc_Ben_9002$maxc, maxc_ds202_L48$maxc)
})

test_that('drug_spec = c("BARCODE", "DRUG_ID_lib") returns a different set of x values for L50 ds101 and L50 ds102',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = c("BARCODE", "DRUG_ID_lib"), 
                                    cell_line_spec = "MASTER_CELL_ID",
                                    conc_col = "CONC") 
  x_L50_101 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  x_L50_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  testthat::expect_false(isTRUE(all.equal(sort(x_L50_101$x), sort(x_L50_102$x))))
})

test_that('drug_spec = c("BARCODE", "DRUG_ID_lib") returns the same set of x values for L50 ds101 and L52 ds102',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = c("BARCODE", "DRUG_ID_lib"), 
                                    cell_line_spec = "MASTER_CELL_ID",
                                    conc_col = "CONC") 
  x_L50_101 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  x_L52_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L52" & 
                                       set_concs_data$DRUG_ID_lib == 9002, "x"])
  testthat::expect_equal(sort(x_L50_101$x), sort(x_L52_102$x))
})


test_that('drug_spec = c("BARCODE", "DRUG_ID_lib") returns a different set of x values for L50 ds102 and L52 ds102',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = c("BARCODE", "DRUG_ID_lib"), 
                                    cell_line_spec = "MASTER_CELL_ID",
                                    conc_col = "CONC") 
  x_L50_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L50" & 
                                       set_concs_data$DRUG_ID_lib == 9002, 
                                     "x"])
  x_L52_102 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                       set_concs_data$lib_drug == "L52" & 
                                       set_concs_data$DRUG_ID_lib == 9002, 
                                     "x"])
  # testthat::expect_equal(sort(x_L50_102$x), sort(x_L52_102$x))
  testthat::expect_false(isTRUE(all.equal(sort(x_L50_102$x), sort(x_L52_102$x))))
})


# test_that("group_conc_ranges=T returns different maxc if dividing cell lines by expansion_id",{
#   set_concs_data <- setGroupsForNlme(norm_data, group_conc_ranges = T, cell_line_spec = "MASTER_CELL_ID")
#   number_of_maxc_master_id <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002 &
#                                                           set_concs_data$MASTER_CELL_ID == 111, "maxc"]))
#   testthat::expect_true(number_of_maxc_master_id == 1)
#   
#   set_concs_data <- setGroupsForNlme(norm_data, group_conc_ranges = T, cell_line_spec = "CELL_ID")
#   number_of_maxc_cell_id<- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002 &
#                                                         set_concs_data$MASTER_CELL_ID == 111, "maxc"]))
#   testthat::expect_true(number_of_maxc_cell_id == 2)
# })


# DRUG_ID 9002 with maxc of:  20 (L50, ds101); 2 & 5 (L50 & L52, ds102); 10 (L48, ds202).
# Bill (CELL_ID = 6) sees ds101 and Bill (CELL_ID  = 211) 102; Ben sees ds202

test_that('drug_spec = "DRUG_ID_lib" and cell_line_spec = "CELL_ID" returns 2 maxc values for cell line Bill and 1 for Ben',{
  set_concs_data <- setGroupsForNlme(norm_data,  
                                    drug_spec = c("DRUG_ID_lib"), 
                                    cell_line_spec = "CELL_ID",
                                    conc_col = "CONC"
  )
  number_of_maxc <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, "maxc"]))
  maxc_ds202_L48 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 202 &
                                            set_concs_data$lib_drug == "L48" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_ds101_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 101 &
                                            set_concs_data$lib_drug == "L50" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_ds102_L50 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                            set_concs_data$lib_drug == "L50" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_ds102_L52 <- unique(set_concs_data[set_concs_data$DRUGSET_ID == 102 &
                                            set_concs_data$lib_drug == "L52" & 
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  
  maxc_Bill_6_9002 <- unique(set_concs_data[set_concs_data$CELL_LINE_NAME == "Bill" & 
                                            set_concs_data$CELL_ID == 6 &
                                            set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_Bill_211_9002 <- unique(set_concs_data[set_concs_data$CELL_LINE_NAME == "Bill" & 
                                             set_concs_data$CELL_ID == 211 &
                                             set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  maxc_Ben_9002 <- unique(set_concs_data[set_concs_data$CELL_LINE_NAME == "Ben" & 
                                           set_concs_data$DRUG_ID_lib == 9002, "maxc"])
  
  
  unique_specifiers <- nrow(unique(set_concs_data[set_concs_data$DRUG_ID_lib == 9002, 
                                                  c("CELL_ID", "DRUG_ID_lib")]))
  
  # Expect there to be different maxc per cell line
  expect_true(number_of_maxc == 3)
  expect_true(number_of_maxc == unique_specifiers)
  
  testthat::expect_true(maxc_Bill_6_9002$maxc == 20)
  testthat::expect_true(maxc_Bill_211_9002$maxc == 5)
  testthat::expect_true(maxc_Ben_9002$maxc == 10)

  testthat::expect_true(maxc_ds101_L50$maxc == 20)
  testthat::expect_true(maxc_ds102_L50$maxc == 5)
  testthat::expect_true(maxc_ds102_L52$maxc == 5)
  testthat::expect_true(maxc_ds202_L48$maxc == 10)
  
  testthat::expect_identical(maxc_Ben_9002$maxc, maxc_ds202_L48$maxc)
})

test_that("If successful, setGroupsForNlme() returns an informative message",{
  testthat::expect_message(
    setGroupsForNlme(norm_data),
    "^The maximum concentration for a dose response will be calculated for each grouping of MASTER_CELL_ID and DRUG_ID_lib.", perl = T)
})

test_that("With an ill formed drug_spec, setGroupsForNlme() returns an informative message",{
  testthat::expect_error(
    setGroupsForNlme(norm_data, drug_spec =c("tweedle", "dee")),
    "^Your normalized data does not contain the columns specified to make the drug column.", perl = T)
})





