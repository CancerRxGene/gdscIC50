context("Test calcTagMean")

test_that("calcTagMean", {
  test_data <- rbind(
    data.frame(SCAN_ID = 101, TAG = "NC1", INTENSITY = c(4648, 9846, 1563, 1519)),
    data.frame(SCAN_ID = 101, TAG = "NC2", INTENSITY = c(4648, 9846, 1563, 1519)),
    data.frame(SCAN_ID = 102, TAG = "NC2", INTENSITY = c(1672, 4366, 1651, 6527)),
    data.frame(SCAN_ID = 103, TAG = "NC1", INTENSITY = c(6175, 1381, 4367, 6345)),
    data.frame(SCAN_ID = 103, TAG = "NC2", INTENSITY = c(NA, NA, NA, NA))
  )
  # No tag at all of "foo"
  expect_error(calcTagMean(test_data, "foo"), 
               "is not present for some or all of the SCAN_IDs in your data")
  # No NC1 for scan 102
  expect_error(calcTagMean(test_data, "NC1"), 
               "is not present for some or all of the SCAN_IDs in your data")
  # Mean of NC2 for scan 103 is NaN
  expect_error(calcTagMean(test_data, "NC2"),
               "has a mean of NaN for some or all of the SCAN_IDs in your data")
})

context("Test normalize data")

test_that("normalizeData returns a data frame for single agent only data", {
  # data("gdsc_example")
  test_data <- readRDS("norm_single_agents_test.rds")
  
  norm_data <- normalizeData(test_data, 
                             neg_control = "NC-1", 
                             pos_control = "PC1-D1-S")
  expect_true(is.data.frame(norm_data))
  
  expected_names <- c("SCAN_ID",
                      "BARCODE",
                      "RESEARCH_PROJECT",
                      "DATE_CREATED",
                      "DRUGSET_ID",  
                      "CELL_LINE_NAME",
                      "CELL_ID",
                      "MASTER_CELL_ID",
                      "COSMIC_ID",
                      "POSITION",
                      "DRUG_ID_lib",
                      "CONC",
                      "INTENSITY",
                      "lib_drug",
                      "dose",
                      "treatment",
                      "DRUG_ID_anch", 
                      "CONC_anch", 
                      "anchor",
                      "NC",
                      "PC",
                      "normalized_intensity",
                      "norm_neg_pos",
                      "time_stamp",
                      "sw_version")
  expect_named(norm_data, expected_names)
  
})

test_that("normalizeData returns a data frame with combination treatment 
          data with the right names", {
            test_data <- readRDS("combo_norm_data_test.rds")
            norm_data <- normalizeData(test_data, 
                                       neg_control = "NC-1", 
                                       pos_control = "B")
            expect_true(is.data.frame(norm_data))
            
            expected_names <- c("SCAN_ID",
                                "BARCODE",
                                "RESEARCH_PROJECT",
                                "DATE_CREATED",
                                "DRUGSET_ID",  
                                "CELL_LINE_NAME",
                                "CELL_ID",
                                "MASTER_CELL_ID",
                                "COSMIC_ID",
                                "POSITION",
                                "DRUG_ID_lib",
                                "CONC",
                                "INTENSITY",
                                "lib_drug",
                                "dose",
                                "treatment",
                                "DRUG_ID_anch",
                                "CONC_anch",
                                "anchor",
                                "NC",
                                "PC",
                                "normalized_intensity",
                                "norm_neg_pos",
                                "time_stamp",
                                "sw_version")
            expect_named(norm_data, expected_names)
            
          })

test_that("normalizeData maintains data mappings to barcode", {
  test_data <- readRDS("norm_single_agents_test.rds")
  norm_data <- normalizeData(test_data)
  expect_identical(unique(test_data[test_data$BARCODE == 1234, "RESEARCH_PROJECT"]), 
                   unique(norm_data[norm_data$BARCODE == 1234, "RESEARCH_PROJECT"]))
  expect_identical(unique(test_data[test_data$BARCODE == 1234, "SCAN_ID"]), 
                   unique(norm_data[norm_data$BARCODE == 1234, "SCAN_ID"]))
  expect_identical(unique(test_data[test_data$BARCODE == 1234, "COSMIC_ID"]), 
                   unique(norm_data[norm_data$BARCODE == 1234, "COSMIC_ID"]))
  expect_identical(unique(test_data[test_data$BARCODE == 1234, "CELL_LINE_NAME"]), 
                   unique(norm_data[norm_data$BARCODE == 1234, "CELL_LINE_NAME"]))
  expect_identical(unique(test_data[test_data$BARCODE == 1234, "CELL_ID"]), 
                   unique(norm_data[norm_data$BARCODE == 1234, "CELL_ID"]))
  expect_identical(unique(test_data[test_data$BARCODE == 1234, "MASTER_CELL_ID"]), 
                   unique(norm_data[norm_data$BARCODE == 1234, "MASTER_CELL_ID"]))
  expect_identical(unique(test_data[test_data$BARCODE == 1234, "DRUGSET_ID"]), 
                   unique(norm_data[norm_data$BARCODE == 1234, "DRUGSET_ID"]))
  
  expect_identical(unique(test_data[test_data$BARCODE == 5678, "RESEARCH_PROJECT"]), 
                   unique(norm_data[norm_data$BARCODE == 5678, "RESEARCH_PROJECT"]))
  expect_identical(unique(test_data[test_data$BARCODE == 5678, "SCAN_ID"]), 
                   unique(norm_data[norm_data$BARCODE == 5678, "SCAN_ID"]))
  expect_identical(unique(test_data[test_data$BARCODE == 5678, "COSMIC_ID"]), 
                   unique(norm_data[norm_data$BARCODE == 5678, "COSMIC_ID"]))
  expect_identical(unique(test_data[test_data$BARCODE == 5678, "CELL_LINE_NAME"]), 
                   unique(norm_data[norm_data$BARCODE == 5678, "CELL_LINE_NAME"]))
  expect_identical(unique(test_data[test_data$BARCODE == 5678, "CELL_ID"]), 
                   unique(norm_data[norm_data$BARCODE == 5678, "CELL_ID"]))
  expect_identical(unique(test_data[test_data$BARCODE == 5678, "MASTER_CELL_ID"]), 
                   unique(norm_data[norm_data$BARCODE == 5678, "MASTER_CELL_ID"]))
  expect_identical(unique(test_data[test_data$BARCODE == 5678, "DRUGSET_ID"]), 
                   unique(norm_data[norm_data$BARCODE == 5678, "DRUGSET_ID"]))
  
  intensity_test <- merge(
    test_data[, c("SCAN_ID", "POSITION", "INTENSITY")],
    norm_data[, c("SCAN_ID", "POSITION", "INTENSITY")],
    by = c("SCAN_ID", "POSITION"))
  expect_identical(intensity_test$INTENSITY.x, intensity_test$INTENSITY.y)
  
})

test_that("normalizeData calculates correct mean values for negative controls", {
  test_data <- readRDS("norm_single_agents_test.rds")
  
  norm_data <- normalizeData(test_data, 
                             neg_control = "NC-1", 
                             pos_control = "PC1-D1-S")
  
  expect_length(unique(norm_data[norm_data$SCAN_ID == 1199, "NC"]), 1)
  
  expect_equal(unique(norm_data[norm_data$SCAN_ID == 1199, "NC"]),
               mean(test_data[test_data$SCAN_ID == 1199 & test_data$TAG == 'NC-1',
                                         "INTENSITY"]
                    )
               )
  expect_length(unique(norm_data[norm_data$SCAN_ID == 1299, "NC"]), 1)
  
  expect_equal(unique(norm_data[norm_data$SCAN_ID == 1299, "NC"]),
               mean(test_data[test_data$SCAN_ID == 1299 & test_data$TAG == 'NC-1',
                              "INTENSITY"]
                    )
               )
  })

test_that("normalizeData calculates correct mean values for positive controls", {
  
  test_data <- readRDS("norm_single_agents_test.rds")
  
  # Test PCs with different values
  norm_data <- normalizeData(test_data, 
                             neg_control = "NC-1", 
                             pos_control = "B")
  
  expect_length(unique(norm_data[norm_data$SCAN_ID == 1199, "PC"]), 1)
  
  expect_equal(unique(norm_data[norm_data$SCAN_ID == 1199, "PC"]),
               mean(test_data[test_data$SCAN_ID == 1199 & test_data$TAG == 'B',
                              "INTENSITY"]
               )
  )
  expect_length(unique(norm_data[norm_data$SCAN_ID == 1299, "PC"]), 1)
  
  expect_equal(unique(norm_data[norm_data$SCAN_ID == 1299, "PC"]),
               mean(test_data[test_data$SCAN_ID == 1299 & test_data$TAG == 'B',
                              "INTENSITY"]
               )
  )
})

test_that("normalizeData calculates the normalized intensity correctly",{
  test_data <- readRDS("norm_single_agents_test.rds")
  
  norm_data <- normalizeData(test_data, 
                             neg_control = "NC-1", 
                             pos_control = "PC1-D1-S",
                             trim = "F")
  
  normalization_test <- (norm_data$INTENSITY - norm_data$PC) / (norm_data$NC - norm_data$PC)
  
  expect_identical(normalization_test, norm_data$normalized_intensity)
})

test_that("normalizeData trims the normalized intensity correctly",{
  test_data <- readRDS("norm_single_agents_test.rds")
  
  norm_data <- normalizeData(test_data, 
                             neg_control = "NC-1", 
                             pos_control = "PC1-D1-S",
                             trim = "T")
  
  normalization_test <- (norm_data$INTENSITY - norm_data$PC) / (norm_data$NC - norm_data$PC)
  
  normalization_test[normalization_test > 1 ] <- 1
  normalization_test[normalization_test < 0 ] <- 0
  
  expect_identical(normalization_test, norm_data$normalized_intensity)
})

test_that("normalizeData outputs library and anchor data for C treatments", {
  test_data <- readRDS("combo_norm_data_test.rds")
  norm_data <- normalizeData(test_data, 
                             neg_control = "NC-1", 
                             pos_control = "B")
  
  c_treatments <- norm_data[norm_data$treatment == 'C',]
  
  expect_false(any(is.na(c_treatments$DRUG_ID_lib)))
  expect_false(any(is.na(c_treatments$CONC)))
  expect_false(any(is.na(c_treatments$lib_drug)))
  expect_false(any(is.na(c_treatments$dose)))
  expect_false(any(is.na(c_treatments$treatment)))
  expect_false(any(is.na(c_treatments$DRUG_ID_anch)))
  expect_false(any(is.na(c_treatments$CONC_anch)))
  expect_false(any(is.na(c_treatments$anchor)))
})

test_that("normalizeComboData points to and gives same results as normalizeData", {
  test_data <- readRDS("combo_norm_data_test.rds")
  norm_data <- normalizeData(test_data, 
                             neg_control = "NC-1", 
                             pos_control = "B")
  
  norm_combo_data <- normalizeComboData(test_data, 
                             neg_control = "NC-1", 
                             pos_control = "B")
  
  norm_data$time_stamp <- NULL
  norm_combo_data$time_stamp <- NULL
  
  expect_identical(norm_data, norm_combo_data)
})
