library(testthat)
library(Seurat)
library(dplyr)
library(readxl)
library(writexl)
library(enrichR)

test_that("avgExp works correctly", {
  markers <- data.frame(Bc = c("CD19", "MS4A1", "CD79B"))
  write.csv(markers, "markers.csv")
  result <- avgExp("Bc", object = pbmc_small, assay = "RNA", slot = "data")
  expect_true(inherits(result, "dgCMatrix"))
  expect_true(all(c("CD19", "MS4A1", "CD79B") %in% rownames(result)))
  unlink("markers.csv")
})

test_that("findMarkersPresto works correctly", {
  result <- findMarkersPresto(ident1 = "0", ident2 = "1", object = pbmc_small, assay = "RNA")
  expect_true(is.data.frame(result))
  expect_true("gene" %in% colnames(result))
  expect_true("avg_log2FC" %in% colnames(result))
})

test_that("abundanceTbl works correctly", {
  abundanceTbl(pbmc_small, row_var = "groups", col_var = "letter.idents", target_dir = ".")
  file_path <- "abundance_tbl_pbmc_small_letter.idents.xlsx"
  result <- readxl::read_xlsx(file_path)
  expect_true(file.exists(file_path))
  expect_true(is.data.frame(result))
  expected_result <- tibble::tibble(cell = c("g1", "g2"), A = c(30, 23), B = c(14, 13))
  expect_equal(result, expected_result)
  unlink(file_path)
})

test_that("propellerCalc works correctly", {
  pbmc_small$condition <- factor(sample(c("diseaseA", "diseaseB"), nrow(pbmc_small), replace = TRUE))
  pbmc_small$cluster <- Idents(pbmc_small)
  pbmc_small$patient <- rep(paste0("P", 0:9), each = 8, length.out = ncol(pbmc_small))
  lookup <- data.frame(patient = paste0("P", 0:9), condition = sample(c("diseaseA", "diseaseB"), 10, replace = TRUE))

  result <- propellerCalc(
    seu_obj1 = pbmc_small,
    condition1 = "diseaseA",
    condition2 = "diseaseB",
    cluster_col = "cluster",
    meta_col = "condition",
    lookup = lookup,
    sample_col = "patient",
    formula = "~ 0 + condition",
    min_cells = 30
  )
  expect_true(is.data.frame(result))
  expect_true("cluster" %in% colnames(result))
})

test_that("enrichrRun works correctly", {
  # Create a dummy DE file
  set.seed(123)
  deg_data <- data.frame(
    gene = rownames(pbmc_small),
    avg_log2FC = rnorm(length(rownames(pbmc_small)), mean = 0, sd = 3),
    p_val_adj = runif(length(rownames(pbmc_small)), min = 0, max = 0.01)
  )

  writexl::write_xlsx(deg_data, path = "test_de.xlsx")
  enrichrRun(
    sheet = "Sheet1",
    dir_input = ".",
    filename = "test_de",
    dbs = c("KEGG_2019_Human"),
    fc_thresh = 1,
    p_thresh = 0.001,
    remove_rp_mt = FALSE
  )

  result_file_pos <- "enrichr_test_de_pos_Sheet1.xlsx"
  result_file_neg <- "enrichr_test_de_neg_Sheet1.xlsx"
  result_pos <- readxl::read_xlsx(result_file_pos)
  result_neg <- readxl::read_xlsx(result_file_neg)

  # Check if the files exist
  expect_true(file.exists(result_file_pos))
  expect_true(file.exists(result_file_neg))

  # Check if the results are data frames
  expect_true(is.data.frame(result_pos))
  expect_true(is.data.frame(result_neg))

  # Check if the results have the expected columns
  expected_columns <- c("Term", "Overlap", "P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes")
  expect_identical(colnames(result_pos), expected_columns)
  expect_identical(colnames(result_neg), expected_columns)

  # Clean up
  unlink("test_de.xlsx")
  unlink("enrichr_test_de_neg_Sheet1.xlsx")
  unlink("enrichr_test_de_pos_Sheet1.xlsx")
})

test_that("ReadCellBender_h5 works correctly", {
  file_name <- system.file("extdata", "raw_feature_bc_matrix_filtered.h5", package = "scMisc")
  result <- ReadCellBender_h5(file_name)
  expect_true(is(result, "dgCMatrix"))
  expect_true(nrow(result) > 0)
  expect_true(ncol(result) > 0)
})
