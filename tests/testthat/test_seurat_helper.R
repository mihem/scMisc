library(testthat)
library(Seurat)
library(dplyr)
library(readxl)

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
  expected_result <- tibble::tibble(cell = c("g1", "g2"), A = c(30, 23), B = c(14,13))
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
