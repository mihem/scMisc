################################################################################
# average expression
################################################################################

#' @title wrapper function for Seurats AverageExpression
#' @description create average expression matrix with Seurat based on `markers.csv`
#' @param par column name in markers.csv
#' @param object Seurat object
#' @param assay which assay to use
#' @param slot which slot to use
#' @param ortho convert to orthologues? Allowed values: `none`, `mouse2human` or `human2mouse`
#' @return matrix with genes as rows and identity clases as columns
#' @examples
#' library(Seurat)
#' markers <- data.frame(Bc = c("CD19", "MS4A1", "CD79B"))
#' write.csv(markers, "markers.csv")
#' avgExp("Bc", object = pbmc_small, assay = "RNA", slot = "data")
#' unlink("markers.csv")
#' @export

avgExp <- function(par, object, assay, slot, ortho = "none") {
  if (!file.exists("markers.csv")) {
    stop("Please make sure that markers.csv file exists")
  }
  if (!inherits(object, "Seurat")) {
    stop("Object must be a Seurat object")
  }
  markers <- readr::read_csv("markers.csv") |>
    as.list(markers) |>
    lapply(function(x) x[!is.na(x)])
  genes <- markers[[par]]
  if (is.null(genes)) {
    stop("No genes were found. Make sure that `par` exists in `markers.csv`")
  }

  if (!(ortho %in% c("none", "mouse2human", "human2mouse"))) {
    stop("ortho must take values: `none`, `mouse2human` or `human2mouse`")
  }
  if (ortho == "mouse2human") {
    genes <- homologene::mouse2human(
      genes,
      db = homologene::homologeneData2
    )$humanGene
  }
  if (ortho == "human2mouse") {
    genes <- homologene::human2mouse(
      genes,
      db = homologene::homologeneData2
    )$mouseGene
  }
  if (ortho == "none") {
    genes <- genes
  }

  Seurat::AverageExpression(
    object,
    features = genes,
    slot = slot,
    assay = assay
  )[[assay]]
}

################################################################################
# find markers presto
################################################################################

#' @title wrapper function for presto find markers
#' @description find significant DE genes using presto
#' @param ident1 cell population 1
#' @param ident2 cell population 2, if NULL use all (default: NULL)
#' @param object Seurat object
#' @param only_pos only return positive markers? (default: FALSE)
#' @param min_pct minimum fraction of cells in either two of the populations (default: 0.1)
#' @param logfc_threshold minimum x-fold difference (default: 0.25)
#' @param assay which assay to use in DE testing (e.g. RNA or SCT)
#' @return data frame with significant DE genes arranged by log2FC
#' @examples
#' library(Seurat)
#' findMarkersPresto(ident1 = "0", ident2 = "1", object = pbmc_small, assay = "RNA")
#' @export
#' @importFrom dplyr desc

findMarkersPresto <- function(
  ident1,
  ident2 = NULL,
  object,
  only_pos = FALSE,
  min_pct = 0.1,
  logfc_threshold = 0.25,
  assay = assay
) {
  if (!inherits(object, "Seurat")) {
    stop("Object must be a Seurat object")
  }
  if (is.null(assay)) {
    stop("Please provide assay information")
  }
  result <- Seurat::FindMarkers(
    object,
    ident.1 = ident1,
    ident.2 = ident2,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    only.pos = only_pos,
    assay = assay
  ) |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(p_val_adj < 0.05) |>
    dplyr::relocate(gene, avg_log2FC, p_val, p_val_adj) |>
    dplyr::arrange(desc(avg_log2FC))
  return(result)
}

################################################################################
# table abundance
################################################################################

#' @title wrapper function to calculate absolute and relative abundance
#' @description wrapper function to calculate absolute and relative abundance
#' @param object Seurat object
#' @param row_var variable in meta data that will represent the rows
#' @param col_var variable in meta data that will represent the columns
#' @param dir_output target directory to save the results (default: .)
#' @examples
#' library(Seurat)
#' abundanceTbl(pbmc_small, row_var = "groups", col_var = "letter.idents")
#' unlink("abundance_tbl_pbmc_small_letter.idents.xlsx")
#' @export

abundanceTbl <- function(object, row_var, col_var, dir_output = ".") {
  if (!inherits(object, "Seurat")) {
    stop("Object must be a Seurat object")
  }
  object_parse <- deparse(substitute(object))
  result_abs <- as.data.frame.matrix(table(
    object@meta.data[[row_var]],
    object@meta.data[[col_var]]
  )) |>
    tibble::rownames_to_column("cell")

  result_pct <- result_abs |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), function(x) {
      x / sum(x) * 100
    })) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), function(x) {
      round(x, 2)
    }))

  file_path <- file.path(
    dir_output,
    glue::glue("abundance_tbl_{object_parse}_{col_var}.xlsx")
  )
  writexl::write_xlsx(
    list("absolute" = result_abs, "percentage" = result_pct),
    file_path
  )
}

################################################################################
# enrichr
################################################################################
#' @title enrichment analysis
#' @description wrapper function to perform enrichment analysis with enrichr and save results
#' @param filename name of deg file file (should be .xlsx) without extension
#' @param dir_input path to the directory where the deg file is located (default: .)
#' @param dir_output path to the directory where the results will be saved (default: .)
#' @param dbs name of the enrichr libraries
#' @param fc_thresh log fc threshold (default 1)
#' @param p_thresh p value threshold (default 0.001)
#' @param sheet sheet name in excel file
#' @param remove_rp_mt remove ribosomal and mitochondrial genes? (boolean value)
#' @return save enrichrment analysis in `dir_output` in two excel sheets (one for positive and one for negative)
#' @examples
#' library(Seurat)
#' library(enrichR)
#' set.seed(123)
#' deg_data <- data.frame(
#'   gene = rownames(pbmc_small),
#'   avg_log2FC = rnorm(length(rownames(pbmc_small)), mean = 0, sd = 3),
#'   p_val_adj = runif(length(rownames(pbmc_small)), min = 0, max = 0.01)
#' )
#'
#' writexl::write_xlsx(deg_data, path = "test_de.xlsx")
#' enrichrRun(
#'   sheet = "Sheet1",
#'   dir_input = ".",
#'   filename = "test_de",
#'   dbs = c("KEGG_2019_Human"),
#'   fc_thresh = 1,
#'   p_thresh = 0.001,
#'   remove_rp_mt = FALSE
#' )
#' unlink("test_de.xlsx")
#' unlink("enrichr_test_de_neg_Sheet1.xlsx")
#' unlink("enrichr_test_de_pos_Sheet1.xlsx")
#' @export

enrichrRun <- function(
  sheet,
  dir_input = ".",
  filename,
  dir_output = ".",
  dbs,
  fc_thresh = 1,
  p_thresh = 0.001,
  remove_rp_mt
) {
  input <- readxl::read_excel(
    file.path(dir_input, glue::glue("{filename}.xlsx")),
    sheet = sheet
  )
  if (remove_rp_mt == TRUE) {
    input <- dplyr::filter(input, !grepl(x = gene, pattern = "(MT-)|(^RP)"))
  }
  input_pos <- input |>
    dplyr::filter(avg_log2FC > fc_thresh) |>
    dplyr::filter(p_val_adj < p_thresh)
  input_neg <- input |>
    dplyr::filter(avg_log2FC < -fc_thresh) |>
    dplyr::filter(p_val_adj < p_thresh)
  input_merge <- list(pos = input_pos, neg = input_neg)
  result <- list()
  for (i in names(input_merge)) {
    if (dim(input_merge[[i]])[[1]] != 0) {
      result[[i]] <- enrichR::enrichr(input_merge[[i]]$gene, dbs)
      for (j in seq_along(result[[i]])) {
        result[[i]][[j]][c("Old.Adjusted.P.value", "Old.P.value")] <- NULL
      } # remove old p value
      names(result[[i]])[
        names(result[[i]]) == "TF_Perturbations_Followed_by_Expression"
      ] <- "TF_Pertubations"
      names(result[[i]])[
        names(result[[i]]) == "Enrichr_Submissions_TF-Gene_Coocurrence"
      ] <- "Enrichr_Submissions_TF"
      writexl::write_xlsx(
        result[[i]],
        file.path(dir_output, glue::glue("enrichr_{filename}_{i}_{sheet}.xlsx"))
      )
    }
  }
}


################################################################################
# helper function for permFDP calculation
################################################################################

#' @title Calculate permFDP-adjusted threshold for propeller results
#' @description Internal helper function to calculate permFDP-adjusted significance threshold for propeller differential abundance testing
#' @param prop_data A dataframe containing propeller results with P.Value column
#' @param seu_obj A Seurat object subset to the two conditions being compared
#' @param sample_col A character representing the name of the sample column
#' @param cluster_col A character representing the name of the cluster column
#' @param meta_col A character representing the name of the condition column
#' @param lookup A dataframe containing sample to condition mapping
#' @param cond1 A character representing the first condition
#' @param cond2 A character representing the second condition
#' @param fdr_threshold The FDR threshold to use for permFDP adjustment
#' @param n_perms Number of permutations for permFDP adjustment (default: 1000)
#' @return A dataframe with added permFDP_sig and permFDP_threshold columns
#' @keywords internal

calculatePropellerPermFDP <- function(
  prop_data,
  seu_obj,
  sample_col,
  cluster_col,
  meta_col,
  lookup,
  cond1,
  cond2,
  fdr_threshold,
  n_perms = 1000
) {
  # Filter to only the two conditions being compared
  seu_subset <- seu_obj[, seu_obj@meta.data[[meta_col]] %in% c(cond1, cond2)]

  # Create proportion matrix manually (clusters x samples)
  # Calculate proportions from cell counts per cluster per sample
  prop_counts <- seu_subset@meta.data |>
    dplyr::group_by(.data[[sample_col]], .data[[cluster_col]]) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::select(.data[[sample_col]], .data[[cluster_col]], prop) |>
    tidyr::pivot_wider(
      names_from = .data[[sample_col]],
      values_from = prop,
      values_fill = 0
    )

  # Convert to matrix with clusters as rows
  prop_matrix <- as.matrix(prop_counts[, -1])
  rownames(prop_matrix) <- prop_counts[[cluster_col]]

  # Filter to clusters present in propeller results
  prop_matrix <- prop_matrix[prop_data$cluster, , drop = FALSE]

  # Remove any rows with NAs
  complete_rows <- complete.cases(prop_matrix)
  n_removed <- sum(!complete_rows)

  if (n_removed > 0) {
    removed_clusters <- prop_data$cluster[!complete_rows]
    warning(
      paste0(
        "Removed ",
        n_removed,
        " cluster(s) with incomplete data for permFDP calculation: ",
        paste(removed_clusters, collapse = ", ")
      )
    )
  }

  prop_matrix <- prop_matrix[complete_rows, , drop = FALSE]
  prop_data_filtered <- prop_data[complete_rows, ]

  if (nrow(prop_matrix) == 0) {
    stop("No complete cases for permFDP calculation. Check your data.")
  }

  # Create group vector from lookup table
  sample_meta <- lookup |>
    dplyr::filter(.data[[sample_col]] %in% colnames(prop_matrix)) |>
    dplyr::arrange(match(.data[[sample_col]], colnames(prop_matrix)))

  group_vector <- ifelse(sample_meta[[meta_col]] == cond1, 1, 2)

  # Validate dimensions match
  if (length(group_vector) != ncol(prop_matrix)) {
    stop(
      paste0(
        "Dimension mismatch: group_vector has ",
        length(group_vector),
        " elements but prop_matrix has ",
        ncol(prop_matrix),
        " columns. "
      )
    )
  }

  # Calculate permFDP adjusted threshold
  message(
    paste0(
      "Running permFDP with ",
      n_perms,
      " permutations for ",
      cond1,
      " vs ",
      cond2,
      " (",
      nrow(prop_matrix),
      " clusters, ",
      ncol(prop_matrix),
      " samples)..."
    )
  )

  corrected_threshold <- permFDP::permFDP.adjust.threshold(
    pVals = prop_data_filtered$P.Value,
    threshold = fdr_threshold,
    myDesign = group_vector,
    intOnly = as.data.frame(prop_matrix),
    nPerms = n_perms
  )

  message(paste0("  -> Corrected threshold: ", round(corrected_threshold, 4)))

  # Add permFDP-adjusted significance to the full data
  prop_data_out <- prop_data |>
    dplyr::mutate(
      permFDP_sig = P.Value < corrected_threshold,
      permFDP_threshold = corrected_threshold
    )

  return(prop_data_out)
}

################################################################################
# propellerCalc
################################################################################

#' @title propellerCalc
#' @description perform a propeller t-test on a Seurat object
#' @param seu_obj1 A seurat object
#' @param condition1 A character representing the first condition
#' @param condition2 A character representing the second condition
#' @param cluster_col A character representing the name of the column which contains cluster information
#' @param meta_col A character representing the name of the column which contains meta information, should be in Seurat object and in lookup table
#' @param lookup A dataframe that contains sample information
#' @param sample_col A character representing the name of the column in seu_obj1 which contains sample information
#' @param formula  A linear model that should be used for the design matrix
#' @param min_cells A numeric value indicating the minimum number of cells in a cluster that should be included in the analysis
#' @param adjustment_method A character indicating the FDR adjustment method: "BH" (default) or "permFDP"
#' @param fdr_threshold The FDR threshold to use for permFDP adjustment (default: NULL)
#' @param n_perms Number of permutations for permFDP adjustment (default: 1000)
#' @return A dataframe containing the results of the t-test
#' @examples
#' set.seed(123)
#' library(Seurat)
#' pbmc_small$condition <- factor(sample(c("diseaseA", "diseaseB"), nrow(pbmc_small), replace = TRUE))
#' pbmc_small$cluster <- Idents(pbmc_small)
#' pbmc_small$patient <- rep(paste0("P", 0:9), each = 8, length.out = ncol(pbmc_small))
#' lookup <- data.frame(
#'   patient = paste0("P", 0:9),
#'   condition = sample(c("diseaseA", "diseaseB"), 10, replace = TRUE)
#' )
#' # Standard BH adjustment
#' propellerCalc(
#'   seu_obj1 = pbmc_small,
#'   condition1 = "diseaseA",
#'   condition2 = "diseaseB",
#'   cluster_col = "cluster",
#'   meta_col = "condition",
#'   lookup = lookup,
#'   sample_col = "patient",
#'   formula = "~ 0 + condition",
#'   min_cells = 30
#' )
#'
#' # With permFDP adjustment (use lower min_cells to have multiple clusters)
#' propellerCalc(
#'   seu_obj1 = pbmc_small,
#'   condition1 = "diseaseA",
#'   condition2 = "diseaseB",
#'   cluster_col = "cluster",
#'   meta_col = "condition",
#'   lookup = lookup,
#'   sample_col = "patient",
#'   formula = "~ 0 + condition",
#'   min_cells = 30,
#'   adjustment_method = "permFDP",
#'   fdr_threshold = 0.1,
#'   n_perms = 100
#' )
#' @export
#' @importFrom stats as.formula model.matrix

propellerCalc <- function(
  seu_obj1,
  condition1,
  condition2,
  cluster_col,
  meta_col,
  lookup,
  sample_col,
  formula,
  min_cells = 30,
  adjustment_method = "BH",
  fdr_threshold = NULL,
  n_perms = 1000
) {
  # filter cells of condition1 and condition2
  seu_obj2 <- seu_obj1[,
    seu_obj1@meta.data[[meta_col]] %in% c(condition1, condition2)
  ]

  # filter clusters with less than min_cells
  cl_interest <-
    as.data.frame.matrix(table(
      seu_obj2@meta.data[[cluster_col]],
      seu_obj2@meta.data[[meta_col]]
    )) |>
    tibble::rownames_to_column("cluster") |>
    dplyr::mutate(group_sum = .data[[condition1]] + .data[[condition2]]) |>
    dplyr::filter(group_sum > min_cells) |>
    dplyr::pull(cluster)

  # transform proportions
  props <- speckle::getTransformedProps(
    clusters = seu_obj2@meta.data[[cluster_col]],
    sample = seu_obj2@meta.data[[sample_col]],
    transform = "logit"
  )

  # create lookup table
  meta_lookup <-
    tibble::tibble(
      !!sample_col := colnames(props$TransformedProps)
    ) |>
    dplyr::left_join(lookup, by = sample_col) |>
    dplyr::distinct(.data[[sample_col]], .keep_all = TRUE)

  my_design <- model.matrix(as.formula(formula), data = meta_lookup)

  my_contrasts <- glue::glue("{meta_col}{condition1}-{meta_col}{condition2}")
  my_args <- list(my_contrasts, levels = my_design)
  my_contr <- do.call(limma::makeContrasts, my_args)

  # calculate propeller t test
  propeller_result <-
    speckle::propeller.ttest(
      prop.list = props,
      design = my_design,
      contrast = my_contr,
      robust = TRUE,
      trend = FALSE,
      sort = TRUE
    ) |>
    tibble::rownames_to_column("cluster") |>
    dplyr::filter(cluster %in% cl_interest) |>
    dplyr::mutate(log2ratio = log2(PropRatio)) |>
    dplyr::mutate(FDR_log = -log10(FDR)) |>
    tibble::tibble()

  # Apply permFDP adjustment if requested
  if (adjustment_method == "permFDP") {
    if (is.null(fdr_threshold)) {
      stop(
        "fdr_threshold must be specified when using adjustment_method = 'permFDP'"
      )
    }
    propeller_result <- calculatePropellerPermFDP(
      prop_data = propeller_result,
      seu_obj = seu_obj2,
      sample_col = sample_col,
      cluster_col = cluster_col,
      meta_col = meta_col,
      lookup = lookup,
      cond1 = condition1,
      cond2 = condition2,
      fdr_threshold = fdr_threshold,
      n_perms = n_perms
    )
  }

  return(propeller_result)
}

################################################################################
# read cellbender data, modified from scCustomize
################################################################################

#' Load CellBender h5 matrices
#'
#' Extract sparse matrix with corrected counts from CellBender h5 output file.
#'
#' @param file_name Path to h5 file.
#' @param use.names Label row names with feature names rather than ID numbers (default TRUE).
#' @param unique.features Make feature names unique (default TRUE).
#'
#' @return sparse matrix
#'
#' @references Code used in function has been modified from `Seurat::Read10X_h5` function of
#' Seurat package \url{https://github.com/satijalab/seurat} (License: GPL-3).
#'
#' @import Matrix
#'
#' @export
#'
#' @examples
#' file_name <- system.file("extdata", "raw_feature_bc_matrix_filtered.h5", package = "scMisc")
#' mat <- ReadCellBender_h5(file_name)
#' @importFrom methods new

ReadCellBender_h5 <- function(
  file_name,
  use.names = TRUE,
  unique.features = TRUE
) {
  # Check hdf5r installed
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  # Check file
  if (!file.exists(file_name)) {
    stop("File not found")
  }

  if (use.names) {
    feature_slot <- "features/name"
  } else {
    feature_slot <- "features/id"
  }

  # Read file
  infile <- hdf5r::H5File$new(filename = file_name, mode = "r")

  counts <- infile[["matrix/data"]]
  indices <- infile[["matrix/indices"]]
  indptr <- infile[["matrix/indptr"]]
  shp <- infile[["matrix/shape"]]
  features <- infile[[paste0("matrix/", feature_slot)]][]
  barcodes <- infile[["matrix/barcodes"]]

  sparse.mat <- Matrix::sparseMatrix(
    i = indices[] + 1,
    p = indptr[],
    x = as.numeric(x = counts[]),
    dims = shp[],
    repr = "T"
  )

  if (unique.features) {
    features <- make.unique(names = features)
  }

  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- Seurat::as.sparse(x = sparse.mat)

  infile$close_all()

  return(sparse.mat)
}
