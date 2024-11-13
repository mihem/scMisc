test_that("theme_rect works as expected", {
    p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
        geom_point()
    p <- p + theme_rect()
    expect_s3_class(p, "ggplot")
    expect_equal(p$theme$panel.border$colour, "black")
    expect_equal(p$theme$panel.border$linewidth, 1)
    expect_equal(p$theme$panel.border$fill, NA)
    expect_equal(p$theme$aspect.ratio, 1)
})

test_that("fPlot works as expected", {
    library(Seurat)
    markers <- data.frame(B = c("MS4A1", "CD79A"))
    write.csv(markers, "markers.csv")
    # run the function
    plot <- fPlot(path = "markers.csv", object = pbmc_small, par = "B", reduction = "tsne", order = TRUE, dir_output = ".")
    # Test 1: Function creates a file in the correct directory
    expect_true(file.exists("fp_pbmc_small_B.png"))
    # Test 2: Function  creates a file that is not empty
    expect_gt(file.info("fp_pbmc_small_B.png")$size, 0)
    # Test 3: Expect error if object is not a Seurat object
    expect_error(
        fPlot(path = "markers.csv", object = data.frame(a = c(1:3)), par = "B", reduction = "tsne", order = TRUE, dir_output = "."),
        "Object must be a Seurat object"
    )
    # Test 4: Check if the function throws an error if no genes are found
    markers <- data.frame(B = c())
    write.csv(markers, "markers.csv")
    expect_error(
        fPlot(path = "markers.csv", object = pbmc_small, par = "B", reduction = "tsne", order = TRUE, dir_output = "."),
        "No genes were found. Make sure that `par` exists in markers.csv"
    )
    # Test 5: Check plot objects
    expect_s3_class(plot, "ggplot")
    expect_equal(plot[[1]]$label$title, "MS4A1")
    expect_equal(plot[[2]]$label$title, "CD79A")

    # Cleanup: Remove the generated file
    unlink("markers.csv")
    unlink("fp_pbmc_small_B.png")
})

test_that("fPlotCustom works as expected", {
    library(Seurat)

    markers <- data.frame(cell_source = c("B", "B"), gene = c("MS4A1", "CD79A"))
    # run the function
    plot <- fPlotCustom(
        object = pbmc_small,
        markers = markers,
        par = "B",
        reduction = "tsne"
    )
    # Test 1: Function creates a file in the correct directory
    expect_true(file.exists("fp_pbmc_small_B.png"))
    # Test 2: Function  creates a file that is not empty
    expect_gt(file.info("fp_pbmc_small_B.png")$size, 0)
    # Test 3: Expect error if object is not a Seurat object
    expect_error(
        fPlotCustom(
            object = data.frame(a = c(1:3)),
            markers = markers,
            par = "B",
            reduction = "tsne"
        ),
        "Object must be a Seurat object"
    )

    # Test 4: Check plot objects
    expect_s3_class(plot, "ggplot")
    expect_equal(plot[[1]]$label$title, "MS4A1")
    expect_equal(plot[[2]]$label$title, "CD79A")
    # Cleanup: Remove the generated file
    unlink("fp_pbmc_small_B.png")
})

test_that("dotPlot works as expected", {
    library(Seurat)
    markers <- data.frame(T = c("CD3E", "CD3D", "IL7R", "TRBC"))
    write.csv(markers, "markers.csv")

    # run the function
    suppressWarnings(
        plot <- dotPlot(
            path = "markers.csv",
            object = pbmc_small,
            par = "T",
            dot_min = 0.01,
            dir_output = "."
        )
    )

    # Test 1: Function creates a file in the correct directory
    expect_true(file.exists("dp_pbmc_small_T.pdf"))
    # Test 2: Function creates a file that is not empty
    expect_gt(file.info("dp_pbmc_small_T.pdf")$size, 0)
    # Test 3: Expect error if object is not a Seurat object
    expect_error(
        dotPlot(
            path = "markers.csv",
            object = data.frame(a = c(1:3)),
            par = "T",
            dot_min = 0.1,
            dir_output = "."
        ),
        "Object must be a Seurat object"
    )
    # Test 4: Check if the function throws an error if no genes are found
    markers <- data.frame(T = c())
    write.csv(markers, "markers.csv")
    expect_error(
        dotPlot(
            path = "markers.csv",
            object = pbmc_small,
            par = "B",
            dot_min = 0.1,
            dir_output = "."
        ),
        "No genes were found. Make sure that `par` exists in markers.csv"
    )
    # Test 5 Check if mouse2human conversion works
    markers <- data.frame(T = c("Cd3e", "Cd3d", "Il7r"))
    write.csv(markers, "markers.csv")
    suppressWarnings(
        dotPlot(
            path = "markers.csv",
            object = pbmc_small,
            par = "T",
            dot_min = 0.01,
            ortho = "mouse2human"
        )
    )
    expect_gt(file.info("dp_pbmc_small_T.pdf")$size, 0)
    # Test 6: Check plot objects
    expect_s3_class(plot, "ggplot")
    p <- ggplot2::ggplot_build(plot)
    expect_equal(
        p$layout$panel_params[[1]]$x$limits,
        c("CD3E", "CD3D", "IL7R")
    )
    # Cleanup: Remove the generated file
    unlink("markers.csv")
    unlink("dp_pbmc_small_T.pdf")
})

test_that("pHeatmap works as expected", {
    # Create a sample matrix
    matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
    rownames(matrix) <- paste0("Gene", 1:10)
    colnames(matrix) <- paste0("Sample", 1:10)

    # Run the function
    pHeatmap(matrix, scale = "row", dir_output = ".")

    # Test 1: Function creates a file in the correct directory
    expect_true(file.exists("hm_matrix.pdf"))
    # Test 2: Function creates a file that is not empty
    expect_gt(file.info("hm_matrix.pdf")$size, 0)

    # Cleanup: Remove the generated file
    unlink("hm_matrix.pdf")
})

test_that("stackedPlot works as expected", {
    library(Seurat)
    set.seed(123)
    pbmc_small$disease <- sample(c("diseaseA", "diseaseB"), ncol(pbmc_small), replace = TRUE)
    pbmc_small$cluster <- sample(c("Cluster1", "Cluster2"), ncol(pbmc_small), replace = TRUE)

    # Run the function
    plot <- stackedPlot(
        object = pbmc_small,
        x_axis = "disease",
        y_axis = "cluster",
        x_order = c("diseaseA", "diseaseB"),
        y_order = c("Cluster1", "Cluster2"),
        color = c("Cluster1" = "blue", "Cluster2" = "red"),
        width = 10,
        dir_output = "."
    )

    # Test 1: Function creates a file in the correct directory
    expect_true(file.exists("stacked_barplot_pbmc_small_disease.pdf"))
    # Test 2: Function creates a file that is not empty
    expect_gt(file.info("stacked_barplot_pbmc_small_disease.pdf")$size, 0)
    # Test 3: Test plot objects
    expect_s3_class(plot, "ggplot")
    p <- ggplot2::ggplot_build(plot)
    expect_equal(
        p$layout$panel_params[[1]]$x$limits,
        c("diseaseA", "diseaseB")
    )
    # Cleanup: Remove the generated file
    unlink("stacked_barplot_pbmc_small_disease.pdf")
})

test_that("abVolPlot works as expected", {
    library(Seurat)
    set.seed(123)
    pbmc_small$predicted.id <- sample(c("Mono", "Tcell"), ncol(pbmc_small), replace = TRUE)
    pbmc_small$sample <- sample(c("CSF_P01", "CSF_P02", "CSF_P03", "CSF_P04"), ncol(pbmc_small), replace = TRUE)
    lookup <-
        data.frame(
            sample = c("CSF_P01", "CSF_P02", "CSF_P03", "CSF_P04"),
            AIE_type = c(rep("control", 2), rep("CASPR2", 2))
        )
    pbmc_small@meta.data <-
        pbmc_small@meta.data |>
        tibble::rownames_to_column("barcode") |>
        dplyr::left_join(lookup, by = "sample") |>
        tibble::column_to_rownames("barcode")


    # Run the function
    suppressWarnings(
        plot <- abVolPlot(
            object = pbmc_small,
            cluster_idents = "predicted.id",
            sample = "sample",
            cluster_order = c("Mono", "Tcell"),
            group_by = "AIE_type",
            group1 = "CASPR2",
            group2 = "control",
            color = c("Mono" = "blue", "Tcell" = "red"),
            width = 5,
            height = 5,
            dir_output = "."
        )
    )

    # Test 1: Function creates a file in the correct directory
    expect_true(file.exists("volcano_plot_predicted.id_pbmc_small_CASPR2_control.pdf"))
    # Test 2: Function creates a file that is not empty
    expect_gt(file.info("volcano_plot_predicted.id_pbmc_small_CASPR2_control.pdf")$size, 0)
    # Test 3: Test plot objects
    expect_s3_class(plot, "ggplot")
    p <- ggplot2::ggplot_build(plot)
    expect_equal(
        p$layout$panel_params[[1]]$x$limits,
        c(-1, 1)
    )

    # Test 4: Check if the function throws an error if object is not a Seurat object
    expect_error(
        abVolPlot(
            object = data.frame(a = c(1:3)),
            cluster_idents = "predicted.id",
            sample = "sample",
            cluster_order = c("Cluster1", "Cluster2"),
            group_by = "AIE_type",
            group1 = "LGI1",
            group2 = "control",
            color = c("Cluster1" = "blue", "Cluster2" = "red"),
            width = 5,
            height = 5,
            dir_output = "."
        ),
        "Object must be a Seurat object"
    )

    # Test 5: Check if the plot has the expected elements
    expect_true("GeomPoint" %in% sapply(plot$layers, function(x) class(x$geom)[1]))

    # Cleanup: Remove the generated file
    unlink("volcano_plot_predicted.id_pbmc_small_CASPR2_control.pdf")
})

test_that("abBoxPlot works as expected", {
    library(Seurat)
    set.seed(123)
    pbmc_small$cluster <- sample(c("Cluster1", "Cluster2"), ncol(pbmc_small), replace = TRUE)
    pbmc_small$sample <- sample(c("CSF_P01", "CSF_P02", "CSF_P03", "CSF_P04"), ncol(pbmc_small), replace = TRUE)
    lookup <-
        data.frame(
            sample = c("CSF_P01", "CSF_P02", "CSF_P03", "CSF_P04"),
            AIE_type = c(rep("control", 2), rep("CASPR2", 2))
        )
    pbmc_small@meta.data <-
        pbmc_small@meta.data |>
        tibble::rownames_to_column("barcode") |>
        dplyr::left_join(lookup, by = "sample") |>
        tibble::column_to_rownames("barcode")

    # Run the function
    suppressWarnings(
        plot <- abBoxPlot(
            object = pbmc_small,
            cluster_idents = "cluster",
            sample = "sample",
            cluster_order = c("Cluster1", "Cluster2"),
            group_by = "AIE_type",
            group_order = c("control", "CASPR2"),
            color = c("control" = "blue", "CASPR2" = "red"),
            width = 9,
            height = 6,
            paired = FALSE,
            number_of_tests = 3,
            dir_output = "."
        )
    )

    # Test 1: Function creates a file in the correct directory
    expect_true(file.exists("boxplot_cluster_pbmc_small_AIE_type.pdf"))
    # Test 2: Function creates a file that is not empty
    expect_gt(file.info("boxplot_cluster_pbmc_small_AIE_type.pdf")$size, 0)
    # Test 3: Test plot objects
    expect_s3_class(plot, "ggplot")
    # Test 4: Check if the function throws an error if object is not a Seurat object
    expect_error(
        abBoxPlot(
            object = data.frame(a = c(1:3)),
            cluster_idents = "cluster",
            sample = "sample",
            cluster_order = c("Cluster1", "Cluster2"),
            group_by = "AIE_type",
            group_order = c("control", "CASPR2", "LGI1"),
            color = c("control" = "blue", "CASPR2" = "red", "LGI1" = "green"),
            width = 9,
            height = 6,
            paired = FALSE,
            number_of_tests = 3,
            dir_output = "."
        ),
        "Object must be a Seurat object"
    )

    # Cleanup: Remove the generated file
    unlink("boxplot_cluster_pbmc_small_AIE_type.pdf")
})

test_that("compStat works as expected", {
    set.seed(123)
    data <- data.frame(
        sample = c(paste0("CSF_P0", 1:9)),
        type = c(rep("control", 4), rep("AIE", 5)),
        cluster1 = c(runif(4, 0, 1), runif(5, 99, 100)),
        cluster2 = c(runif(4, 0, 70), runif(5, 20, 100))
    )

    # Test 1: Check if the function returns a data frame
    result <- compStat(x_var = c("cluster1", "cluster2"), group = "type", data = data, paired = FALSE)
    expect_s3_class(result, "data.frame")

    # Test 2: Check if the function returns the correct columns
    expect_true(all(c(".y.", "group1", "group2", "p", "p.adj", "p.adj.signif") %in% colnames(result)))

    # Test 3: Check if the function returns significant values
    expect_true(any(result$p.adj < 0.05))
})
