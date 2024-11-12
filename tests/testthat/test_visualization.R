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
    fPlot(path = "markers.csv", object = pbmc_small, par = "B", reduction = "tsne", order = TRUE, dir_output = ".")
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
    # Cleanup: Remove the generated file
    unlink("markers.csv")
    unlink("fp_pbmc_small_B.png")
})

test_that("fPlotCustom works as expected", {
    library(Seurat)

    markers <- data.frame(cell_source = c("B", "B"), gene = c("MS4A1", "CD79A"))
    # run the function
    fPlotCustom(
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
    # Cleanup: Remove the generated file
    unlink("fp_pbmc_small_B.png")
})

test_that("dotPlot works as expected", {
    library(Seurat)
    markers <- data.frame(T = c("CD3E", "CD3D", "IL7R", "TRBC"))
    write.csv(markers, "markers.csv")

    # run the function
    suppressWarnings(
        dotPlot(
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
