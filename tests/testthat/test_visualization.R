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
    # Cleanup: Remove the generated file
    unlink("markers.csv")
    unlink("fp_pbmc_small_B.png")
})
