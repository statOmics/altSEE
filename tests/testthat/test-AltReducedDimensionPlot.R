# Tests for AltReducedDimensionPlot
# library(testthat); library(altSEE); source("setup-sce.R"); source("test-AltReducedDimensionPlot.R")

test_that("AltReducedDimensionPlot constructor works", {
    x <- AltReducedDimensionPlot()
    expect_s4_class(x, "AltReducedDimensionPlot")
    expect_identical(x[["Experiment"]], NA_character_)
    expect_identical(x[["Type"]], NA_character_)
    expect_identical(x[["XAxis"]], 1L)
    expect_identical(x[["YAxis"]], 2L)
})

test_that("AltReducedDimensionPlot constructor accepts slot values", {
    x <- AltReducedDimensionPlot(Experiment = "peptides", Type = "PCA", XAxis = 1L, YAxis = 2L)
    expect_identical(x[["Experiment"]], "peptides")
    expect_identical(x[["Type"]], "PCA")
})

test_that("AltReducedDimensionPlot validity rejects non-positive axis indices", {
    expect_error(AltReducedDimensionPlot(XAxis = 0L), "must be a single positive integer")
    expect_error(AltReducedDimensionPlot(YAxis = -1L), "must be a single positive integer")
})

test_that(".fullName and .panelColor work for AltReducedDimensionPlot", {
    x <- AltReducedDimensionPlot()
    expect_identical(iSEE:::.fullName(x), "Alt reduced dimension plot")
    expect_type(iSEE:::.panelColor(x), "character")
})

test_that(".cacheCommonInfo caches altExp names", {
    x   <- AltReducedDimensionPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    se2 <- iSEE:::.cacheCommonInfo(x, se2)  # second call returns early

    cache <- iSEE:::.getCachedCommonInfo(se2, "AltReducedDimensionPlot")
    expect_true("peptides" %in% cache$valid.altExp.names)
    expect_true("(Main)" %in% cache$valid.altExp.names)
})

test_that(".refineParameters resolves NA Experiment and Type", {
    x   <- AltReducedDimensionPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[["Experiment"]], "peptides")
    expect_identical(y[["Type"]], "PCA")
    expect_identical(y[["XAxis"]], 1L)
    expect_identical(y[["YAxis"]], 2L)
})

test_that(".refineParameters works for (Main) sentinel", {
    x   <- AltReducedDimensionPlot(Experiment = "(Main)")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[["Experiment"]], "(Main)")
    expect_identical(y[["Type"]], "PCA")
})

test_that(".refineParameters sets Type to NA when no reducedDims available", {
    sce0 <- sce
    # Remove reducedDims from peptides altExp
    reducedDims(altExp(sce0, "peptides")) <- SimpleList()
    x   <- AltReducedDimensionPlot(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce0)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_true(is.na(y[["Type"]]))
})

test_that(".defineDataInterface produces a list without error", {
    x   <- AltReducedDimensionPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    expect_error(iSEE:::.defineDataInterface(x, se2, list()), NA)
})

test_that(".generateDotPlotData produces X/Y from reducedDim", {
    x   <- AltReducedDimensionPlot(Experiment = "peptides", Type = "PCA")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)

    expect_true(is.data.frame(envir$plot.data))
    expect_identical(nrow(envir$plot.data), ncol(sce))
    rd <- reducedDim(altExp(se2, "peptides"), "PCA")
    expect_equal(unname(envir$plot.data$X), unname(rd[, 1]))
    expect_equal(unname(envir$plot.data$Y), unname(rd[, 2]))
    expect_identical(names(out), c("commands", "labels"))
})

test_that(".generateDotPlotData returns NA-filled data when Type is NA", {
    x   <- AltReducedDimensionPlot(Experiment = "peptides")
    sce0 <- sce
    reducedDims(altExp(sce0, "peptides")) <- SimpleList()
    se2 <- iSEE:::.cacheCommonInfo(x, sce0)
    x   <- iSEE:::.refineParameters(x, se2)

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)
    expect_true(all(is.na(envir$plot.data$X)))
    expect_true(all(is.na(envir$plot.data$Y)))
})

test_that(".generateDotPlot renders a ggplot for valid data", {
    library(ggplot2)
    x   <- AltReducedDimensionPlot(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)

    envir           <- new.env()
    envir$se        <- se2
    envir$colormap  <- iSEE::ExperimentColorMap()
    envir$plot.type <- "scatter"
    out <- iSEE:::.generateDotPlotData(x, envir)
    iSEE:::.generateDotPlot(x, out$labels, envir)
    expect_s3_class(envir$dot.plot, "ggplot")
})

test_that(".generateDotPlot renders error placeholder when Type is NA", {
    library(ggplot2)
    x   <- AltReducedDimensionPlot(Experiment = "peptides")
    sce0 <- sce
    reducedDims(altExp(sce0, "peptides")) <- SimpleList()
    se2 <- iSEE:::.cacheCommonInfo(x, sce0)
    x   <- iSEE:::.refineParameters(x, se2)

    envir          <- new.env()
    envir$se       <- se2
    envir$colormap <- iSEE::ExperimentColorMap()
    out <- iSEE:::.generateDotPlotData(x, envir)
    result <- iSEE:::.generateDotPlot(x, out$labels, envir)
    expect_s3_class(result$plot, "ggplot")
})

test_that(".definePanelTour returns a valid data.frame", {
    x   <- AltReducedDimensionPlot()
    tour <- iSEE:::.definePanelTour(x)
    expect_s3_class(tour, "data.frame")
    expect_identical(colnames(tour), c("element", "intro"))
    expect_true(any(grepl("alternative experiment", tour$intro)))
})
