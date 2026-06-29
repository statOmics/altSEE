# Tests for AltVolcanoPlot
# library(testthat); library(altSEE); source("setup-sce.R"); source("test-AltVolcanoPlot.R")

test_that("AltVolcanoPlot constructor works", {
    x <- AltVolcanoPlot()
    expect_s4_class(x, "AltVolcanoPlot")
    expect_identical(x[["Experiment"]], NA_character_)
})

test_that("AltVolcanoPlot constructor accepts slot values", {
    x <- AltVolcanoPlot(Experiment = "peptides")
    expect_identical(x[["Experiment"]], "peptides")
})

test_that(".fullName and .panelColor work for AltVolcanoPlot", {
    x <- AltVolcanoPlot()
    expect_identical(iSEE:::.fullName(x), "Alt volcano plot")
    expect_type(iSEE:::.panelColor(x), "character")
})

test_that(".cacheCommonInfo caches p-value and logFC fields per altExp", {
    x   <- AltVolcanoPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    se2 <- iSEE:::.cacheCommonInfo(x, se2)  # second call returns early

    cache <- iSEE:::.getCachedCommonInfo(se2, "AltVolcanoPlot")
    expect_true("peptides" %in% names(cache$valid.p.fields.by.altExp))
    expect_true(length(cache$valid.p.fields.by.altExp[["peptides"]]) > 0L)
    expect_true(length(cache$valid.lfc.fields.by.altExp[["peptides"]]) > 0L)
    expect_true(length(cache$valid.p.fields.by.altExp[["(Main)"]]) > 0L)
})

test_that(".refineParameters resolves NA Experiment and sets p-value/logFC fields", {
    x   <- AltVolcanoPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[["Experiment"]], "peptides")
    expect_false(is.na(y[[iSEE:::.rowDataYAxis]]))
    expect_false(is.na(y[[iSEE:::.rowDataXAxisRowData]]))
    expect_identical(y[["XAxis"]], "Row data")
})

test_that(".refineParameters sets YAxis to NA when no p-value fields exist", {
    sce0 <- sce
    rowData(altExp(sce0, "peptides"))$PValue <- NULL
    x   <- AltVolcanoPlot(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce0)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_true(is.na(y[[iSEE:::.rowDataYAxis]]))
})

test_that(".refineParameters works for (Main) experiment", {
    x   <- AltVolcanoPlot(Experiment = "(Main)")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[["Experiment"]], "(Main)")
    expect_false(is.na(y[[iSEE:::.rowDataYAxis]]))
})

test_that(".defineDataInterface produces a list without error", {
    x   <- AltVolcanoPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    expect_error(iSEE:::.defineDataInterface(x, se2, list()), NA)
})

test_that(".generateDotPlotData produces X/Y from altExp rowData", {
    x   <- AltVolcanoPlot(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)

    expect_true(is.data.frame(envir$plot.data))
    expect_identical(nrow(envir$plot.data), nrow(altExp(se2, "peptides")))
    expect_true("IsSig" %in% colnames(envir$plot.data))
    # Y should be -log10(p-value)
    pval_col <- x[[iSEE:::.rowDataYAxis]]
    raw_pvals <- rowData(altExp(se2, "peptides"))[[pval_col]]
    expect_equal(envir$plot.data$Y, -log10(raw_pvals))
    expect_identical(names(out), c("commands", "labels"))
})

test_that(".generateDotPlotData produces valid data for (Main) experiment", {
    x   <- AltVolcanoPlot(Experiment = "(Main)")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)
    expect_identical(nrow(envir$plot.data), nrow(se2))
})

test_that(".generateDotPlotData returns error sentinel when YAxis is NA", {
    x   <- AltVolcanoPlot(Experiment = "peptides")
    sce0 <- sce
    rowData(altExp(sce0, "peptides"))$PValue <- NULL
    se2 <- iSEE:::.cacheCommonInfo(x, sce0)
    x   <- iSEE:::.refineParameters(x, se2)

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)
    expect_identical(out$labels$X, "")
})

test_that(".generateDotPlot produces a ggplot for valid data", {
    library(ggplot2)
    x   <- AltVolcanoPlot(Experiment = "peptides")
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

test_that(".generateDotPlot renders error placeholder when YAxis is NA", {
    library(ggplot2)
    x   <- AltVolcanoPlot(Experiment = "peptides")
    sce0 <- sce
    rowData(altExp(sce0, "peptides"))$PValue <- NULL
    se2 <- iSEE:::.cacheCommonInfo(x, sce0)
    x   <- iSEE:::.refineParameters(x, se2)

    envir          <- new.env()
    envir$se       <- se2
    envir$colormap <- iSEE::ExperimentColorMap()
    out    <- iSEE:::.generateDotPlotData(x, envir)
    result <- iSEE:::.generateDotPlot(x, out$labels, envir)
    expect_s3_class(result$plot, "ggplot")
})

test_that(".definePanelTour returns a valid data.frame", {
    x    <- AltVolcanoPlot()
    tour <- iSEE:::.definePanelTour(x)
    expect_s3_class(tour, "data.frame")
    expect_identical(colnames(tour), c("element", "intro"))
    expect_true(any(grepl("alternative experiment", tour$intro)))
})
