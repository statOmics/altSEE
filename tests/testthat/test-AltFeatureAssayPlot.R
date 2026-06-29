# Tests for AltFeatureAssayPlot
# library(testthat); library(altSEE); source("setup-sce.R"); source("test-AltFeatureAssayPlot.R")

test_that("AltFeatureAssayPlot constructor works", {
    x <- AltFeatureAssayPlot()
    expect_s4_class(x, "AltFeatureAssayPlot")
    expect_identical(x[["Experiment"]], NA_character_)
})

test_that("AltFeatureAssayPlot constructor accepts slot values", {
    x <- AltFeatureAssayPlot(Experiment = "peptides")
    expect_identical(x[["Experiment"]], "peptides")
})

test_that(".fullName and .panelColor work for AltFeatureAssayPlot", {
    x <- AltFeatureAssayPlot()
    expect_identical(iSEE:::.fullName(x), "Alt feature assay plot")
    expect_type(iSEE:::.panelColor(x), "character")
})

test_that(".cacheCommonInfo caches assay names per altExp", {
    x   <- AltFeatureAssayPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    se2 <- iSEE:::.cacheCommonInfo(x, se2)  # second call returns early

    cache <- iSEE:::.getCachedCommonInfo(se2, "AltFeatureAssayPlot")
    expect_true("peptides_norm" %in% cache$valid.assay.names.by.altExp[["peptides"]])
    expect_true("proteins" %in% cache$valid.assay.names.by.altExp[["(Main)"]])
})

test_that(".refineParameters resolves NA Experiment and sets first feature", {
    x   <- AltFeatureAssayPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[["Experiment"]], "peptides")
    expect_true(y[[iSEE:::.featAssayYAxisFeatName]] %in% rownames(altExp(se2, "peptides")))
})

test_that(".refineParameters resolves to main experiment via sentinel", {
    x   <- AltFeatureAssayPlot(Experiment = "(Main)")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[["Experiment"]], "(Main)")
    expect_true(y[[iSEE:::.featAssayYAxisFeatName]] %in% rownames(se2))
})

test_that(".refineParameters preserves valid Y-axis feature from altExp", {
    x   <- AltFeatureAssayPlot(Experiment = "peptides")
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PEP3"
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[[iSEE:::.featAssayYAxisFeatName]], "PEP3")
})

test_that(".refineParameters returns NULL if altExp has no rows", {
    sce0 <- sce
    altExp(sce0, "peptides") <- altExp(sce0, "peptides")[integer(0), ]
    x   <- AltFeatureAssayPlot(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce0)
    expect_warning(out <- iSEE:::.refineParameters(x, se2), "no rows")
    expect_null(out)
})

test_that(".defineDataInterface produces a list without error", {
    x   <- AltFeatureAssayPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    expect_error(iSEE:::.defineDataInterface(x, se2, list(single = list(feature = character(0)))), NA)
})

test_that(".generateDotPlotData builds plot.data for peptides altExp", {
    x   <- AltFeatureAssayPlot(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)

    expect_true(is.data.frame(envir$plot.data))
    expect_identical(nrow(envir$plot.data), ncol(sce))
    expect_true("Y" %in% colnames(envir$plot.data))
    expect_identical(names(out), c("commands", "labels"))
})

test_that(".generateDotPlotData works for (Main) experiment", {
    x   <- AltFeatureAssayPlot(Experiment = "(Main)")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)
    expect_identical(nrow(envir$plot.data), ncol(sce))
})

test_that(".generateDotPlotData handles column-data x-axis", {
    x   <- AltFeatureAssayPlot(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayXAxis]]       <- iSEE:::.featAssayXAxisColDataTitle
    x[[iSEE:::.featAssayXAxisColData]] <- "Treatment"

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)
    expect_true(!is.null(envir$plot.data$X))
    expect_identical(levels(envir$plot.data$X), c("A", "B"))
})

test_that(".generateDotPlotData returns error sentinel on bad feature", {
    x   <- AltFeatureAssayPlot(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "NOT_A_FEATURE"

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)
    expect_identical(out$labels$Y, "")
})

test_that(".generateDotPlot produces a ggplot object", {
    library(ggplot2)
    x   <- AltFeatureAssayPlot(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)

    envir           <- new.env()
    envir$se        <- se2
    envir$colormap  <- iSEE::ExperimentColorMap()
    envir$plot.type <- "violin"
    out <- iSEE:::.generateDotPlotData(x, envir)
    iSEE:::.generateDotPlot(x, out$labels, envir)
    expect_s3_class(envir$dot.plot, "ggplot")
})

test_that(".generateDotPlot renders error plot when labels$Y is empty", {
    library(ggplot2)
    x   <- AltFeatureAssayPlot(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "NOT_A_FEATURE"

    envir    <- new.env()
    envir$se <- se2
    out      <- iSEE:::.generateDotPlotData(x, envir)
    result   <- iSEE:::.generateDotPlot(x, out$labels, envir)
    expect_s3_class(result$plot, "ggplot")
})

test_that(".definePanelTour returns a valid data.frame", {
    x   <- AltFeatureAssayPlot()
    tour <- iSEE:::.definePanelTour(x)
    expect_s3_class(tour, "data.frame")
    expect_identical(colnames(tour), c("element", "intro"))
    expect_true(any(grepl("alternative experiment", tour$intro)))
})
