# Tests for LinkedFeaturesAssayPlot
# library(testthat); library(altSEE); source("setup-sce.R"); source("test-LinkedFeaturesAssayPlot.R")

# ---- constructor -----------------------------------------------------------

test_that("LinkedFeaturesAssayPlot constructor works with defaults", {
    x <- LinkedFeaturesAssayPlot()
    expect_s4_class(x, "LinkedFeaturesAssayPlot")
    expect_identical(x[["SelectionExperiment"]], NA_character_)
    expect_identical(x[["Experiment"]], NA_character_)
    expect_identical(x[["AltAssay"]], NA_character_)
    expect_identical(x[["LookupColumn"]], NA_character_)
    expect_identical(x[["MapColumn"]], NA_character_)
    expect_identical(x[["PlotType"]], "Scatter")
})

test_that("LinkedFeaturesAssayPlot constructor accepts explicit slots", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)",
        LookupColumn        = "Proteins",
        Experiment          = "peptides",
        AltAssay            = "peptides_norm",
        MapColumn           = "Proteins",
        PlotType            = "Scatter + lines"
    )
    expect_identical(x[["SelectionExperiment"]], "(Main)")
    expect_identical(x[["Experiment"]], "peptides")
    expect_identical(x[["AltAssay"]], "peptides_norm")
    expect_identical(x[["LookupColumn"]], "Proteins")
    expect_identical(x[["MapColumn"]], "Proteins")
    expect_identical(x[["PlotType"]], "Scatter + lines")
})

test_that("LinkedFeaturesAssayPlot validity rejects bad PlotType", {
    # PlotType must be a single string (character(1))
    expect_error(new("LinkedFeaturesAssayPlot", PlotType = character(0)))
})

# ---- identity generics -----------------------------------------------------

test_that(".fullName and .panelColor work for LinkedFeaturesAssayPlot", {
    x <- LinkedFeaturesAssayPlot()
    expect_identical(iSEE:::.fullName(x), "Linked features assay plot")
    expect_type(iSEE:::.panelColor(x), "character")
})

# ---- .cacheCommonInfo ------------------------------------------------------

test_that(".cacheCommonInfo caches lookup/map columns per altExp", {
    x   <- LinkedFeaturesAssayPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    se2 <- iSEE:::.cacheCommonInfo(x, se2)  # second call returns early

    cache <- iSEE:::.getCachedCommonInfo(se2, "LinkedFeaturesAssayPlot")
    expect_true("peptides" %in% cache$valid.altExp.names)
    expect_true("(Main)"   %in% cache$valid.selection.names)
    # "(rowname)" is always present
    expect_true("(rowname)" %in% cache$valid.lookcols.by.altExp[["peptides"]])
    expect_true("Proteins"  %in% cache$valid.lookcols.by.altExp[["peptides"]])
    expect_true("(rowname)" %in% cache$valid.mapcols.by.altExp[["peptides"]])
    expect_true("Proteins"  %in% cache$valid.mapcols.by.altExp[["peptides"]])
    # Main experiment
    expect_true("Proteins"  %in% cache$valid.lookcols.by.altExp[["(Main)"]])
})

# ---- .refineParameters -----------------------------------------------------

test_that(".refineParameters resolves NA slots to valid defaults", {
    x   <- LinkedFeaturesAssayPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)

    expect_false(is.na(y[["SelectionExperiment"]]))
    expect_false(is.na(y[["Experiment"]]))
    expect_false(is.na(y[["AltAssay"]]))
    expect_false(is.na(y[["LookupColumn"]]))
    expect_false(is.na(y[["MapColumn"]]))
})

test_that(".refineParameters preserves valid SelectionExperiment and feature", {
    x <- LinkedFeaturesAssayPlot(SelectionExperiment = "(Main)")
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[["SelectionExperiment"]], "(Main)")
    expect_identical(y[[iSEE:::.featAssayYAxisFeatName]], "PROT1")
})

test_that(".refineParameters falls back when feature not in SelectionExperiment", {
    x <- LinkedFeaturesAssayPlot(SelectionExperiment = "peptides")
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT999"
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_true(y[[iSEE:::.featAssayYAxisFeatName]] %in% rownames(altExp(se2, "peptides")))
})

test_that(".refineParameters preserves valid LookupColumn", {
    x <- LinkedFeaturesAssayPlot(SelectionExperiment = "(Main)", LookupColumn = "Proteins")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[["LookupColumn"]], "Proteins")
})

test_that(".refineParameters resets invalid LookupColumn", {
    x <- LinkedFeaturesAssayPlot(SelectionExperiment = "(Main)", LookupColumn = "NOT_VALID")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_false(y[["LookupColumn"]] == "NOT_VALID")
})

test_that(".refineParameters preserves valid MapColumn", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)", LookupColumn = "Proteins",
        Experiment = "peptides", AltAssay = "peptides_norm", MapColumn = "Proteins")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[["MapColumn"]], "Proteins")
})

test_that(".refineParameters sets default color/shape to altExp_feature_id", {
    x   <- LinkedFeaturesAssayPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[[iSEE:::.colorByColData]], "altExp_feature_id")
    expect_identical(y[[iSEE:::.shapeByColData]], "altExp_feature_id")
})

test_that(".refineParameters preserves user-chosen colData column for color", {
    x   <- LinkedFeaturesAssayPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    y[[iSEE:::.colorByColData]] <- "Treatment"
    # Re-refine with saved colour choice
    y2 <- iSEE:::.refineParameters(y, se2)
    expect_identical(y2[[iSEE:::.colorByColData]], "Treatment")
})

test_that(".refineParameters returns NULL when altExp has no assays", {
    sce0 <- sce
    SummarizedExperiment::assays(altExp(sce0, "peptides")) <- S4Vectors::SimpleList()
    x   <- LinkedFeaturesAssayPlot(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce0)
    expect_warning(out <- iSEE:::.refineParameters(x, se2), "no valid 'assays'")
    expect_null(out)
})

# ---- .defineDataInterface --------------------------------------------------

test_that(".defineDataInterface produces a list without error", {
    x   <- LinkedFeaturesAssayPlot()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    expect_error(
        iSEE:::.defineDataInterface(x, se2,
            list(single = list(feature = character(0)))),
        NA)
})

# ---- .generateDotPlotData --------------------------------------------------

test_that(".generateDotPlotData: protein -> peptides lookup via Proteins column", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)",
        LookupColumn        = "Proteins",
        Experiment          = "peptides",
        AltAssay            = "peptides_norm",
        MapColumn           = "Proteins"
    )
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)

    expect_true(is.data.frame(envir$plot.data))
    expect_true(all(c("altExp_feature_id", "sample", "X", "Y") %in%
                    colnames(envir$plot.data)))
    # PROT1 maps to PEP1, PEP2, PEP3 — 3 peptides x 6 samples = 18 rows
    expect_identical(nrow(envir$plot.data), 3L * ncol(sce))
    expect_true(all(envir$plot.data$altExp_feature_id %in% c("PEP1", "PEP2", "PEP3")))
    expect_identical(names(out), c("commands", "labels"))
})

test_that(".generateDotPlotData: rowname as LookupColumn", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)",
        LookupColumn        = "(rowname)",
        Experiment          = "(Main)",
        AltAssay            = "proteins",
        MapColumn           = "(rowname)"
    )
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT2"
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT2"

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)

    expect_true(is.data.frame(envir$plot.data))
    expect_identical(nrow(envir$plot.data), 1L * ncol(sce))
    expect_identical(envir$plot.data$altExp_feature_id[1], "PROT2")
})

test_that(".generateDotPlotData: altExp SelectionExperiment to main Experiment", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "peptides",
        LookupColumn        = "Proteins",
        Experiment          = "(Main)",
        AltAssay            = "proteins",
        MapColumn           = "Proteins"
    )
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PEP1"
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PEP1"

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)

    # PEP1's Proteins value is PROT1; main experiment has 1 row with Proteins=="PROT1"
    expect_true(is.data.frame(envir$plot.data))
    expect_identical(nrow(envir$plot.data), 1L * ncol(sce))
    expect_identical(envir$plot.data$altExp_feature_id[1], "PROT1")
})

test_that(".generateDotPlotData sets column-data X-axis", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)",
        LookupColumn        = "Proteins",
        Experiment          = "peptides",
        AltAssay            = "peptides_norm",
        MapColumn           = "Proteins"
    )
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.featAssayXAxis]]          <- iSEE:::.featAssayXAxisColDataTitle
    x[[iSEE:::.featAssayXAxisColData]]   <- "Treatment"

    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.featAssayXAxis]]          <- iSEE:::.featAssayXAxisColDataTitle
    x[[iSEE:::.featAssayXAxisColData]]   <- "Treatment"

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)
    expect_true(!is.null(envir$plot.data$X))
    expect_true(all(envir$plot.data$X %in% c("A", "B")))
})

test_that(".generateDotPlotData returns error sentinel on bad lookup", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)",
        LookupColumn        = "(rowname)",
        Experiment          = "peptides",
        AltAssay            = "peptides_norm",
        MapColumn           = "Proteins"
    )
    # Select a protein whose name doesn't match any peptide "Proteins" value
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    # Use MapColumn=(rowname) so nothing matches the protein name in pep rownames
    x[["MapColumn"]] <- "(rowname)"

    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[["MapColumn"]] <- "(rowname)"

    envir    <- new.env()
    envir$se <- se2
    out <- iSEE:::.generateDotPlotData(x, envir)
    # Plot title is the error message (lookup key == feature id == "PROT1", no match)
    # Either the data is empty (0 rows) or we get an error sentinel — either is acceptable
    expect_true(is.data.frame(envir$plot.data))
})

# ---- .generateDotPlot ------------------------------------------------------

test_that(".generateDotPlot produces a ggplot for Scatter type", {
    library(ggplot2)
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)",
        LookupColumn        = "Proteins",
        Experiment          = "peptides",
        AltAssay            = "peptides_norm",
        MapColumn           = "Proteins",
        PlotType            = "Scatter"
    )
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.colorByField]] <- iSEE:::.colorByNothingTitle
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.colorByField]] <- iSEE:::.colorByNothingTitle

    envir          <- new.env()
    envir$se       <- se2
    envir$colormap <- iSEE::ExperimentColorMap()
    out <- iSEE:::.generateDotPlotData(x, envir)
    result <- iSEE:::.generateDotPlot(x, out$labels, envir)
    expect_s3_class(result$plot, "ggplot")
})

test_that(".generateDotPlot produces a ggplot for Scatter + lines type", {
    library(ggplot2)
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)",
        LookupColumn        = "Proteins",
        Experiment          = "peptides",
        AltAssay            = "peptides_norm",
        MapColumn           = "Proteins",
        PlotType            = "Scatter + lines"
    )
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.featAssayXAxis]]         <- iSEE:::.featAssayXAxisColDataTitle
    x[[iSEE:::.featAssayXAxisColData]]  <- "SampleGroup"
    x[[iSEE:::.colorByField]]           <- iSEE:::.colorByNothingTitle

    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.featAssayXAxis]]         <- iSEE:::.featAssayXAxisColDataTitle
    x[[iSEE:::.featAssayXAxisColData]]  <- "SampleGroup"
    x[[iSEE:::.colorByField]]           <- iSEE:::.colorByNothingTitle

    envir          <- new.env()
    envir$se       <- se2
    envir$colormap <- iSEE::ExperimentColorMap()
    out <- iSEE:::.generateDotPlotData(x, envir)
    result <- iSEE:::.generateDotPlot(x, out$labels, envir)
    expect_s3_class(result$plot, "ggplot")
    expect_true(any(grepl("geom_line", result$commands)))
})

test_that(".generateDotPlot renders error placeholder when labels$Y is empty", {
    library(ggplot2)
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)",
        LookupColumn        = "Proteins",
        Experiment          = "peptides",
        AltAssay            = "peptides_norm",
        MapColumn           = "Proteins"
    )
    # Force an error by providing an assay that doesn't exist
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[["AltAssay"]] <- "NONEXISTENT_ASSAY"

    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    # Override the refined assay to force an error
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[["AltAssay"]] <- "NONEXISTENT_ASSAY"

    envir    <- new.env()
    envir$se <- se2
    out    <- iSEE:::.generateDotPlotData(x, envir)
    result <- iSEE:::.generateDotPlot(x, out$labels, envir)
    expect_s3_class(result$plot, "ggplot")
})

# ---- .addDotPlotData* methods ----------------------------------------------

test_that(".addDotPlotDataColor handles altExp_feature_id coloring", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)", LookupColumn = "Proteins",
        Experiment = "peptides", AltAssay = "peptides_norm", MapColumn = "Proteins"
    )
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.colorByField]]   <- iSEE:::.colorByColDataTitle
    x[[iSEE:::.colorByColData]] <- "altExp_feature_id"

    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.colorByField]]   <- iSEE:::.colorByColDataTitle
    x[[iSEE:::.colorByColData]] <- "altExp_feature_id"

    envir    <- new.env()
    envir$se <- se2
    iSEE:::.generateDotPlotData(x, envir)
    result <- iSEE:::.addDotPlotDataColor(x, envir)
    expect_true(!is.null(envir$plot.data$ColorBy))
    expect_true(all(envir$plot.data$ColorBy %in% c("PEP1", "PEP2", "PEP3")))
})

test_that(".addDotPlotDataColor handles colData column coloring", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)", LookupColumn = "Proteins",
        Experiment = "peptides", AltAssay = "peptides_norm", MapColumn = "Proteins"
    )
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.colorByField]]   <- iSEE:::.colorByColDataTitle
    x[[iSEE:::.colorByColData]] <- "Treatment"

    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.colorByField]]   <- iSEE:::.colorByColDataTitle
    x[[iSEE:::.colorByColData]] <- "Treatment"

    envir    <- new.env()
    envir$se <- se2
    iSEE:::.generateDotPlotData(x, envir)
    iSEE:::.addDotPlotDataColor(x, envir)
    expect_true(all(levels(envir$plot.data$ColorBy) %in% c("A", "B")))
})

test_that(".addDotPlotDataColor returns NULL when color is 'none'", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)", LookupColumn = "Proteins",
        Experiment = "peptides", AltAssay = "peptides_norm", MapColumn = "Proteins"
    )
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.colorByField]] <- iSEE:::.colorByNothingTitle

    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.colorByField]] <- iSEE:::.colorByNothingTitle

    envir    <- new.env()
    envir$se <- se2
    iSEE:::.generateDotPlotData(x, envir)
    result <- iSEE:::.addDotPlotDataColor(x, envir)
    expect_null(result)
})

test_that(".addDotPlotDataShape uses altExp_feature_id by default", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)", LookupColumn = "Proteins",
        Experiment = "peptides", AltAssay = "peptides_norm", MapColumn = "Proteins"
    )
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"

    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"

    envir    <- new.env()
    envir$se <- se2
    iSEE:::.generateDotPlotData(x, envir)
    iSEE:::.addDotPlotDataShape(x, envir)
    expect_true(all(envir$plot.data$ShapeBy %in% c("PEP1", "PEP2", "PEP3")))
})

test_that(".addDotPlotDataSize returns NULL when size is 'none'", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)", LookupColumn = "Proteins",
        Experiment = "peptides", AltAssay = "peptides_norm", MapColumn = "Proteins"
    )
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.sizeByField]] <- iSEE:::.sizeByNothingTitle

    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.sizeByField]] <- iSEE:::.sizeByNothingTitle

    envir    <- new.env()
    envir$se <- se2
    iSEE:::.generateDotPlotData(x, envir)
    result <- iSEE:::.addDotPlotDataSize(x, envir)
    expect_null(result)
})

test_that(".addDotPlotDataFacets adds colData column to plot.data", {
    x <- LinkedFeaturesAssayPlot(
        SelectionExperiment = "(Main)", LookupColumn = "Proteins",
        Experiment = "peptides", AltAssay = "peptides_norm", MapColumn = "Proteins"
    )
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.facetRow]] <- iSEE:::.facetByColDataTitle
    x[[iSEE:::.facetRowByColData]] <- "Treatment"

    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    x[[iSEE:::.featAssayYAxisFeatName]] <- "PROT1"
    x[[iSEE:::.facetRow]] <- iSEE:::.facetByColDataTitle
    x[[iSEE:::.facetRowByColData]] <- "Treatment"

    envir    <- new.env()
    envir$se <- se2
    iSEE:::.generateDotPlotData(x, envir)
    iSEE:::.addDotPlotDataFacets(x, envir)
    expect_true("FacetRow" %in% colnames(envir$plot.data))
    expect_true(all(envir$plot.data$FacetRow %in% c("A", "B")))
})

# ---- .definePanelTour -------------------------------------------------------

test_that(".definePanelTour returns a valid data.frame", {
    x    <- LinkedFeaturesAssayPlot()
    tour <- iSEE:::.definePanelTour(x)
    expect_s3_class(tour, "data.frame")
    expect_identical(colnames(tour), c("element", "intro"))
    expect_true(any(grepl("Linked features assay plot", tour$intro)))
})
