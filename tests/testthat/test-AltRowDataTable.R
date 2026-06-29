# Tests for AltRowDataTable
# library(testthat); library(altSEE); source("setup-sce.R"); source("test-AltRowDataTable.R")

test_that("AltRowDataTable constructor works", {
    x <- AltRowDataTable()
    expect_s4_class(x, "AltRowDataTable")
    expect_identical(x[["Experiment"]], NA_character_)
})

test_that("AltRowDataTable constructor accepts slot values", {
    x <- AltRowDataTable(Experiment = "peptides")
    expect_identical(x[["Experiment"]], "peptides")
})

test_that(".fullName and .panelColor work for AltRowDataTable", {
    x <- AltRowDataTable()
    expect_identical(iSEE:::.fullName(x), "Alt row data table")
    expect_type(iSEE:::.panelColor(x), "character")
})

test_that(".cacheCommonInfo caches altExp names and rowData columns", {
    x  <- AltRowDataTable()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    se2 <- iSEE:::.cacheCommonInfo(x, se2)  # second call returns early

    cache <- iSEE:::.getCachedCommonInfo(se2, "AltRowDataTable")
    expect_identical(cache$valid.altExp.names, c("peptides", "(Main)"))
    expect_true("Proteins" %in% cache$valid.rowData.names.by.altExp[["peptides"]])
    expect_true("Proteins" %in% cache$valid.rowData.names.by.altExp[["(Main)"]])
})

test_that(".refineParameters resolves NA Experiment to first altExp", {
    x   <- AltRowDataTable()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[["Experiment"]], "peptides")
})

test_that(".refineParameters preserves valid Experiment slot", {
    x   <- AltRowDataTable(Experiment = "(Main)")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[["Experiment"]], "(Main)")
})

test_that(".refineParameters preserves valid Selected feature", {
    x   <- AltRowDataTable(Experiment = "peptides", Selected = "PEP1")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_identical(y[["Selected"]], "PEP1")
})

test_that(".refineParameters drops invalid Selected feature", {
    x   <- AltRowDataTable(Experiment = "peptides", Selected = "NOT_A_FEATURE")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    y   <- iSEE:::.refineParameters(x, se2)
    expect_false(y[["Selected"]] == "NOT_A_FEATURE")
})

test_that(".defineDataInterface produces a list without error", {
    x   <- AltRowDataTable()
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)
    expect_error(iSEE:::.defineDataInterface(x, se2, list()), NA)
})

test_that(".generateTable returns peptide rowData for peptides altExp", {
    x   <- AltRowDataTable(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)

    envir    <- new.env()
    envir$se <- se2
    iSEE:::.generateTable(x, envir)
    expect_true(is.data.frame(envir$tab))
    expect_identical(rownames(envir$tab), rownames(altExp(se2, "peptides")))
    expect_true("Proteins" %in% colnames(envir$tab))
})

test_that(".generateTable returns main experiment rowData for (Main)", {
    x   <- AltRowDataTable(Experiment = "(Main)")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)

    envir    <- new.env()
    envir$se <- se2
    iSEE:::.generateTable(x, envir)
    expect_identical(rownames(envir$tab), rownames(se2))
})

test_that(".generateTable respects row_selected from another panel", {
    x   <- AltRowDataTable(Experiment = "peptides")
    se2 <- iSEE:::.cacheCommonInfo(x, sce)
    x   <- iSEE:::.refineParameters(x, se2)

    envir              <- new.env()
    envir$se           <- se2
    envir$row_selected <- list(active = c("PEP1", "PEP2"))
    iSEE:::.generateTable(x, envir)
    expect_true(all(rownames(envir$tab) %in% c("PEP1", "PEP2")))
})

test_that(".definePanelTour returns a valid data.frame", {
    x   <- AltRowDataTable()
    tour <- iSEE:::.definePanelTour(x)
    expect_s3_class(tour, "data.frame")
    expect_identical(colnames(tour), c("element", "intro"))
    expect_true(any(grepl("alternative experiment", tour$intro)))
})
