# Tests for .scatter_plot2 internal helper
# library(testthat); library(altSEE); source("setup-sce.R"); source("test-scatter_plot2.R")

.eval_cmds <- function(cmds, envir) {
    eval(parse(text = paste(cmds, collapse = "\n")), envir = envir)
}

.make_scatter2_inputs <- function(x_numeric = TRUE) {
    n_pep  <- 3L
    n_samp <- 6L
    pep_ids  <- paste0("PEP", seq_len(n_pep))
    samp_ids <- paste0("S",   seq_len(n_samp))
    plot_data <- data.frame(
        altExp_feature_id = rep(pep_ids, each = n_samp),
        sample  = rep(samp_ids, times = n_pep),
        X = if (x_numeric) rnorm(n_pep * n_samp) else
            rep(factor(c("ctrl","trt"), levels = c("ctrl","trt")), n_pep * n_samp / 2),
        Y = rnorm(n_pep * n_samp),
        stringsAsFactors = FALSE
    )
    param <- LinkedFeaturesAssayPlot()
    se2   <- iSEE:::.cacheCommonInfo(param, sce)
    param <- iSEE:::.refineParameters(param, se2)
    list(plot_data = plot_data, param = param)
}

test_that(".scatter_plot2 returns a named character vector", {
    inp   <- .make_scatter2_inputs()
    result <- altSEE:::.scatter_plot2(
        inp$plot_data, inp$param,
        x_lab = "X", y_lab = "Y",
        color_lab = NULL, shape_lab = NULL, size_lab = NULL,
        title = "Test")
    expect_type(result, "character")
    expect_true(length(result) > 0L)
    expect_true(!is.null(names(result)))
})

test_that(".scatter_plot2 builds ggplot when evaluated (no color)", {
    library(ggplot2)
    inp   <- .make_scatter2_inputs()
    pd    <- inp$plot_data
    param <- inp$param
    param[[iSEE:::.colorByField]] <- iSEE:::.colorByNothingTitle

    cmds <- altSEE:::.scatter_plot2(
        pd, param,
        x_lab = "X", y_lab = "Y",
        color_lab = NULL, shape_lab = NULL, size_lab = NULL,
        title = "Test")

    envir           <- new.env()
    envir$plot.data <- pd
    envir$colormap  <- iSEE::ExperimentColorMap()
    .eval_cmds(cmds, envir)
    expect_s3_class(envir$dot.plot, "ggplot")
})

test_that(".scatter_plot2 works with character X axis (no color)", {
    library(ggplot2)
    inp   <- .make_scatter2_inputs(x_numeric = FALSE)
    pd    <- inp$plot_data
    param <- inp$param
    param[[iSEE:::.colorByField]] <- iSEE:::.colorByNothingTitle

    cmds <- altSEE:::.scatter_plot2(
        pd, param,
        x_lab = "group", y_lab = "Y",
        color_lab = NULL, shape_lab = NULL, size_lab = NULL,
        title = "Test")

    envir           <- new.env()
    envir$plot.data <- pd
    envir$colormap  <- iSEE::ExperimentColorMap()
    .eval_cmds(cmds, envir)
    expect_s3_class(envir$dot.plot, "ggplot")
})

test_that(".scatter_plot2 returns guard plot for > 1000 distinct features", {
    pd <- data.frame(
        altExp_feature_id = paste0("F", seq_len(1001)),
        id     = paste0("F", seq_len(1001)),
        sample = "S1",
        X      = rnorm(1001),
        Y      = rnorm(1001),
        stringsAsFactors = FALSE
    )
    param <- LinkedFeaturesAssayPlot()
    se2   <- iSEE:::.cacheCommonInfo(param, sce)
    param <- iSEE:::.refineParameters(param, se2)

    result <- altSEE:::.scatter_plot2(
        pd, param,
        x_lab = "X", y_lab = "Y",
        color_lab = NULL, shape_lab = NULL, size_lab = NULL,
        title = "Test")
    expect_true(any(grepl("More than 1000", result)))
})

test_that(".scatter_plot2 handles numeric ColorBy in plot.data", {
    library(ggplot2)
    inp  <- .make_scatter2_inputs()
    pd   <- inp$plot_data
    # Use a numeric ColorBy (e.g., a simulated expression value per row)
    pd$ColorBy <- rnorm(nrow(pd))

    param <- inp$param
    param[[iSEE:::.colorByField]]   <- iSEE:::.colorByColDataTitle
    param[[iSEE:::.colorByColData]] <- "AveExpr"

    cmds <- altSEE:::.scatter_plot2(
        pd, param,
        x_lab = "X", y_lab = "Y",
        color_lab = "AveExpr", shape_lab = NULL, size_lab = NULL,
        title = "Test")

    envir           <- new.env()
    envir$plot.data <- pd
    envir$colormap  <- iSEE::ExperimentColorMap()
    .eval_cmds(cmds, envir)
    expect_s3_class(envir$dot.plot, "ggplot")
})
