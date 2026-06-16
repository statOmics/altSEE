############################################################
# AltExpFeatureAssayPlot
############################################################

# Extends FeatureAssayPlot to plot assay values from an alternative experiment.
# The Y-axis feature is typically driven by an AltExpRowDataTable via
# YAxisFeatureSource, so the feature name received is an altExp rowname.
# The assay is taken from the selected altExp rather than the main experiment.

#' AltExpFeatureAssayPlot: feature assay plot for alternative experiments
#'
#' An S4 class extending \code{\linkS4class{FeatureAssayPlot}} to plot
#' assay values from a feature in an \emph{alternative experiment}
#' (\code{\link[SingleCellExperiment]{altExp}}) stored inside a
#' \code{\linkS4class{SingleCellExperiment}}.
#'
#' The Y-axis feature is typically received from an
#' \code{\linkS4class{AltExpRowDataTable}} panel via \code{YAxisFeatureSource},
#' so that clicking a row in the table immediately plots the corresponding
#' altExp feature across all samples.
#'
#' @slot AltExp \code{character(1)}.  Name of the alternative experiment to
#'   visualise (must be present in \code{altExpNames(se)}).
#'   The \code{Assay} slot (inherited from \code{\linkS4class{FeatureAssayPlot}})
#'   refers to an assay within this altExp.
#'
#' @section Inherited slots:
#' All slots from \code{\linkS4class{FeatureAssayPlot}} are inherited,
#' including \code{Assay}, \code{XAxis}, \code{XAxisColumnData},
#' \code{YAxisFeatureName}, \code{YAxisFeatureSource}, and those controlling
#' colour, shape, size, faceting, and multi-selection behaviour.
#'
#' @seealso
#' \code{\link{AltExpFeatureAssayPlot}} for the constructor.
#' \code{\linkS4class{FeatureAssayPlot}} for the parent class.
#' \code{\linkS4class{AltExpRowDataTable}} for a compatible row-selection source.
#'
#' @name AltExpFeatureAssayPlot-class
#' @rdname AltExpFeatureAssayPlot-class
#' @exportClass AltExpFeatureAssayPlot
NULL

setClass("AltExpFeatureAssayPlot",
  contains = "FeatureAssayPlot",
  slots = c(AltExp = "character"))

#' @importFrom iSEE .emptyDefault
setMethod("initialize", "AltExpFeatureAssayPlot", function(.Object, ...) {
  args <- list(...)
  args <- .emptyDefault(args, "AltExp", NA_character_)
  do.call(callNextMethod, c(list(.Object), args))
})

#' Construct an AltExpFeatureAssayPlot panel
#'
#' Creates an instance of \code{\linkS4class{AltExpFeatureAssayPlot}} for use
#' as a panel in an iSEE application.
#'
#' @param ... Named arguments corresponding to slots of
#'   \code{\linkS4class{AltExpFeatureAssayPlot}} or its parent classes.
#'
#' @return An \code{\linkS4class{AltExpFeatureAssayPlot}} object.
#'
#' @seealso \code{\linkS4class{AltExpFeatureAssayPlot}} for slot details.
#' @export
AltExpFeatureAssayPlot <- function(...) new("AltExpFeatureAssayPlot", ...)

setMethod(".fullName",   "AltExpFeatureAssayPlot", function(x) "AltExp feature assay plot")
setMethod(".panelColor", "AltExpFeatureAssayPlot", function(x) "#AA5500")

#' @importFrom iSEE .getCachedCommonInfo .setCachedCommonInfo
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames
setMethod(".cacheCommonInfo", "AltExpFeatureAssayPlot", function(x, se) {
  if (!is.null(.getCachedCommonInfo(se, "AltExpFeatureAssayPlot"))) return(se)
  se <- callNextMethod()

  valid_assays_by_ae <- lapply(altExpNames(se), function(ae_name) {
    assayNames(altExp(se, ae_name))
  })
  names(valid_assays_by_ae) <- altExpNames(se)

  .setCachedCommonInfo(se, "AltExpFeatureAssayPlot",
    valid.altExp.names         = altExpNames(se),
    valid.assay.names.by.altExp = valid_assays_by_ae)
})

#' @importFrom iSEE .replaceMissingWithFirst
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames
setMethod(".refineParameters", "AltExpFeatureAssayPlot", function(x, se) {
  # Save before callNextMethod() validates against the main experiment rownames
  # and resets any altExp feature name that was transmitted from a source panel.
  saved_y <- slot(x, iSEE:::.featAssayYAxisFeatName)
  saved_x <- slot(x, iSEE:::.featAssayXAxisFeatName)

  x <- callNextMethod()
  if (is.null(x)) return(NULL)

  x <- .replaceMissingWithFirst(x, "AltExp", altExpNames(se))
  ae <- altExp(se, slot(x, "AltExp"))

  if (nrow(ae) == 0L) {
    warning(sprintf("no rows in altExp '%s'", slot(x, "AltExp")))
    return(NULL)
  }

  all_ae_assays <- assayNames(ae)
  if (length(all_ae_assays) == 0L) {
    warning(sprintf("no assays in altExp '%s'", slot(x, "AltExp")))
    return(NULL)
  }

  if (!slot(x, iSEE:::.featAssayAssay) %in% all_ae_assays) {
    slot(x, iSEE:::.featAssayAssay) <- all_ae_assays[1]
  }

  ae_rownames <- rownames(ae)

  # Restore saved feature names if they are valid altExp features; otherwise
  # fall back to the first altExp feature.  This preserves selections
  # transmitted from AltExpRowDataTable (altExp rownames) that callNextMethod()
  # would otherwise discard because they are absent from rownames(se).
  slot(x, iSEE:::.featAssayYAxisFeatName) <-
    if (!is.na(saved_y) && nzchar(saved_y) && saved_y %in% ae_rownames) {
      saved_y
    } else {
      ae_rownames[1]
    }

  slot(x, iSEE:::.featAssayXAxisFeatName) <-
    if (!is.na(saved_x) && nzchar(saved_x) && saved_x %in% ae_rownames) {
      saved_x
    } else {
      ae_rownames[1]
    }

  x
})

#' @importFrom iSEE .getEncodedName .getCachedCommonInfo .selectizeInput.iSEE
#'   .radioButtons.iSEE .conditionalOnRadio .addSpecificTour
#' @importFrom SingleCellExperiment altExpNames altExp
#' @importFrom SummarizedExperiment assayNames
#' @importFrom shiny selectInput checkboxInput selectizeInput
setMethod(".defineDataInterface", "AltExpFeatureAssayPlot", function(x, se, select_info) {
  panel_name <- .getEncodedName(x)
  .input_FUN <- function(field) paste0(panel_name, "_", field)

  all_altexps     <- altExpNames(se)
  current_ae      <- slot(x, "AltExp")
  all_ae_assays   <- assayNames(altExp(se, current_ae))
  column_covariates <- .getCachedCommonInfo(se, "ColumnDotPlot")$valid.colData.names
  tab_by_row      <- select_info$single$feature

  xaxis_choices <- iSEE:::.featAssayXAxisNothingTitle
  if (length(column_covariates)) {
    xaxis_choices <- c(xaxis_choices, iSEE:::.featAssayXAxisColDataTitle)
  }
  xaxis_choices <- c(xaxis_choices, iSEE:::.featAssayXAxisFeatNameTitle)

  .addSpecificTour(class(x)[1], "AltExp", function(plot_name) {
    data.frame(rbind(
      c(element = paste0("#", plot_name, "_AltExp + .selectize-control"),
        intro = "Select the alternative experiment whose assay values are
plotted.  Changing this updates the available assays and — when driven by
an <em>AltExp row data table</em> — the available features."),
      c(element = paste0("#", plot_name, "_", iSEE:::.featAssayAssay,
                         " + .selectize-control"),
        intro = "Select the assay within the chosen alternative experiment
to display on the y-axis — for example <code>logcounts</code> for
log-normalised values.")
    ), stringsAsFactors = FALSE)
  })

  list(
    # AltExp selector
    selectInput(.input_FUN("AltExp"),
      label    = "AltExp:",
      choices  = all_altexps,
      selected = iSEE:::.choose_link(current_ae, all_altexps)
    ),
    # Y-axis feature (driven by AltExpRowDataTable via YAxisFeatureSource)
    .selectizeInput.iSEE(x, iSEE:::.featAssayYAxisFeatName,
      label    = "Y-axis feature:",
      choices  = NULL, selected = NULL, multiple = FALSE
    ),
    selectInput(.input_FUN(iSEE:::.featAssayYAxisRowTable),
      label    = NULL,
      choices  = tab_by_row,
      selected = iSEE:::.choose_link(slot(x, iSEE:::.featAssayYAxisRowTable),
                                     tab_by_row)
    ),
    checkboxInput(.input_FUN(iSEE:::.featAssayYAxisFeatDynamic),
      label = "Use dynamic feature selection for the y-axis",
      value = slot(x, iSEE:::.featAssayYAxisFeatDynamic)
    ),
    # Assay from the altExp
    selectInput(.input_FUN(iSEE:::.featAssayAssay),
      label    = "Assay:",
      choices  = all_ae_assays,
      selected = iSEE:::.choose_link(slot(x, iSEE:::.featAssayAssay),
                                     all_ae_assays)
    ),
    # X-axis
    .radioButtons.iSEE(x, iSEE:::.featAssayXAxis,
      label    = "X-axis:",
      inline   = TRUE,
      choices  = xaxis_choices,
      selected = slot(x, iSEE:::.featAssayXAxis)
    ),
    .conditionalOnRadio(.input_FUN(iSEE:::.featAssayXAxis),
      iSEE:::.featAssayXAxisColDataTitle,
      selectInput(.input_FUN(iSEE:::.featAssayXAxisColData),
        label    = "X-axis column data:",
        choices  = column_covariates,
        selected = slot(x, iSEE:::.featAssayXAxisColData)
      )
    ),
    .conditionalOnRadio(.input_FUN(iSEE:::.featAssayXAxis),
      iSEE:::.featAssayXAxisFeatNameTitle,
      selectizeInput(.input_FUN(iSEE:::.featAssayXAxisFeatName),
        label = "X-axis feature:", choices = NULL, selected = NULL,
        multiple = FALSE
      ),
      selectInput(.input_FUN(iSEE:::.featAssayXAxisRowTable),
        label    = NULL,
        choices  = tab_by_row,
        selected = slot(x, iSEE:::.featAssayXAxisRowTable)
      ),
      checkboxInput(.input_FUN(iSEE:::.featAssayXAxisFeatDynamic),
        label = "Use dynamic feature selection for the x-axis",
        value = slot(x, iSEE:::.featAssayXAxisFeatDynamic)
      )
    )
  )
})

#' @importFrom iSEE .getEncodedName .createProtectedParameterObservers
#'   .trackSingleSelection .trackRelinkedSelection
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment assayNames
#' @importFrom shiny isolate observe observeEvent updateSelectInput updateSelectizeInput
setMethod(".createObservers", "AltExpFeatureAssayPlot",
    function(x, se, input, session, pObjects, rObjects) {
  callNextMethod()

  plot_name    <- .getEncodedName(x)
  ae_field     <- paste0(plot_name, "_AltExp")
  assay_field  <- paste0(plot_name, "_", iSEE:::.featAssayAssay)
  feat_field   <- paste0(plot_name, "_", iSEE:::.featAssayYAxisFeatName)
  source_field <- paste0(plot_name, "_", iSEE:::.featAssayYAxisRowTable)

  # Helper: repopulate the Y-axis feature selectize with altExp rownames.
  # Called via session$onFlushed so it fires AFTER all reactive observers in
  # the current flush have run — including FeatureAssayPlot's own observer that
  # calls updateSelectizeInput(choices = rownames(se)).  Because Shiny sends
  # queued messages in order, our message arrives at the browser after iSEE's
  # and overwrites it with the correct altExp choices.
  repopulate_feat_choices <- function() {
    current_ae  <- pObjects$memory[[plot_name]][["AltExp"]]
    ae_rownames <- rownames(altExp(se, current_ae))
    updateSelectizeInput(session, feat_field,
      choices  = ae_rownames,
      selected = pObjects$memory[[plot_name]][[iSEE:::.featAssayYAxisFeatName]],
      server   = TRUE
    )
  }

  # Initial population on app startup
  session$onFlushed(repopulate_feat_choices, once = TRUE)

  # Intercept single-selection from the source panel and write directly into
  # pObjects$memory, bypassing iSEE's updateSelectizeInput(choices=rownames(se))
  # gate which silently drops altExp feature names not in the main experiment.
  observe({
    source_panel <- input[[source_field]]
    if (!length(source_panel) || !nzchar(source_panel) ||
        identical(source_panel, iSEE:::.noSelection)) return()

    iSEE:::.trackSingleSelection(source_panel, rObjects)

    isolate({
      new_feat <- iSEE:::.singleSelectionValue(
        pObjects$memory[[source_panel]],
        pObjects$contents[[source_panel]]
      )
      if (is.null(new_feat) || !nzchar(new_feat)) return()

      current_ae <- pObjects$memory[[plot_name]][["AltExp"]]
      ae_rn      <- rownames(altExp(se, current_ae))
      if (!new_feat %in% ae_rn) return()

      old_feat <- pObjects$memory[[plot_name]][[iSEE:::.featAssayYAxisFeatName]]
      if (identical(new_feat, old_feat)) return()

      pObjects$memory[[plot_name]][[iSEE:::.featAssayYAxisFeatName]] <- new_feat
      updateSelectizeInput(session, feat_field,
        choices  = ae_rn,
        selected = new_feat,
        server   = TRUE
      )
      iSEE:::.requestCleanUpdate(plot_name, pObjects, rObjects)
    })
  })

  # Re-population whenever the source panel transmits a new feature OR the
  # source panel itself is changed (relink).  We observe the same reactives
  # that FeatureAssayPlot's inherited observer watches, then schedule our fix
  # via onFlushed so it runs AFTER that observer within the same flush.
  observe({
    source_panel <- input[[source_field]]
    if (length(source_panel) && nzchar(source_panel)) {
      .trackSingleSelection(source_panel, rObjects)
    }
    .trackRelinkedSelection(plot_name, rObjects)
    session$onFlushed(repopulate_feat_choices, once = TRUE)
  }, priority = -1L)

  # nocov start
  observeEvent(input[[ae_field]], {
    matched_ae <- input[[ae_field]]
    if (identical(matched_ae, pObjects$memory[[plot_name]][["AltExp"]])) return(NULL)
    pObjects$memory[[plot_name]][["AltExp"]] <- matched_ae

    new_assays <- assayNames(altExp(se, matched_ae))
    cur_assay  <- pObjects$memory[[plot_name]][[iSEE:::.featAssayAssay]]
    new_assay  <- if (cur_assay %in% new_assays) cur_assay else new_assays[1]
    pObjects$memory[[plot_name]][[iSEE:::.featAssayAssay]] <- new_assay

    new_rn   <- rownames(altExp(se, matched_ae))
    cur_feat <- pObjects$memory[[plot_name]][[iSEE:::.featAssayYAxisFeatName]]
    new_feat <- if (!is.na(cur_feat) && nzchar(cur_feat) && cur_feat %in% new_rn) {
      cur_feat
    } else {
      new_rn[1]
    }
    pObjects$memory[[plot_name]][[iSEE:::.featAssayYAxisFeatName]] <- new_feat

    updateSelectInput(session, assay_field,
      choices  = new_assays,
      selected = new_assay)
    updateSelectizeInput(session, feat_field,
      choices  = new_rn,
      selected = new_feat,
      server   = TRUE)

    iSEE:::.requestCleanUpdate(plot_name, pObjects, rObjects)
  }, ignoreInit = TRUE)
  # nocov end

  invisible(NULL)
})

#' @importFrom iSEE .getEncodedName .getPanelColor .addTourStep .dataParamBoxOpen
setMethod(".definePanelTour", "AltExpFeatureAssayPlot", function(x) {
  collated <- rbind(
    c(
      element = paste0("#", .getEncodedName(x)),
      intro   = sprintf(
        "The <font color=\"%s\">AltExp feature assay plot</font> shows assay
values for a single feature from an <em>alternative experiment</em> stored
inside a <code>SingleCellExperiment</code>.  Each point is a sample; the
y-axis shows the assay value of the selected altExp feature.
<br><br>
The feature is typically chosen by clicking a row in an
<em>AltExp row data table</em> panel linked via the
<em>Y-axis feature source</em>.", .getPanelColor(x))
    ),
    .addTourStep(x, .dataParamBoxOpen,
      "The <i>Data parameters</i> box controls which alternative experiment,
assay, and feature to display, as well as the x-axis variable.
<br><br><strong>Action:</strong> click to expand.")
  )

  parent_tour <- callNextMethod()

  # Drop the generic FeatureAssayPlot opening step — ours replaces it
  parent_tour <- parent_tour[
    !grepl("Feature assay plot", parent_tour$intro, fixed = TRUE), ]

  rbind(
    data.frame(element = collated[, 1], intro = collated[, 2],
               stringsAsFactors = FALSE),
    parent_tour
  )
})

#' @importFrom iSEE .textEval
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment assay colData
setMethod(".generateDotPlotData", "AltExpFeatureAssayPlot", function(x, envir) {
  data_cmds <- list()

  ae_name         <- slot(x, "AltExp")
  gene_selected_y <- slot(x, iSEE:::.featAssayYAxisFeatName)
  assay_choice    <- slot(x, iSEE:::.featAssayAssay)

  plot_title <- gene_selected_y
  y_lab      <- sprintf("%s (%s)", gene_selected_y, assay_choice)

  err_msg <- "Select another altExp with the correct assay"
  make_error_result <- function() {
    err_cmd <- paste0(
      "plot.data <- data.frame(X = rep(NA_real_, ncol(se)),",
      " Y = rep(NA_real_, ncol(se)), row.names = colnames(se));")
    .textEval(err_cmd, envir)
    list(commands = err_cmd, labels = list(title = err_msg, X = "", Y = ""))
  }

  data_cmds[["y"]] <- c(
    sprintf("ae <- altExp(se, %s);", deparse(ae_name)),
    sprintf(
      "plot.data <- data.frame(Y = assay(ae, %s)[%s, ], row.names = colnames(se));",
      deparse(assay_choice), deparse(gene_selected_y))
  )

  x_choice <- slot(x, iSEE:::.featAssayXAxis)

  if (x_choice == iSEE:::.featAssayXAxisColDataTitle) {
    x_lab      <- slot(x, iSEE:::.featAssayXAxisColData)
    plot_title <- paste(plot_title, "vs", x_lab)
    data_cmds[["x"]] <- sprintf(
      "plot.data$X <- colData(se)[, %s];", deparse(x_lab))

  } else if (x_choice == iSEE:::.featAssayXAxisFeatNameTitle) {
    gene_selected_x <- slot(x, iSEE:::.featAssayXAxisFeatName)
    plot_title      <- paste(plot_title, "vs", gene_selected_x)
    x_lab           <- sprintf("%s (%s)", gene_selected_x, assay_choice)
    data_cmds[["x"]] <- sprintf(
      "plot.data$X <- assay(ae, %s)[%s, ];",
      deparse(assay_choice), deparse(gene_selected_x))

  } else if (x_choice == iSEE:::.featAssayXAxisSelectionsTitle) {
    x_lab      <- "Column selection"
    plot_title <- paste(plot_title, "vs column selection")
    target <- if (exists("col_selected", envir = envir, inherits = FALSE)) {
      "col_selected"
    } else {
      "list()"
    }
    data_cmds[["x"]] <- sprintf(
      "plot.data$X <- iSEE::multiSelectionToFactor(%s, colnames(se));", target)

  } else {
    x_lab            <- ""
    data_cmds[["x"]] <- "plot.data$X <- factor(character(ncol(se)))"
  }

  data_cmds <- unlist(data_cmds)
  ok <- tryCatch({ .textEval(data_cmds, envir); TRUE }, error = function(e) FALSE)
  if (!ok) return(make_error_result())

  list(commands = data_cmds,
       labels   = list(title = plot_title, X = x_lab, Y = y_lab))
})

#' @importFrom iSEE .generateDotPlot .textEval
setMethod(".generateDotPlot", "AltExpFeatureAssayPlot", function(x, labels, envir) {
  if (identical(labels$Y, "")) {
    plot_cmds <- c(
      "dot.plot <- ggplot(plot.data, aes(x = X, y = Y)) +",
      "  scale_x_continuous(limits = c(0, 1)) +",
      "  scale_y_continuous(limits = c(0, 1)) +",
      sprintf("  labs(title = %s, x = NULL, y = NULL) +", deparse(labels$title)),
      "  theme_bw() +",
      "  theme(",
      "    plot.title = element_text(hjust = 0.5, color = 'red', size = 12),",
      "    axis.text  = element_blank(),",
      "    axis.ticks = element_blank()",
      "  );"
    )
    .textEval(plot_cmds, envir)
    return(list(plot = envir$dot.plot, commands = plot_cmds))
  }
  callNextMethod()
})
