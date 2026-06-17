############################################################
# AltVolcanoPlot
############################################################

# Extends VolcanoPlot from iSEEu to display rowData fields from an alternative
# experiment or the main experiment (via the "(Main)" sentinel).
# The key override is in .generateDotPlotData, which substitutes altExp(se, ae)
# or se for all rowData lookups.

#' AltVolcanoPlot: volcano plot for alternative experiments
#'
#' An S4 class extending \code{\linkS4class{VolcanoPlot}} from
#' \pkg{iSEEu} to display a volcano plot for features of an
#' \emph{alternative experiment}
#' (\code{\link[SingleCellExperiment]{altExp}}) stored inside a
#' \code{\linkS4class{SingleCellExperiment}}.
#'
#' The Y-axis shows \eqn{-\log_{10}(p\text{-value})} and the X-axis shows the
#' log fold-change, both taken from \code{rowData} of the selected altExp.
#' Points are coloured by significance status using the inherited
#' \code{PValueThreshold}, \code{LogFCThreshold}, and
#' \code{PValueCorrection} slots.
#'
#' @slot Experiment \code{character(1)}.  Name of the alternative experiment whose
#'   \code{rowData} supplies the p-value and log fold-change columns, or the
#'   sentinel \code{"(Main)"} to use the main experiment's
#'   \code{rowData}.  Defaults to \code{NA}, resolved at runtime to
#'   the first available alternative experiment.
#'
#' @section Inherited slots:
#' All slots from \code{\linkS4class{VolcanoPlot}} are inherited, including
#' \code{YAxis} (p-value column), \code{XAxisRowData} (logFC column),
#' \code{PValueThreshold}, \code{LogFCThreshold}, and
#' \code{PValueCorrection}.
#'
#' @seealso
#' \code{\link{AltVolcanoPlot}} for the constructor.
#' \code{\linkS4class{VolcanoPlot}} for the parent class in \pkg{iSEEu}.
#' \code{\linkS4class{AltRowDataTable}} for a compatible row-selection source.
#'
#' @name AltVolcanoPlot-class
#' @rdname AltVolcanoPlot-class
#' @exportClass AltVolcanoPlot
NULL

#' @importClassesFrom iSEEu VolcanoPlot
setClass("AltVolcanoPlot",
  contains = "VolcanoPlot",
  slots = c(Experiment = "character"))

#' @importFrom iSEE .emptyDefault
setMethod("initialize", "AltVolcanoPlot", function(.Object, ...) {
  args <- list(...)
  args <- .emptyDefault(args, "Experiment", NA_character_)
  do.call(callNextMethod, c(list(.Object), args))
})

#' Construct an AltVolcanoPlot panel
#'
#' Creates an instance of \code{\linkS4class{AltVolcanoPlot}} for use as a
#' panel in an iSEE application.
#'
#' @param ... Named arguments corresponding to slots of
#'   \code{\linkS4class{AltVolcanoPlot}} or its parent classes.
#'
#' @return An \code{\linkS4class{AltVolcanoPlot}} object.
#'
#' @seealso \code{\linkS4class{AltVolcanoPlot}} for slot details.
#' @export
AltVolcanoPlot <- function(...) new("AltVolcanoPlot", ...)

setMethod(".fullName",   "AltVolcanoPlot", function(x) "Alt volcano plot")
setMethod(".panelColor", "AltVolcanoPlot", function(x) "#6B2D8B")

#' @importFrom iSEE .getCachedCommonInfo .setCachedCommonInfo .findAtomicFields
#'   .whichNumeric
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment rowData
setMethod(".cacheCommonInfo", "AltVolcanoPlot", function(x, se) {
  if (!is.null(.getCachedCommonInfo(se, "AltVolcanoPlot"))) return(se)
  se <- callNextMethod()

  all_names <- c(altExpNames(se), .selectionMainExpTitle)

  valid_p_by_ae       <- vector("list", length(all_names))
  valid_lfc_by_ae     <- vector("list", length(all_names))
  valid_rd_by_ae      <- vector("list", length(all_names))
  valid_numeric_by_ae <- vector("list", length(all_names))
  names(valid_p_by_ae) <- names(valid_lfc_by_ae) <-
    names(valid_rd_by_ae) <- names(valid_numeric_by_ae) <- all_names

  for (ae_name in all_names) {
    ae_obj      <- if (identical(ae_name, .selectionMainExpTitle)) se else altExp(se, ae_name)
    rd          <- rowData(ae_obj)
    displayable <- .findAtomicFields(rd)
    continuous  <- displayable[.whichNumeric(rd[, displayable, drop = FALSE])]

    valid_rd_by_ae[[ae_name]]      <- displayable
    valid_numeric_by_ae[[ae_name]] <- continuous
    valid_p_by_ae[[ae_name]]       <- iSEEu:::.matchPValueFields(se, continuous)
    valid_lfc_by_ae[[ae_name]]     <- iSEEu:::.matchLogFCFields(se, continuous)
  }

  .setCachedCommonInfo(se, "AltVolcanoPlot",
    valid.altExp.names           = all_names,
    valid.rowData.by.altExp      = valid_rd_by_ae,
    valid.numeric.by.altExp      = valid_numeric_by_ae,
    valid.p.fields.by.altExp     = valid_p_by_ae,
    valid.lfc.fields.by.altExp   = valid_lfc_by_ae)
})

#' @importFrom iSEE .replaceMissingWithFirst
#' @importFrom SingleCellExperiment altExpNames
setMethod(".refineParameters", "AltVolcanoPlot", function(x, se) {
  x <- callNextMethod()
  if (is.null(x)) return(NULL)

  x <- .replaceMissingWithFirst(x, "Experiment", c(altExpNames(se), .selectionMainExpTitle))
  ae_name <- slot(x, "Experiment")

  cache <- .getCachedCommonInfo(se, "AltVolcanoPlot")

  valid_p       <- cache$valid.p.fields.by.altExp[[ae_name]]
  valid_lfc     <- cache$valid.lfc.fields.by.altExp[[ae_name]]
  valid_rd      <- cache$valid.rowData.by.altExp[[ae_name]]
  valid_numeric <- cache$valid.numeric.by.altExp[[ae_name]]

  # Require named p-value AND logFC columns — fall back to any numeric field is
  # intentionally NOT done here.  Without both, the volcano math (abs, -log10)
  # is meaningless, so we mark the panel invalid and let generateDotPlotData
  # produce a readable error plot.
  if (length(valid_p) == 0L || length(valid_lfc) == 0L) {
    slot(x, iSEE:::.rowDataYAxis) <- NA_character_
    return(x)
  }

  x <- .replaceMissingWithFirst(x, iSEE:::.rowDataYAxis,        valid_p)
  x <- .replaceMissingWithFirst(x, iSEE:::.rowDataXAxisRowData, valid_lfc)
  x[["XAxis"]] <- "Row data"
  x
})

#' @importFrom iSEE .getEncodedName .getCachedCommonInfo .selectInput.iSEE
#'   .addSpecificTour
#' @importFrom SingleCellExperiment altExpNames
#' @importFrom shiny selectInput hr
setMethod(".defineDataInterface", "AltVolcanoPlot", function(x, se, select_info) {
  panel_name <- .getEncodedName(x)

  all_altexps <- c(altExpNames(se), .selectionMainExpTitle)
  ae_name       <- slot(x, "Experiment")
  cache         <- .getCachedCommonInfo(se, "AltVolcanoPlot")
  valid_p       <- cache$valid.p.fields.by.altExp[[ae_name]]
  valid_lfc     <- cache$valid.lfc.fields.by.altExp[[ae_name]]
  valid_rd      <- cache$valid.rowData.by.altExp[[ae_name]]
  valid_numeric <- cache$valid.numeric.by.altExp[[ae_name]]

  if (length(valid_p) == 0L || length(valid_lfc) == 0L) {
    p_choices   <- character(0)
    lfc_choices <- character(0)
  } else {
    p_choices   <- valid_p
    lfc_choices <- valid_lfc
  }

  .addSpecificTour(class(x)[1], "Experiment", function(plot_name) {
    data.frame(rbind(
      c(element = paste0("#", plot_name, "_Experiment + .selectize-control"),
        intro = "Select the alternative experiment whose <code>rowData</code>
supplies the p-values and log fold-changes for the volcano plot."),
      c(element = paste0("#", plot_name, "_YAxis + .selectize-control"),
        intro = "Choose the <code>rowData</code> column of the selected
alternative experiment that contains the raw (uncorrected, untransformed)
p-values.  The panel takes care of the −log<sub>10</sub> transform."),
      c(element = paste0("#", plot_name, "_XAxisRowData + .selectize-control"),
        intro = "Choose the <code>rowData</code> column of the selected
alternative experiment that contains the log fold-changes.")
    ), stringsAsFactors = FALSE)
  })

  c(
    list(
      selectInput(
        paste0(panel_name, "_Experiment"),
        label    = "Experiment:",
        choices  = all_altexps,
        selected = iSEE:::.choose_link(ae_name, all_altexps)
      ),
      .selectInput.iSEE(x, iSEE:::.rowDataYAxis,
        label    = "P-value field:",
        choices  = p_choices,
        selected = slot(x, iSEE:::.rowDataYAxis)
      ),
      .selectInput.iSEE(x, iSEE:::.rowDataXAxisRowData,
        label    = "Log FC field:",
        choices  = lfc_choices,
        selected = slot(x, iSEE:::.rowDataXAxisRowData)
      )
    ),
    list(hr()),
    iSEEu:::.define_gene_sig_ui(x)
  )
})

#' @importFrom iSEE .getEncodedName
#' @importFrom SingleCellExperiment altExpNames altExp
#' @importFrom SummarizedExperiment rowData
#' @importFrom shiny observeEvent updateSelectInput
setMethod(".createObservers", "AltVolcanoPlot",
    function(x, se, input, session, pObjects, rObjects) {
  callNextMethod()

  plot_name <- .getEncodedName(x)
  ae_field  <- paste0(plot_name, "_Experiment")
  y_field   <- paste0(plot_name, "_", iSEE:::.rowDataYAxis)
  x_field   <- paste0(plot_name, "_", iSEE:::.rowDataXAxisRowData)

  # nocov start
  # Note: we deliberately do NOT call .createUnprotectedParameterObservers for
  # "Experiment" here.  That helper would update pObjects$memory[["Experiment"]] before
  # this observeEvent fires, causing the identical() guard to short-circuit and
  # leaving the p-value / logFC dropdowns stale.  We own the full update cycle.
  observeEvent(input[[ae_field]], {
    matched_ae <- input[[ae_field]]
    if (identical(matched_ae, pObjects$memory[[plot_name]][["Experiment"]])) return(NULL)
    pObjects$memory[[plot_name]][["Experiment"]] <- matched_ae

    cache         <- .getCachedCommonInfo(se, "AltVolcanoPlot")
    valid_p       <- cache$valid.p.fields.by.altExp[[matched_ae]]
    valid_lfc     <- cache$valid.lfc.fields.by.altExp[[matched_ae]]
    valid_rd      <- cache$valid.rowData.by.altExp[[matched_ae]]
    valid_numeric <- cache$valid.numeric.by.altExp[[matched_ae]]

    if (length(valid_p) == 0L || length(valid_lfc) == 0L) {
      pObjects$memory[[plot_name]][[iSEE:::.rowDataYAxis]] <- NA_character_
      updateSelectInput(session, y_field, choices = character(0), selected = NULL)
      updateSelectInput(session, x_field, choices = character(0), selected = NULL)
      iSEE:::.requestCleanUpdate(plot_name, pObjects, rObjects)
      return(NULL)
    }

    p_choices   <- valid_p
    lfc_choices <- valid_lfc

    cur_y <- pObjects$memory[[plot_name]][[iSEE:::.rowDataYAxis]]
    new_y <- if (!is.na(cur_y) && cur_y %in% p_choices) cur_y else p_choices[1]
    pObjects$memory[[plot_name]][[iSEE:::.rowDataYAxis]] <- new_y

    cur_x <- pObjects$memory[[plot_name]][[iSEE:::.rowDataXAxisRowData]]
    new_x <- if (cur_x %in% lfc_choices) cur_x else lfc_choices[1]
    pObjects$memory[[plot_name]][[iSEE:::.rowDataXAxisRowData]] <- new_x

    updateSelectInput(session, y_field, choices = p_choices,   selected = new_y)
    updateSelectInput(session, x_field, choices = lfc_choices, selected = new_x)
    iSEE:::.requestCleanUpdate(plot_name, pObjects, rObjects)
  }, ignoreInit = TRUE)
  # nocov end

  invisible(NULL)
})

#' @importFrom iSEE .textEval
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment rowData
setMethod(".generateDotPlotData", "AltVolcanoPlot", function(x, envir) {
  ae_name <- slot(x, "Experiment")
  y_lab   <- slot(x, iSEE:::.rowDataYAxis)
  x_lab   <- slot(x, iSEE:::.rowDataXAxisRowData)

  err_msg <- "Select an experiment with p-values and log fold changes in the rowData"

  # Helper: populate envir with a safe NA-filled plot.data and return the error
  # result list.  nrow(se) rows keeps downstream colour/shape/size methods from
  # failing on a length mismatch when they attach rowData(se) columns.
  # IsSig is included because the parent VolcanoPlot renderer expects it.
  make_error_result <- function() {
    cmds <- paste0(
      "plot.data <- data.frame(",
      "X = rep(NA_real_, nrow(se)), Y = rep(NA_real_, nrow(se)),",
      " IsSig = rep('none', nrow(se)),",
      " row.names = rownames(se));"
    )
    .textEval(cmds, envir)
    list(commands = cmds, labels = list(title = err_msg, X = "", Y = ""))
  }

  # Fast path: .refineParameters already detected no valid fields and set the
  # NA sentinel.  Also guards against NA_character_ reaching rowData()[, NA].
  if (is.na(y_lab)) return(make_error_result())

  data_cmds <- c(
    if (identical(ae_name, .selectionMainExpTitle)) "ae <- se;" else
      sprintf("ae <- altExp(se, %s);", deparse(ae_name)),
    sprintf("plot.data <- data.frame(Y = rowData(ae)[, %s], row.names = rownames(ae));",
            deparse(y_lab)),
    sprintf("plot.data$X <- rowData(ae)[, %s];", deparse(x_lab))
  )

  extra_cmds <- c(
    iSEEu:::.define_de_status(x, lfc = "plot.data$X", pval = "plot.data$Y"),
    "plot.data$IsSig <- c('down', 'none', 'up')[.de_status];",
    "plot.data$Y <- -log10(plot.data$Y)"
  )

  all_cmds <- c(data_cmds, extra_cmds)

  # Defensive catch: during AltExp transitions pObjects$memory may still carry
  # a stale YAxis field from the previous altExp before .refineParameters has
  # had a chance to run.  Any subscript / non-numeric error is caught here so
  # the panel shows the guidance message rather than freezing the app.
  ok <- tryCatch({ .textEval(all_cmds, envir); TRUE }, error = function(e) FALSE)
  if (!ok) return(make_error_result())

  list(
    commands = all_cmds,
    labels   = list(
      title = sprintf("%s %s vs %s", ae_name, y_lab, x_lab),
      X     = x_lab,
      Y     = sprintf("-Log10[%s]", y_lab)
    )
  )
})

#' @importFrom iSEE .generateDotPlot .textEval
setMethod(".generateDotPlot", "AltVolcanoPlot", function(x, labels, envir) {
  # Two conditions both indicate the error state:
  #   1. The NA sentinel set by .refineParameters (normal path)
  #   2. Empty X label — set by make_error_result() in .generateDotPlotData
  #      when a stale field name caused .textEval to fail mid-transition.
  if (is.na(slot(x, iSEE:::.rowDataYAxis)) || identical(labels$X, "")) {
    # Use plot.data with aes(x = X, y = Y) so coordinfo$mapping$x is "X".
    # Without this, Shiny's nearPoints() cannot infer xvar and freezes the panel.
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

#' @importFrom iSEE .getEncodedName .getPanelColor .addTourStep .dataParamBoxOpen
setMethod(".definePanelTour", "AltVolcanoPlot", function(x) {
  collated <- rbind(
    c(
      element = paste0("#", .getEncodedName(x)),
      intro   = sprintf(
        "The <font color=\"%s\">Alt volcano plot</font> shows the
−log<sub>10</sub>(p-value) against the log fold-change for every feature in
the selected alternative experiment.  Points are coloured <em>down</em> /
<em>none</em> / <em>up</em> based on the configured significance thresholds,
just like the standard iSEEu Volcano plot, but using <code>rowData</code>
from the altExp instead of the main experiment.",
        .getPanelColor(x))
    ),
    .addTourStep(x, .dataParamBoxOpen,
      "The <i>Data parameters</i> box lets you choose the alternative
experiment and the <code>rowData</code> columns to use for the p-value and
log fold-change axes.<br><br><strong>Action:</strong> click to expand.")
  )

  parent_tour <- callNextMethod()
  skip <- grep("VisualBoxOpen$", parent_tour$element)
  if (length(skip) > 0L) parent_tour <- parent_tour[-seq_len(skip - 1L), ]

  rbind(
    data.frame(element = collated[, 1], intro = collated[, 2],
               stringsAsFactors = FALSE),
    parent_tour
  )
})
