############################################################
# AltReducedDimensionPlot
############################################################

# Extends ColumnDotPlot (not ReducedDimensionPlot) to avoid observer conflicts:
# ReducedDimensionPlot's Type observer calls reducedDim(se, ...) which would
# fail for altExp-only reducedDim names. We reimplement the reducedDim logic
# from scratch, substituting altExp(se, ae) for se throughout.

#' AltReducedDimensionPlot: reduced dimension plot for alternative experiments
#'
#' An S4 class extending \code{\linkS4class{ColumnDotPlot}} to display
#' dimensionality reduction results stored in the
#' \code{\link[SingleCellExperiment]{reducedDims}} of an
#' \emph{alternative experiment}
#' (\code{\link[SingleCellExperiment]{altExp}}) inside a
#' \code{\linkS4class{SingleCellExperiment}}.
#'
#' @slot Experiment \code{character(1)}.  Name of the alternative experiment
#'   (must be present in \code{altExpNames(se)}), or the sentinel
#'   \code{"(Main)"} to use the main experiment's reducedDims.
#' @slot Type \code{character(1)}.  Name of the dimensionality reduction result
#'   within the selected \code{AltExp} (must be present in
#'   \code{reducedDimNames(altExp(se, AltExp))}).  Defaults to \code{NA},
#'   resolved at runtime to the first available result.
#' @slot XAxis \code{integer(1)}.  Component to plot on the x-axis. Default 1.
#' @slot YAxis \code{integer(1)}.  Component to plot on the y-axis. Default 2.
#'
#' @section Inherited slots:
#' All slots from \code{\linkS4class{ColumnDotPlot}} are inherited, including
#' those controlling colour, shape, size, faceting, and multi-selection
#' behaviour.
#'
#' @seealso
#' \code{\link{AltReducedDimensionPlot}} for the constructor.
#' \code{\linkS4class{ColumnDotPlot}} for the parent class.
#' \code{\linkS4class{ReducedDimensionPlot}} for the analogous main-experiment panel.
#'
#' @name AltReducedDimensionPlot-class
#' @rdname AltReducedDimensionPlot-class
#' @exportClass AltReducedDimensionPlot
NULL

setClass("AltReducedDimensionPlot",
  contains = "ColumnDotPlot",
  slots = c(
    Experiment = "character",
    Type   = "character",
    XAxis  = "integer",
    YAxis  = "integer"
  ))

setValidity2("AltReducedDimensionPlot", function(object) {
  # Type is intentionally allowed to be NA at construction time: like AltExp,
  # it is resolved at runtime by .refineParameters (first available
  # reducedDim result), so it must not be enforced as a non-NA string here.
  msg <- character(0)
  for (field in c("XAxis", "YAxis")) {
    val <- slot(object, field)
    if (length(val) != 1L || is.na(val) || val <= 0L)
      msg <- c(msg, sprintf("'%s' must be a single positive integer", field))
  }
  if (length(msg)) return(msg)
  TRUE
})

#' @importFrom iSEE .emptyDefault
setMethod("initialize", "AltReducedDimensionPlot", function(.Object, ...) {
  args <- list(...)
  args <- .emptyDefault(args, "Experiment", NA_character_)
  args <- .emptyDefault(args, "Type",   NA_character_)
  args <- .emptyDefault(args, "XAxis",  1L)
  args <- .emptyDefault(args, "YAxis",  2L)
  do.call(callNextMethod, c(list(.Object), args))
})

#' Construct an AltReducedDimensionPlot panel
#'
#' Creates an instance of \code{\linkS4class{AltReducedDimensionPlot}} for
#' use as a panel in an iSEE application.
#'
#' @param ... Named arguments corresponding to slots of
#'   \code{\linkS4class{AltReducedDimensionPlot}} or its parent classes.
#'
#' @return An \code{\linkS4class{AltReducedDimensionPlot}} object.
#'
#' @seealso \code{\linkS4class{AltReducedDimensionPlot}} for slot details.
#' @export
AltReducedDimensionPlot <- function(...) new("AltReducedDimensionPlot", ...)

setMethod(".fullName",   "AltReducedDimensionPlot", function(x) "Alt reduced dimension plot")
setMethod(".panelColor", "AltReducedDimensionPlot", function(x) "#006666")

#' @importFrom iSEE .getCachedCommonInfo .setCachedCommonInfo
#' @importFrom SingleCellExperiment altExpNames
setMethod(".cacheCommonInfo", "AltReducedDimensionPlot", function(x, se) {
  if (!is.null(.getCachedCommonInfo(se, "AltReducedDimensionPlot"))) return(se)
  se <- callNextMethod()
  .setCachedCommonInfo(se, "AltReducedDimensionPlot",
    valid.altExp.names = c(altExpNames(se), .selectionMainExpTitle))
})

#' @importFrom iSEE .replaceMissingWithFirst
#' @importFrom SingleCellExperiment altExp altExpNames reducedDim reducedDimNames
#' @importFrom methods is
setMethod(".refineParameters", "AltReducedDimensionPlot", function(x, se) {
  x <- callNextMethod()
  if (is.null(x)) return(NULL)

  x <- .replaceMissingWithFirst(x, "Experiment", c(altExpNames(se), .selectionMainExpTitle))
  ae_name <- slot(x, "Experiment")
  ae <- if (identical(ae_name, .selectionMainExpTitle)) se else altExp(se, ae_name)

  if (!is(ae, "SingleCellExperiment")) {
    slot(x, "Type") <- NA_character_
    return(x)
  }

  available <- reducedDimNames(ae)
  available  <- available[vapply(available,
    function(n) ncol(reducedDim(ae, n)) > 0L, logical(1))]

  if (length(available) == 0L) {
    slot(x, "Type") <- NA_character_
    return(x)
  }

  chosen <- slot(x, "Type")
  if (!is.na(chosen) && chosen %in% available &&
      slot(x, "XAxis") <= ncol(reducedDim(ae, chosen)) &&
      slot(x, "YAxis") <= ncol(reducedDim(ae, chosen))) {
    # all valid, nothing to do
  } else {
    y <- available[1]
    slot(x, "Type")  <- y
    slot(x, "XAxis") <- 1L
    slot(x, "YAxis") <- min(ncol(reducedDim(ae, y)), 2L)
  }
  x
})

#' @importFrom iSEE .getEncodedName .getCachedCommonInfo .selectInput.iSEE .addSpecificTour
#' @importFrom SingleCellExperiment altExp altExpNames reducedDim reducedDimNames
#' @importFrom shiny selectInput
setMethod(".defineDataInterface", "AltReducedDimensionPlot", function(x, se, select_info) {
  panel_name <- .getEncodedName(x)
  .input_FUN <- function(field) paste0(panel_name, "_", field)

  all_altexps <- c(altExpNames(se), .selectionMainExpTitle)
  current_ae  <- slot(x, "Experiment")
  ae <- if (identical(current_ae, .selectionMainExpTitle)) se else altExp(se, current_ae)

  if (is(ae, "SingleCellExperiment")) {
    available <- reducedDimNames(ae)
    available  <- available[vapply(available,
      function(n) ncol(reducedDim(ae, n)) > 0L, logical(1))]
  } else {
    available <- character(0)
  }

  cur_type <- slot(x, "Type")
  max_dim  <- if (!is.na(cur_type) && is(ae, "SingleCellExperiment") &&
                  cur_type %in% reducedDimNames(ae)) {
    ncol(reducedDim(ae, cur_type))
  } else 2L

  .addSpecificTour(class(x)[1], "Experiment", function(plot_name) {
    data.frame(rbind(
      c(element = paste0("#", plot_name, "_Experiment + .selectize-control"),
        intro = "Select which alternative experiment to visualise.
Only alternative experiments that contain at least one non-empty
<code>reducedDim</code> result are useful here."),
      c(element = paste0("#", plot_name, "_Type + .selectize-control"),
        intro = "Choose the dimensionality reduction result from the
selected alternative experiment — for example <code>PCA</code> or
<code>UMAP</code>.  Only results with at least one dimension are listed."),
      c(element = paste0("#", plot_name, "_XAxis + .selectize-control"),
        intro = "The component of the chosen result to place on the x-axis."),
      c(element = paste0("#", plot_name, "_YAxis + .selectize-control"),
        intro = "The component of the chosen result to place on the y-axis.")
    ), stringsAsFactors = FALSE)
  })

  list(
    selectInput(.input_FUN("Experiment"),
      label    = "Experiment:",
      choices  = all_altexps,
      selected = iSEE:::.choose_link(current_ae, all_altexps)
    ),
    .selectInput.iSEE(x, "Type",
      label    = "Type:",
      choices  = available,
      selected = cur_type
    ),
    .selectInput.iSEE(x, "XAxis",
      label    = "Dimension 1:",
      choices  = seq_len(max_dim),
      selected = slot(x, "XAxis")
    ),
    .selectInput.iSEE(x, "YAxis",
      label    = "Dimension 2:",
      choices  = seq_len(max_dim),
      selected = slot(x, "YAxis")
    )
  )
})

#' @importFrom iSEE .getEncodedName .createProtectedParameterObservers
#' @importFrom SingleCellExperiment altExp reducedDim reducedDimNames
#' @importFrom shiny observeEvent updateSelectInput
setMethod(".createObservers", "AltReducedDimensionPlot",
    function(x, se, input, session, pObjects, rObjects) {
  callNextMethod()

  plot_name  <- .getEncodedName(x)
  ae_field   <- paste0(plot_name, "_Experiment")
  type_field <- paste0(plot_name, "_Type")
  x_field    <- paste0(plot_name, "_XAxis")
  y_field    <- paste0(plot_name, "_YAxis")

  .createProtectedParameterObservers(plot_name,
    fields   = c("XAxis", "YAxis"),
    input    = input,
    pObjects = pObjects,
    rObjects = rObjects)

  # nocov start
  # Type: update XAxis/YAxis max when the reduction type changes
  observeEvent(input[[type_field]], {
    matched <- input[[type_field]]
    if (is.null(matched) || !nzchar(matched)) return(NULL)
    if (identical(matched, pObjects$memory[[plot_name]][["Type"]])) return(NULL)
    pObjects$memory[[plot_name]][["Type"]] <- matched

    cur_ae <- pObjects$memory[[plot_name]][["Experiment"]]
    ae <- if (identical(cur_ae, .selectionMainExpTitle)) se else altExp(se, cur_ae)
    new_max <- ncol(reducedDim(ae, matched))
    cap_X   <- pmin(new_max, pObjects$memory[[plot_name]][["XAxis"]])
    cap_Y   <- pmin(new_max, pObjects$memory[[plot_name]][["YAxis"]])
    pObjects$memory[[plot_name]][["XAxis"]] <- cap_X
    pObjects$memory[[plot_name]][["YAxis"]] <- cap_Y

    updateSelectInput(session, x_field, choices = seq_len(new_max), selected = cap_X)
    updateSelectInput(session, y_field, choices = seq_len(new_max), selected = cap_Y)
    iSEE:::.requestCleanUpdate(plot_name, pObjects, rObjects)
  }, ignoreInit = TRUE)

  # AltExp: update Type choices and reset axes when the experiment changes
  observeEvent(input[[ae_field]], {
    matched_ae <- input[[ae_field]]
    if (identical(matched_ae, pObjects$memory[[plot_name]][["Experiment"]])) return(NULL)
    pObjects$memory[[plot_name]][["Experiment"]] <- matched_ae

    ae <- if (identical(matched_ae, .selectionMainExpTitle)) se else altExp(se, matched_ae)
    if (!is(ae, "SingleCellExperiment")) {
      pObjects$memory[[plot_name]][["Type"]] <- NA_character_
      iSEE:::.requestCleanUpdate(plot_name, pObjects, rObjects)
      return(NULL)
    }
    valid <- reducedDimNames(ae)
    valid <- valid[vapply(valid, function(n) ncol(reducedDim(ae, n)) > 0L, logical(1))]

    if (length(valid) == 0L) {
      pObjects$memory[[plot_name]][["Type"]] <- NA_character_
      iSEE:::.requestCleanUpdate(plot_name, pObjects, rObjects)
      return(NULL)
    }

    cur_type <- pObjects$memory[[plot_name]][["Type"]]
    new_type <- if (!is.na(cur_type) && cur_type %in% valid) cur_type else valid[1]
    pObjects$memory[[plot_name]][["Type"]] <- new_type

    new_max <- ncol(reducedDim(ae, new_type))
    pObjects$memory[[plot_name]][["XAxis"]] <- 1L
    pObjects$memory[[plot_name]][["YAxis"]] <- min(new_max, 2L)

    updateSelectInput(session, type_field, choices = valid,          selected = new_type)
    updateSelectInput(session, x_field,    choices = seq_len(new_max), selected = 1L)
    updateSelectInput(session, y_field,    choices = seq_len(new_max), selected = min(new_max, 2L))
    iSEE:::.requestCleanUpdate(plot_name, pObjects, rObjects)
  }, ignoreInit = TRUE)
  # nocov end

  invisible(NULL)
})

#' @importFrom iSEE .textEval
#' @importFrom SingleCellExperiment altExp reducedDim
setMethod(".generateDotPlotData", "AltReducedDimensionPlot", function(x, envir) {
  ae_name <- slot(x, "Experiment")
  type    <- slot(x, "Type")
  xaxis   <- slot(x, "XAxis")
  yaxis   <- slot(x, "YAxis")

  # Guard: altExp is not an SCE or has no reducedDims — Type is NA_character_.
  # Produce NA-filled plot.data with ncol(se) rows so downstream colour/shape/
  # size methods can still attach colData columns without length mismatches.
  if (is.na(type)) {
    data_cmds <- paste0(
      "plot.data <- data.frame(",
      "X = rep(NA_real_, ncol(se)), Y = rep(NA_real_, ncol(se)),",
      " row.names = colnames(se));"
    )
    .textEval(data_cmds, envir)
    return(list(
      commands = data_cmds,
      labels   = list(
        title = paste0("Select another altExp that has a reducedDims slot",
                       " with dimension reduction results"),
        X = "", Y = ""
      )
    ))
  }

  data_cmds <- c(
    if (identical(ae_name, .selectionMainExpTitle)) "ae <- se;" else
      sprintf("ae <- altExp(se, %s);", deparse(ae_name)),
    sprintf("red.dim <- reducedDim(ae, %s);", deparse(type)),
    sprintf(
      "plot.data <- data.frame(X = red.dim[, %i], Y = red.dim[, %i], row.names = colnames(se));",
      xaxis, yaxis)
  )

  .textEval(data_cmds, envir)

  list(
    commands = data_cmds,
    labels   = list(
      title = sprintf("%s (%s)", ae_name, type),
      X     = sprintf("Dimension %s", xaxis),
      Y     = sprintf("Dimension %s", yaxis)
    )
  )
})

#' @importFrom iSEE .generateDotPlot .textEval
setMethod(".generateDotPlot", "AltReducedDimensionPlot", function(x, labels, envir) {
  if (is.na(slot(x, "Type"))) {
    # Use plot.data with aes(x = X, y = Y) so that coordinfo$mapping$x is "X".
    # Without this mapping, Shiny's nearPoints() cannot infer xvar and throws an
    # error whenever the user clicks or hovers over the panel.  Explicit limits
    # ensure the coordinate system is well-defined even though all X/Y are NA.
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
setMethod(".definePanelTour", "AltReducedDimensionPlot", function(x) {
  collated <- rbind(
    c(
      element = paste0("#", .getEncodedName(x)),
      intro   = sprintf(
        "The <font color=\"%s\">Alt reduced dimension plot</font> panel
displays dimensionality reduction results (PCA, UMAP, etc.) stored in a
<code>reducedDim</code> of an <em>alternative experiment</em> inside a
<code>SingleCellExperiment</code>.  Each point is a cell; coordinates come
from the selected altExp rather than the main experiment.",
        .getPanelColor(x))
    ),
    .addTourStep(x, .dataParamBoxOpen,
      "The <i>Data parameters</i> box lets you choose the alternative
experiment, the dimensionality reduction type, and which two dimensions
to display.<br><br><strong>Action:</strong> click to expand.")
  )

  parent_tour <- callNextMethod()

  rbind(
    data.frame(element = collated[, 1], intro = collated[, 2],
               stringsAsFactors = FALSE),
    parent_tour
  )
})
