############################################################
# LinkedFeaturesAssayPlot
############################################################

# Extends FeatureAssayPlot from iSEE to support linked feature visualisation
# across the main experiment and altExps. Decouples feature *selection*
# (from SelectionExperiment) from feature *visualisation* (from Experiment); either
# side can independently be the main experiment or any altExp.  LookupColumn
# lives in SelectionExperiment's rowData; MapColumn lives in AltExp's rowData.
#
# Works with stock iSEE: .addDotPlotDataColor/.addDotPlotDataShape/
# .addDotPlotDataSize/.addDotPlotDataFacets are internal (non-exported) S4
# generics already defined on ColumnDotPlot in iSEE itself. setMethod() needs
# the generic resolvable by name in this namespace, which normally requires
# an export; instead we bind these names locally to iSEE's internal generics
# (below) and keep the usual setMethod("<name>", ...) string form -- passing
# `iSEE:::.addDotPlotDataColor` directly as setMethod()'s first argument
# breaks roxygen2::document(), which evaluates that argument while parsing
# and expects a plain string. S4 dispatch on a generic is global to the
# generic regardless of export status, so methods attached this way are
# still found when iSEE's own internal code calls these generics on a
# LinkedFeaturesAssayPlot instance.
.addDotPlotDataColor  <- iSEE:::.addDotPlotDataColor
.addDotPlotDataShape  <- iSEE:::.addDotPlotDataShape
.addDotPlotDataSize   <- iSEE:::.addDotPlotDataSize
.addDotPlotDataFacets <- iSEE:::.addDotPlotDataFacets

#' LinkedFeaturesAssayPlot: linked feature-selection and visualisation
#'
#' An S4 class extending \code{\linkS4class{FeatureAssayPlot}} that decouples
#' \emph{feature selection} from \emph{feature visualisation} across the main
#' experiment and the alternative experiments stored in a
#' \code{\linkS4class{SingleCellExperiment}}.
#'
#' A feature is selected from \code{SelectionExperiment} — an altExp, or the main
#' experiment via the \code{"(Main)"} sentinel — typically via a
#' \code{\linkS4class{AltRowDataTable}} (or any panel transmitting a row
#' selection) targeting that source. The value in \code{LookupColumn} for that
#' feature is used as a join key to retrieve matching rows in \code{Experiment} —
#' itself an altExp, or the main experiment via the same sentinel — via
#' \code{MapColumn}.  All matching rows are displayed as separate traces per
#' sample.
#'
#' @slot SelectionExperiment \code{character(1)}.  Name of the alternative
#'   experiment from which the feature is selected and in whose
#'   \code{rowData} \code{LookupColumn} is looked up.  May also be the
#'   sentinel \code{"(Main)"} to select features from — and look
#'   up \code{LookupColumn} in — the main experiment instead of an altExp.
#' @slot LookupColumn \code{character(1)}.  Name of a \code{character} or
#'   \code{factor} column in \code{rowData} of \code{SelectionExperiment} (or the
#'   sentinel \code{"(rowname)"} to use the rowname itself as the join key).
#' @slot Experiment \code{character(1)}.  Name of the alternative experiment whose
#'   feature values are visualised.  May also be the sentinel
#'   \code{"(Main)"} to visualise main experiment assay values
#'   instead of an altExp's.
#' @slot AltAssay \code{character(1)}.  Assay within \code{Experiment} to plot.
#' @slot MapColumn \code{character(1)}.  Name of a \code{character} or
#'   \code{factor} column in \code{rowData} of \code{Experiment} (or
#'   \code{"(rowname)"}).  Rows whose value matches the lookup key are shown.
#' @slot PlotType \code{character(1)}.  One of \code{"Auto"},
#'   \code{"Scatter"}, or \code{"Scatter + lines"}.  Defaults to
#'   \code{"Scatter"} because the long-format \code{plot.data} produced by
#'   this panel is not compatible with iSEE's violin/square renderers.
#'
#' @section Inherited slots:
#' All slots from \code{\linkS4class{FeatureAssayPlot}} are inherited,
#' including those controlling x-axis choice, colour, shape, size, faceting,
#' and multi-selection behaviour.
#'
#' @seealso
#' \code{\link{LinkedFeaturesAssayPlot}} for the constructor.
#' \code{\linkS4class{AltRowDataTable}} for a compatible row-selection source.
#'
#' @name LinkedFeaturesAssayPlot-class
#' @rdname LinkedFeaturesAssayPlot-class
#' @exportClass LinkedFeaturesAssayPlot
NULL

# Sentinel value for the `SelectionExperiment` slot/UI indicating that the
# feature should be selected from the main experiment's rownames/rowData
# instead of from a named altExp.
.selectionMainExpTitle <- "(Main)"

# Returns assay names for `x`, falling back to integer-as-character indices
# ("1", "2", …) when assays exist but carry no names.
#' @importFrom SummarizedExperiment assayNames assays
.safe_assay_names <- function(x) {
  nms <- assayNames(x)
  if (length(nms) > 0L) nms else as.character(seq_len(length(assays(x, withDimnames = FALSE))))
}

setClass("LinkedFeaturesAssayPlot",
  contains = "FeatureAssayPlot",
  slots    = c(
    SelectionExperiment = "character",
    Experiment      = "character",
    AltAssay        = "character",
    LookupColumn    = "character",
    MapColumn       = "character",
    PlotType        = "character"
  ))

#' @importFrom SingleCellExperiment altExp
.altExpPlotLinked_selectionRownames <- function(se, sel_name) {
  if (identical(sel_name, .selectionMainExpTitle)) {
    rownames(se)
  } else {
    rownames(altExp(se, sel_name))
  }
}

#' @importFrom iSEE .validStringError
setValidity2("LinkedFeaturesAssayPlot", function(object) {
  msg <- character(0)
  msg <- .validStringError(msg, object, "PlotType")
  if (length(msg)) return(msg)
  TRUE
})

#' @importFrom iSEE .emptyDefault
setMethod("initialize", "LinkedFeaturesAssayPlot", function(.Object, ...) {
  args <- list(...)
  args <- .emptyDefault(args, "SelectionExperiment", NA_character_)
  args <- .emptyDefault(args, "Experiment",          NA_character_)
  args <- .emptyDefault(args, "AltAssay",        NA_character_)
  args <- .emptyDefault(args, "LookupColumn",    NA_character_)
  args <- .emptyDefault(args, "MapColumn",       NA_character_)
  args <- .emptyDefault(args, "PlotType",        "Scatter")
  args <- .emptyDefault(args, iSEE:::.shapeByField,   iSEE:::.shapeByColDataTitle)
  args <- .emptyDefault(args, iSEE:::.shapeByColData, "altExp_feature_id")
  do.call(callNextMethod, c(list(.Object), args))
})

#' Construct an LinkedFeaturesAssayPlot panel
#'
#' Creates an instance of \code{\linkS4class{LinkedFeaturesAssayPlot}} for use as a
#' panel in an iSEE application.
#'
#' @param ... Named arguments corresponding to slots of
#'   \code{\linkS4class{LinkedFeaturesAssayPlot}} or its parent classes.
#'
#' @return An \code{\linkS4class{LinkedFeaturesAssayPlot}} object.
#'
#' @examples
#' \dontrun{
#' # Minimal panel using runtime defaults
#' LinkedFeaturesAssayPlot()
#'
#'
#' # Select a protein from the main experiment and visualise its peptides
#' LinkedFeaturesAssayPlot(
#'   SelectionExperiment    = "(Main)",
#'   LookupColumn       = "Proteins",
#'   Experiment            = "peptides",
#'   AltAssay           = "peptides_norm",
#'   MapColumn          = "Proteins"
#' )
#'
#' # Select a peptide from an altExp and visualise its protein in the main
#' # experiment
#' LinkedFeaturesAssayPlot(
#'   SelectionExperiment    = "peptides",
#'   LookupColumn       = "Proteins",
#'   Experiment            = "(Main)",
#'   AltAssay           = "proteins",
#'   MapColumn          = "Proteins"
#' )
#' 
#' #######################
#' # Full example with app
#' #######################
#' 
#' library(iSEEu)
#' # Define ReducedDimensionPlot, VolcanoPlot, RowDataTable
#' # LinkedFeatureAssayPlot Panels for the main assay
#' rdp <- ReducedDimensionPlot()
#' vp <- VolcanoPlot()
#' rdt <- RowDataTable(RowSelectionSource = "VolcanoPlot1")
#' mlfap <- LinkedFeaturesAssayPlot(
#' SelectionExperiment = "(Main)",
#' LookupColumn = "(rowname)",
#' Experiment= "(Main)",
#' AltAssay = "proteins",
#' MapColumn = "(rowname)",
#' YAxisFeatureSource = "RowDataTable1",
#' PlotType = "Scatter + lines",
#' XAxis = "Column data",
#' XAxisColumnData = "sampleId"
#' )
#' 
#' # Define LinkedFeatureAssayPlot for peptides altExp assay  
#' alfap <- LinkedFeaturesAssayPlot(
#' YAxisFeatureSource = "RowDataTable1",
#' SelectionExperiment = "(Main)",
#' LookupColumn = "Proteins",
#' Experiment= "peptides",
#' AltAssay = "peptides_norm",
#' PlotType = "Scatter + lines",
#' XAxis = "Column data",
#' XAxisColumnData = "sampleId",
#' MapColumn = "Proteins"
#' )
#' 
#' # Run app 
#' data(sceProteinsPeptides)
#' app <- iSEE(
#' sceProteinsPeptides,
#' initial=list(rdp, vp, rdt, mlfap, alfap)
#' )
#' 
#' if (interactive()) {
#' shiny::runApp(app, launch.browser = TRUE)
#' }
#' ### End Full example 
#' }
#'
#' @seealso \code{\linkS4class{LinkedFeaturesAssayPlot}} for slot details.
#' @export
LinkedFeaturesAssayPlot <- function(...) new("LinkedFeaturesAssayPlot", ...)

setMethod(".fullName",   "LinkedFeaturesAssayPlot", function(x) "Linked features assay plot")
setMethod(".panelColor", "LinkedFeaturesAssayPlot", function(x) "#5500AA")

############################################################
# Cache
############################################################

#' @importFrom iSEE .getCachedCommonInfo .setCachedCommonInfo
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment rowData
setMethod(".cacheCommonInfo", "LinkedFeaturesAssayPlot", function(x, se) {
  if (!is.null(.getCachedCommonInfo(se, "LinkedFeaturesAssayPlot"))) return(se)
  se <- callNextMethod()

  valid_lookcols_by_ae <- lapply(altExpNames(se), function(ae_name) {
    rd_cls <- sapply(rowData(altExp(se, ae_name)), class)
    c("(rowname)", names(rd_cls)[rd_cls %in% c("character", "factor")])
  })
  names(valid_lookcols_by_ae) <- altExpNames(se)

  # The main experiment is also a valid selection source AND visualisation
  # target, via the .selectionMainExpTitle sentinel.
  main_rd_cls <- sapply(rowData(se), class)
  main_valid_cols <- c("(rowname)",
    names(main_rd_cls)[main_rd_cls %in% c("character", "factor")])
  valid_lookcols_by_ae[[.selectionMainExpTitle]] <- main_valid_cols

  valid_mapcols_by_ae <- lapply(altExpNames(se), function(ae_name) {
    rd_cls <- sapply(rowData(altExp(se, ae_name)), class)
    c("(rowname)", names(rd_cls)[rd_cls %in% c("character", "factor")])
  })
  names(valid_mapcols_by_ae) <- altExpNames(se)
  valid_mapcols_by_ae[[.selectionMainExpTitle]] <- main_valid_cols

  .setCachedCommonInfo(se, "LinkedFeaturesAssayPlot",
    valid.altExp.names       = altExpNames(se),
    valid.selection.names    = c(altExpNames(se), .selectionMainExpTitle),
    valid.lookcols.by.altExp = valid_lookcols_by_ae,
    valid.mapcols.by.altExp  = valid_mapcols_by_ae)
})

############################################################
# Parameter refinement
############################################################

#' @importFrom iSEE .replaceMissingWithFirst .getCachedCommonInfo .findAtomicFields
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames colData
setMethod(".refineParameters", "LinkedFeaturesAssayPlot", function(x, se) {
  # Save slots that callNextMethod() (-> FeatureAssayPlot -> ColumnDotPlot) may
  # clobber with its own defaults before we get a chance to act.
  saved_y         <- slot(x, iSEE:::.featAssayYAxisFeatName)
  saved_color_col <- slot(x, iSEE:::.colorByColData)
  saved_shape_col <- slot(x, iSEE:::.shapeByColData)

  x <- callNextMethod()
  if (is.null(x)) return(NULL)

  # SelectionExperiment and its feature name (the main experiment is also a
  # valid selection source, via the .selectionMainExpTitle sentinel)
  x <- .replaceMissingWithFirst(x, "SelectionExperiment",
                                 c(altExpNames(se), .selectionMainExpTitle))
  sel_ae_name <- slot(x, "SelectionExperiment")
  sel_rn      <- .altExpPlotLinked_selectionRownames(se, sel_ae_name)

  slot(x, iSEE:::.featAssayYAxisFeatName) <-
    if (!is.na(saved_y) && nzchar(saved_y) && saved_y %in% sel_rn) {
      saved_y
    } else {
      sel_rn[1]
    }

  # Validate LookupColumn against SelectionExperiment rowData
  cache    <- .getCachedCommonInfo(se, "LinkedFeaturesAssayPlot")
  valid_lc <- cache$valid.lookcols.by.altExp[[sel_ae_name]]
  cur_lc   <- slot(x, "LookupColumn")
  if (is.na(cur_lc) || !nzchar(cur_lc) || !cur_lc %in% valid_lc) {
    slot(x, "LookupColumn") <- valid_lc[1]
  }

  # Validate Experiment and AltAssay (the main experiment is also a valid
  # visualisation target, via the .selectionMainExpTitle sentinel)
  x <- .replaceMissingWithFirst(x, "Experiment", c(altExpNames(se), .selectionMainExpTitle))
  vis_ae_name <- slot(x, "Experiment")

  all_assays <- if (identical(vis_ae_name, .selectionMainExpTitle)) {
    .safe_assay_names(se)
  } else {
    .safe_assay_names(altExp(se, vis_ae_name))
  }
  if (length(all_assays) == 0L) {
    warning(sprintf("no valid 'assays' for plotting '%s'", class(x)[1]))
    return(NULL)
  }
  if (!slot(x, "AltAssay") %in% all_assays) {
    x <- .replaceMissingWithFirst(x, "AltAssay", all_assays)
  }

  # Validate MapColumn against visualisation Experiment rowData
  valid_mc <- cache$valid.mapcols.by.altExp[[vis_ae_name]]
  cur_mc   <- slot(x, "MapColumn")
  if (is.na(cur_mc) || !nzchar(cur_mc) || !cur_mc %in% valid_mc) {
    slot(x, "MapColumn") <- valid_mc[1]
  }

  # Shape/colour: "altExp_feature_id" is not a real colData column, so
  # ColumnDotPlot.refineParameters (called above via callNextMethod) replaces
  # it with the first real colData column.  We use the value saved BEFORE
  # callNextMethod to distinguish two cases:
  #
  #   saved == NA / "" / "altExp_feature_id"  --> startup or still at our
  #     default; restore "altExp_feature_id" so features are coloured/shaped
  #     by their altExp identity by default.
  #
  #   saved == a specific colData column the user chose --> keep it (and only
  #     reset if that column has since disappeared from the data).
  cd_cls <- sapply(colData(se), class)
  valid_shape_cols <- c("altExp_feature_id",
                        names(cd_cls)[cd_cls %in% c("factor", "character")])
  if (is.na(saved_shape_col) || !nzchar(saved_shape_col) ||
      identical(saved_shape_col, "altExp_feature_id") ||
      !saved_shape_col %in% valid_shape_cols) {
    slot(x, iSEE:::.shapeByField)   <- iSEE:::.shapeByColDataTitle
    slot(x, iSEE:::.shapeByColData) <- "altExp_feature_id"
  } else {
    slot(x, iSEE:::.shapeByField)   <- iSEE:::.shapeByColDataTitle
    slot(x, iSEE:::.shapeByColData) <- saved_shape_col
  }

  valid_color_cols <- c("altExp_feature_id", .findAtomicFields(colData(se)))
  if (is.na(saved_color_col) || !nzchar(saved_color_col) ||
      identical(saved_color_col, "altExp_feature_id") ||
      !saved_color_col %in% valid_color_cols) {
    slot(x, iSEE:::.colorByField)   <- iSEE:::.colorByColDataTitle
    slot(x, iSEE:::.colorByColData) <- "altExp_feature_id"
  } else {
    slot(x, iSEE:::.colorByField)   <- iSEE:::.colorByColDataTitle
    slot(x, iSEE:::.colorByColData) <- saved_color_col
  }

  x
})

############################################################
# Data interface
############################################################

#' @importFrom iSEE .getEncodedName .getCachedCommonInfo .selectizeInput.iSEE
#'   .radioButtons.iSEE .conditionalOnRadio .addSpecificTour .findAtomicFields
#' @importFrom SingleCellExperiment altExpNames altExp
#' @importFrom SummarizedExperiment assayNames rowData colData
#' @importFrom shiny selectInput checkboxInput
setMethod(".defineDataInterface", "LinkedFeaturesAssayPlot", function(x, se, select_info) {
  panel_name <- .getEncodedName(x)
  .input_FUN <- function(field) paste0(panel_name, "_", field)

  all_altexps     <- altExpNames(se)
  all_sel_sources <- c(all_altexps, .selectionMainExpTitle)
  all_vis_sources <- c(all_altexps, .selectionMainExpTitle)
  current_sel_ae  <- slot(x, "SelectionExperiment")
  current_ae      <- slot(x, "Experiment")

  all_assays <- if (identical(current_ae, .selectionMainExpTitle)) {
    .safe_assay_names(se)
  } else if (!is.na(current_ae) && current_ae %in% all_altexps) {
    .safe_assay_names(altExp(se, current_ae))
  } else {
    c(.safe_assay_names(se), unlist(lapply(all_altexps, function(nm) .safe_assay_names(altExp(se, nm)))))
  }

  cache    <- .getCachedCommonInfo(se, "LinkedFeaturesAssayPlot")
  valid_lc <- cache$valid.lookcols.by.altExp[[current_sel_ae]]
  valid_mc <- cache$valid.mapcols.by.altExp[[current_ae]]

  tab_by_row <- select_info$single$feature
  # Read colData columns directly from se so X-axis choices are always
  # populated regardless of ColumnDotPlot cache state.
  column_covariates <- .findAtomicFields(colData(se))

  xaxis_choices <- iSEE:::.featAssayXAxisNothingTitle
  if (length(column_covariates)) {
    xaxis_choices <- c(xaxis_choices, iSEE:::.featAssayXAxisColDataTitle)
  }

  .addSpecificTour(class(x)[1], iSEE:::.featAssayYAxisFeatName, function(plot_name) {
    data.frame(rbind(
      c(element = paste0("#", plot_name, "_", iSEE:::.featAssayYAxisFeatName,
                         " + .selectize-control"),
        intro = "Choose the feature from the <em>selection source</em> (an
altExp, or the main experiment) whose row annotation will be used as the
lookup key to retrieve matching rows in the visualisation source."),
      c(element = paste0("#", plot_name, "_", iSEE:::.featAssayYAxisRowTable,
                         " + .selectize-control"),
        intro = "The feature selection can be driven automatically by a row
selection made in another panel, for example an
<em>AltExp row data table</em>.")
    ), stringsAsFactors = FALSE)
  })

  .addSpecificTour(class(x)[1], "SelectionExperiment", function(plot_name) {
    data.frame(rbind(
      c(element = paste0("#", plot_name, "_SelectionExperiment + .selectize-control"),
        intro = "Choose where a feature is selected from: an alternative
experiment, or the main experiment itself (<em>Main experiment</em>).  The
<em>Lookup column</em> must be a column in this source's <code>rowData</code>."),
      c(element = paste0("#", plot_name, "_LookupColumn + .selectize-control"),
        intro = "Choose a <code>character</code> or <code>factor</code> column
from the <em>selection source</em>'s <code>rowData</code>.  The value stored
for the selected feature in this column is the join key used to retrieve rows
from the visualisation source.")
    ), stringsAsFactors = FALSE)
  })

  .addSpecificTour(class(x)[1], "Experiment", function(plot_name) {
    data.frame(rbind(
      c(element = paste0("#", plot_name, "_Experiment + .selectize-control"),
        intro = "Choose where feature values are displayed from: an
alternative experiment, or the main experiment itself (<em>Main
experiment</em>).  This can be the same as, or different from, the selection
source."),
      c(element = paste0("#", plot_name, "_AltAssay + .selectize-control"),
        intro = "Choose which assay within the visualisation source to display
on the y-axis."),
      c(element = paste0("#", plot_name, "_MapColumn + .selectize-control"),
        intro = "Choose a column from the <em>visualisation source</em>'s
<code>rowData</code>.  Rows whose value matches the key from the Lookup column
are displayed."),
      c(element = paste0("#", plot_name, "_PlotType"),
        intro = "<strong>Scatter</strong> produces a dot plot.
<strong>Scatter + lines</strong> additionally connects dots for each altExp
feature with a line — useful for ordered x-axes such as sample IDs.")
    ), stringsAsFactors = FALSE)
  })

  list(
    # ---- Feature selection (from SelectionExperiment) -----------------------
    .selectizeInput.iSEE(x, iSEE:::.featAssayYAxisFeatName,
      label    = "Feature in selection experiment:",
      choices  = NULL, selected = NULL, multiple = FALSE),
    selectInput(.input_FUN(iSEE:::.featAssayYAxisRowTable),
      label    = NULL,
      choices  = tab_by_row,
      selected = iSEE:::.choose_link(slot(x, iSEE:::.featAssayYAxisRowTable),
                                     tab_by_row)),
    checkboxInput(.input_FUN(iSEE:::.featAssayYAxisFeatDynamic),
      label = "Use dynamic feature selection for the y-axis",
      value = slot(x, iSEE:::.featAssayYAxisFeatDynamic)),
    # ---- Selection source (an altExp, or the main experiment) -----------
    selectInput(.input_FUN("SelectionExperiment"),
      label    = "Selection source:",
      choices  = all_sel_sources,
      selected = iSEE:::.choose_link(current_sel_ae, all_sel_sources)),
    selectInput(.input_FUN("LookupColumn"),
      label    = "Lookup column (selection source):",
      choices  = valid_lc,
      selected = iSEE:::.choose_link(slot(x, "LookupColumn"), valid_lc)),
    # ---- Visualisation source (an altExp, or the main experiment) -------
    selectInput(.input_FUN("Experiment"),
      label    = "Visualisation source:",
      choices  = all_vis_sources,
      selected = iSEE:::.choose_link(current_ae, all_vis_sources)),
    selectInput(.input_FUN("AltAssay"),
      label    = "Visualisation assay:",
      choices  = all_assays,
      selected = iSEE:::.choose_link(slot(x, "AltAssay"), all_assays)),
    selectInput(.input_FUN("MapColumn"),
      label    = "Map column (visualisation source):",
      choices  = valid_mc,
      selected = iSEE:::.choose_link(slot(x, "MapColumn"), valid_mc)),
    selectInput(.input_FUN("PlotType"),
      label    = "Plot type:",
      choices  = c("Auto", "Scatter", "Scatter + lines"),
      selected = slot(x, "PlotType")),
    # ---- X-axis ---------------------------------------------------------
    .radioButtons.iSEE(x, iSEE:::.featAssayXAxis,
      label    = "X-axis:",
      inline   = TRUE,
      choices  = xaxis_choices,
      selected = slot(x, iSEE:::.featAssayXAxis)),
    .conditionalOnRadio(.input_FUN(iSEE:::.featAssayXAxis),
      iSEE:::.featAssayXAxisColDataTitle,
      selectInput(.input_FUN(iSEE:::.featAssayXAxisColData),
        label    = "X-axis column data:",
        choices  = column_covariates,
        selected = slot(x, iSEE:::.featAssayXAxisColData)))
  )
})

############################################################
# Observers
############################################################

#' @importFrom iSEE .getEncodedName .createUnprotectedParameterObservers
#'   .createProtectedParameterObservers .trackSingleSelection .trackRelinkedSelection
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment assayNames
#' @importFrom shiny isolate observe observeEvent updateSelectInput
#'   updateSelectizeInput
setMethod(".createObservers", "LinkedFeaturesAssayPlot",
    function(x, se, input, session, pObjects, rObjects) {
  callNextMethod()

  plot_name    <- .getEncodedName(x)
  ae_field     <- paste0(plot_name, "_Experiment")
  assay_field  <- paste0(plot_name, "_AltAssay")
  sel_ae_field <- paste0(plot_name, "_SelectionExperiment")
  feat_field   <- paste0(plot_name, "_", iSEE:::.featAssayYAxisFeatName)
  source_field <- paste0(plot_name, "_", iSEE:::.featAssayYAxisRowTable)

  # Unprotected observers for custom slots that need no special side-effects.
  # Experiment is handled by its own observeEvent below so the AltAssay dropdown
  # is updated atomically — keeping both in an unprotected observer causes the
  # same race-condition freeze seen in AltFeatureAssayPlot.
  .createUnprotectedParameterObservers(plot_name,
    fields   = c("AltAssay", "LookupColumn", "MapColumn", "PlotType"),
    input    = input,
    pObjects = pObjects,
    rObjects = rObjects)

  # Stock FeatureAssayPlot relies on .create_dimname_observers (tied to the
  # X-axis *feature name* single-selection slot) to also detect changes to the
  # XAxis mode radio itself; that mechanism reads input[[..._XAxisRowTable]]
  # and bails out silently when that input doesn't exist. Our
  # .defineDataInterface never renders an X-axis feature-name UI (only "None"
  # and "Column data" are offered), so that input never exists and the radio
  # change is never picked up. Wire it explicitly here instead.
  .createProtectedParameterObservers(plot_name,
    fields   = iSEE:::.featAssayXAxis,
    input    = input,
    pObjects = pObjects,
    rObjects = rObjects)

  # Start Experiment change: reset AltAssay to first valid assay of the new
  # visualisation source (an altExp, or the main experiment) so that
  # .generateDotPlotData never receives a stale assay name.
  observeEvent(input[[ae_field]], {
    new_ae <- input[[ae_field]]
    if (identical(new_ae, pObjects$memory[[plot_name]][["Experiment"]])) return(NULL)
    pObjects$memory[[plot_name]][["Experiment"]] <- new_ae

    new_assays <- if (identical(new_ae, .selectionMainExpTitle)) {
      .safe_assay_names(se)
    } else {
      .safe_assay_names(altExp(se, new_ae))
    }
    cur_assay  <- pObjects$memory[[plot_name]][["AltAssay"]]
    new_assay  <- if (!is.na(cur_assay) && nzchar(cur_assay) &&
                       cur_assay %in% new_assays) cur_assay else new_assays[1]
    pObjects$memory[[plot_name]][["AltAssay"]] <- new_assay
    updateSelectInput(session, assay_field, choices = new_assays, selected = new_assay)
    iSEE:::.requestCleanUpdate(plot_name, pObjects, rObjects)
  }, ignoreInit = TRUE)
  # end

  # Repopulate the feature selectize with SelectionExperiment rownames.
  # Scheduled via session$onFlushed so it fires AFTER FeatureAssayPlot's own
  # updateSelectizeInput(choices=rownames(se)), overwriting it in the browser.
  repopulate_feat_choices <- function() {
    current_sel_ae <- pObjects$memory[[plot_name]][["SelectionExperiment"]]
    sel_ae_rn <- .altExpPlotLinked_selectionRownames(se, current_sel_ae)
    updateSelectizeInput(session, feat_field,
      choices  = sel_ae_rn,
      selected = pObjects$memory[[plot_name]][[iSEE:::.featAssayYAxisFeatName]],
      server   = TRUE)
  }

  session$onFlushed(repopulate_feat_choices, once = TRUE)

  # Bypass: intercept single-selection from the source panel and write directly
  # into pObjects$memory, bypassing the choices=rownames(se) gate in iSEE's
  # inherited observer that silently drops altExp feature names.
  observe({
    source_panel <- input[[source_field]]
    if (!length(source_panel) || !nzchar(source_panel) ||
        identical(source_panel, iSEE:::.noSelection)) return()

    iSEE:::.trackSingleSelection(source_panel, rObjects)

    isolate({
      new_feat <- iSEE:::.singleSelectionValue(
        pObjects$memory[[source_panel]],
        pObjects$contents[[source_panel]])
      if (is.null(new_feat) || !nzchar(new_feat)) return()

      current_sel_ae <- pObjects$memory[[plot_name]][["SelectionExperiment"]]
      sel_ae_rn <- .altExpPlotLinked_selectionRownames(se, current_sel_ae)
      if (!new_feat %in% sel_ae_rn) return()

      old_feat <- pObjects$memory[[plot_name]][[iSEE:::.featAssayYAxisFeatName]]
      if (identical(new_feat, old_feat)) return()

      pObjects$memory[[plot_name]][[iSEE:::.featAssayYAxisFeatName]] <- new_feat
      updateSelectizeInput(session, feat_field,
        choices = sel_ae_rn, selected = new_feat, server = TRUE)
      iSEE:::.requestCleanUpdate(plot_name, pObjects, rObjects)
    })
  })

  # Re-populate after source panel relinking so the selectize always shows
  # SelectionExperiment choices, not rownames(se).
  observe({
    source_panel <- input[[source_field]]
    if (length(source_panel) && nzchar(source_panel)) {
      .trackSingleSelection(source_panel, rObjects)
    }
    .trackRelinkedSelection(plot_name, rObjects)
    session$onFlushed(repopulate_feat_choices, once = TRUE)
  }, priority = -1L)

  # start SelectionExperiment change: update feature selectize + trigger full re-render
  # so that LookupColumn choices in defineDataInterface reflect the new altExp.
  #  
  observeEvent(input[[sel_ae_field]], {
    new_sel_ae <- input[[sel_ae_field]]
    if (identical(new_sel_ae, pObjects$memory[[plot_name]][["SelectionExperiment"]])) {
      return(NULL)
    }
    pObjects$memory[[plot_name]][["SelectionExperiment"]] <- new_sel_ae

    sel_ae_rn <- .altExpPlotLinked_selectionRownames(se, new_sel_ae)
    cur_feat  <- pObjects$memory[[plot_name]][[iSEE:::.featAssayYAxisFeatName]]
    new_feat  <- if (!is.na(cur_feat) && nzchar(cur_feat) && cur_feat %in% sel_ae_rn) {
      cur_feat
    } else {
      sel_ae_rn[1]
    }
    pObjects$memory[[plot_name]][[iSEE:::.featAssayYAxisFeatName]] <- new_feat
    updateSelectizeInput(session, feat_field,
      choices = sel_ae_rn, selected = new_feat, server = TRUE)
    iSEE:::.requestCleanUpdate(plot_name, pObjects, rObjects)
  }, ignoreInit = TRUE)
  # end

  invisible(NULL)
})

############################################################
# Data generation
############################################################

#' @importFrom iSEE .textEval
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment rowData assay colData
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
setMethod(".generateDotPlotData", "LinkedFeaturesAssayPlot", function(x, envir) {
  data_cmds <- list()

  gene_selected_y <- slot(x, iSEE:::.featAssayYAxisFeatName)
  lookup_col      <- slot(x, "LookupColumn")
  map_col         <- slot(x, "MapColumn")
  alt_exp_name    <- slot(x, "Experiment")
  alt_assay_name  <- slot(x, "AltAssay")
  sel_ae_name     <- slot(x, "SelectionExperiment")

  err_msg <- "Select other selection source/lookup column or visualisation source/Map column/assay"
  make_error_result <- function() {
    err_cmd <- paste0(
      "plot.data <- data.frame(",
      "altExp_feature_id = NA_character_, sample = NA_character_,",
      " X = NA_real_, Y = NA_real_);")
    .textEval(err_cmd, envir)
    list(commands = err_cmd, labels = list(title = err_msg, X = "", Y = ""))
  }

  # Lookup key from the selection source's rowData: either the main
  # experiment (.selectionMainExpTitle sentinel) or a named altExp.
  pg_cmd <- if (identical(lookup_col, "(rowname)")) {
    sprintf("pg <- %s;", deparse(gene_selected_y))
  } else if (identical(sel_ae_name, .selectionMainExpTitle)) {
    sprintf("pg <- rowData(se)[%s, %s];",
            deparse(gene_selected_y), deparse(lookup_col))
  } else {
    c(
      sprintf("sel_ae <- altExp(se, %s);", deparse(sel_ae_name)),
      sprintf("pg <- rowData(sel_ae)[%s, %s];",
              deparse(gene_selected_y), deparse(lookup_col))
    )
  }

  # Visualisation source: either the main experiment (.selectionMainExpTitle
  # sentinel) or a named altExp.
  ae_cmd <- if (identical(alt_exp_name, .selectionMainExpTitle)) {
    "ae <- se;"
  } else {
    sprintf("ae <- altExp(se, %s);", deparse(alt_exp_name))
  }

  data_cmds[["y"]] <- c(
    ae_cmd,
    pg_cmd,
    paste0("if (length(pg) == 0 || all(is.na(pg))) ",
           "stop('No value found in the Lookup column for the selected feature.",
           " Check the Lookup column setting.');"),
    if (identical(map_col, "(rowname)")) {
      "alt <- ae[rownames(ae) %in% pg, ];"
    } else {
      sprintf("alt <- ae[rowData(ae)[[%s]] %%in%% pg, ];", deparse(map_col))
    },
    sprintf(paste0(
      "plot.data <- as.data.frame.matrix(assay(alt, %s)) |>",
      " tibble::rownames_to_column('altExp_feature_id') |>",
      " tidyr::pivot_longer(names_to = 'sample', values_to = 'Y',",
      " -'altExp_feature_id') |>",
      " as.data.frame();"
    ), local({
      int_idx <- suppressWarnings(as.integer(alt_assay_name))
      if (!is.na(int_idx)) paste0(int_idx, "L") else deparse(alt_assay_name)
    }))
  )

  x_choice <- slot(x, iSEE:::.featAssayXAxis)
  if (x_choice == iSEE:::.featAssayXAxisColDataTitle) {
    x_lab <- slot(x, iSEE:::.featAssayXAxisColData)
    data_cmds[["x"]] <- sprintf(
      "plot.data$X <- colData(se)[plot.data$sample, %s];",
      deparse(x_lab))
  } else {
    x_lab <- ""
    data_cmds[["x"]] <- "plot.data$X <- as.factor('');"
  }

  data_cmds <- unlist(data_cmds)
  ok <- tryCatch({ .textEval(data_cmds, envir); TRUE }, error = function(e) FALSE)
  if (!ok) return(make_error_result())

  pg_value   <- paste(get("pg", envir = envir), collapse = ", ")
  plot_title <- if (pg_value == gene_selected_y) {
    pg_value
  } else {
    sprintf("%s (%s)", pg_value, gene_selected_y)
  }

  list(
    commands = data_cmds,
    labels   = list(title = plot_title, X = x_lab, Y = alt_assay_name)
  )
})

############################################################
# Plot generation
############################################################

#' @importFrom iSEE .buildAes .addFacets .addCustomLabelsCommands
#'   .addLabelCentersCommands .addMultiSelectionPlotCommands .textEval
#' @importFrom dplyr n_distinct
#' @importFrom ggplot2 ggplot aes geom_line theme_bw theme element_text
#'   scale_x_continuous scale_y_continuous labs
setMethod(".generateDotPlot", "LinkedFeaturesAssayPlot", function(x, labels, envir) {
  # Error sentinel: .generateDotPlotData signals failure via labels$Y == "".
  # Render a friendly red-title message instead of letting the standard renderer
  # hit missing columns in plot.data and freeze.
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

  plot_data <- envir$plot.data

  is_subsetted   <- exists("plot.data.all", envir = envir, inherits = FALSE)
  is_downsampled <- exists("plot.data.pre", envir = envir, inherits = FALSE)

  add_lines <- grepl("lines", slot(x, "PlotType"), fixed = TRUE)
  plot_type  <- if (slot(x, "PlotType") == "Auto") {
    envir$plot.type
  } else {
    "scatter2"
  }

  args <- list(
    plot_data,
    param_choices  = x,
    x_lab          = labels$X,
    y_lab          = labels$Y,
    color_lab      = labels$ColorBy,
    shape_lab      = labels$ShapeBy,
    size_lab       = labels$SizeBy,
    title          = labels$title,
    is_subsetted   = is_subsetted,
    is_downsampled = is_downsampled
  )

  plot_cmds <- switch(plot_type,
    square            = do.call(iSEE:::.square_plot,  args),
    violin            = do.call(iSEE:::.violin_plot,  args),
    violin_horizontal = do.call(iSEE:::.violin_plot,  c(args, list(horizontal = TRUE))),
    scatter           = do.call(iSEE:::.scatter_plot, args),
    scatter2          = do.call(.scatter_plot2,        args)
  )

  if (slot(x, iSEE:::.shapeByField) != iSEE:::.shapeByNothingTitle) {
    n_shapes <- dplyr::n_distinct(plot_data$ShapeBy)
    if (n_shapes < 26) {
      N <- length(plot_cmds)
      plot_cmds[[N]]     <- paste(plot_cmds[[N]], "+")
      plot_cmds[[N + 1]] <- sprintf("scale_shape_manual(values = seq_len(%d))", n_shapes)
    }
  }

  if (add_lines) {
    N <- length(plot_cmds)
    color_set <- !is.null(plot_data$ColorBy)
    aes_line  <- iSEE:::.buildAes(
      color = color_set, group = TRUE,
      alt   = c(color = iSEE:::.set_colorby_when_none(x)))
    aes_line <- gsub("GroupBy", "altExp_feature_id", aes_line, fixed = TRUE)
    plot_cmds[[N]] <- paste(plot_cmds[[N]], "+")
    plot_cmds[["geom_line"]] <- sprintf(
      "geom_line(data = plot.data, mapping = %s)", aes_line)
  }

  facet_cmd <- iSEE:::.addFacets(x)
  if (length(facet_cmd)) {
    N <- length(plot_cmds)
    plot_cmds[[N]] <- paste(plot_cmds[[N]], "+")
    plot_cmds <- c(plot_cmds, facet_cmd)
  }

  plot_cmds <- iSEE:::.addCustomLabelsCommands(x, commands = plot_cmds,
                                                plot_type = plot_type)

  if (plot_type == "scatter") {
    plot_cmds <- iSEE:::.addLabelCentersCommands(x, commands = plot_cmds)
  }

  plot_cmds <- iSEE:::.addMultiSelectionPlotCommands(x,
    flip     = (plot_type == "violin_horizontal"),
    envir    = envir,
    commands = plot_cmds)

  if (dplyr::n_distinct(plot_data$altExp_feature_id) > 10) {
    N <- length(plot_cmds)
    plot_cmds[[N]] <- paste(plot_cmds[[N]], "+")
    plot_cmds[["no_legend"]] <- "theme(legend.position = 'none')"
  }

  list(plot = iSEE:::.textEval(plot_cmds, envir), commands = plot_cmds)
})

############################################################
# addDotPlotData* — index colData by plot.data$sample
############################################################

#' @importFrom iSEE .textEval
#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr n_distinct
setMethod(".addDotPlotDataColor", "LinkedFeaturesAssayPlot", function(x, envir) {
  color_choice <- slot(x, iSEE:::.colorByField)

  if (color_choice == iSEE:::.colorByColDataTitle) {
    covariate_name <- slot(x, iSEE:::.colorByColData)
    label <- covariate_name
    cmds  <- if (covariate_name == "altExp_feature_id") {
      "plot.data$ColorBy <- plot.data$altExp_feature_id;"
    } else {
      sprintf("plot.data$ColorBy <- colData(se)[plot.data$sample, %s];",
              deparse(covariate_name))
    }

  } else if (color_choice == iSEE:::.colorByFeatNameTitle) {
    chosen_gene  <- slot(x, iSEE:::.colorByFeatName)
    assay_choice <- slot(x, iSEE:::.colorByFeatNameAssay)
    label <- sprintf("%s\n(%s)", chosen_gene, assay_choice)
    cmds  <- sprintf(
      "plot.data$ColorBy <- rep(assay(se, %s)[%s, ], dplyr::n_distinct(plot.data$altExp_feature_id));",
      deparse(assay_choice), deparse(chosen_gene))

  } else if (color_choice == iSEE:::.colorBySampNameTitle) {
    chosen_sample <- slot(x, iSEE:::.colorBySampName)
    label <- chosen_sample
    cmds  <- sprintf(
      "plot.data$ColorBy <- logical(nrow(plot.data));\nplot.data[plot.data$sample == %s, 'ColorBy'] <- TRUE;",
      deparse(chosen_sample))

  } else if (color_choice == iSEE:::.colorByColSelectionsTitle) {
    label  <- "Column selection"
    target <- if (exists("col_selected", envir = envir, inherits = FALSE)) {
      "col_selected"
    } else {
      "list()"
    }
    cmds <- sprintf(
      "plot.data$ColorBy <- iSEE::multiSelectionToFactor(%s, colnames(se));",
      target)

  } else {
    return(NULL)
  }

  iSEE:::.textEval(cmds, envir)
  list(commands = cmds, labels = list(ColorBy = label))
})

#' @importFrom iSEE .textEval
#' @importFrom SummarizedExperiment colData
setMethod(".addDotPlotDataShape", "LinkedFeaturesAssayPlot", function(x, envir) {
  shape_choice   <- slot(x, iSEE:::.shapeByField)
  covariate_name <- slot(x, iSEE:::.shapeByColData)

  if (shape_choice == iSEE:::.shapeByColDataTitle) {
    if (covariate_name == "altExp_feature_id") {
      label <- "altExp_feature_id"
      cmds  <- "plot.data$ShapeBy <- plot.data$altExp_feature_id;"
    } else if (covariate_name %in% colnames(colData(envir$se))) {
      label <- covariate_name
      cmds  <- sprintf("plot.data$ShapeBy <- colData(se)[plot.data$sample, %s];",
                       deparse(covariate_name))
    } else {
      label <- "altExp_feature_id"
      cmds  <- "plot.data$ShapeBy <- plot.data$altExp_feature_id;"
    }
  } else {
    label <- "altExp_feature_id"
    cmds  <- "plot.data$ShapeBy <- plot.data$altExp_feature_id;"
  }

  iSEE:::.textEval(cmds, envir)
  list(commands = cmds, labels = list(ShapeBy = label))
})

#' @importFrom iSEE .textEval
#' @importFrom SummarizedExperiment colData
setMethod(".addDotPlotDataSize", "LinkedFeaturesAssayPlot", function(x, envir) {
  size_choice <- slot(x, iSEE:::.sizeByField)

  if (size_choice == iSEE:::.sizeByColDataTitle) {
    covariate_name <- slot(x, iSEE:::.sizeByColData)
    label <- covariate_name
    cmds  <- sprintf("plot.data$SizeBy <- colData(se)[plot.data$sample, %s];",
                     deparse(covariate_name))
  } else {
    return(NULL)
  }

  iSEE:::.textEval(cmds, envir)
  list(commands = cmds, labels = list(SizeBy = label))
})

#' @importFrom iSEE .textEval
#' @importFrom SummarizedExperiment colData
setMethod(".addDotPlotDataFacets", "LinkedFeaturesAssayPlot", function(x, envir) {
  facet_cmds <- NULL
  labels      <- list()

  params <- list(
    list(iSEE:::.facetRow,    "FacetRow",    iSEE:::.facetRowByColData),
    list(iSEE:::.facetColumn, "FacetColumn", iSEE:::.facetColumnByColData)
  )

  for (f in seq_len(2)) {
    current     <- params[[f]]
    param_field <- current[[1]]
    pd_field    <- current[[2]]
    facet_mode  <- slot(x, param_field)

    if (facet_mode == iSEE:::.facetByColDataTitle) {
      facet_data <- x[[current[[3]]]]
      facet_cmds[pd_field] <- sprintf(
        "plot.data$%s <- colData(se)[plot.data$sample, %s];",
        pd_field, deparse(facet_data))
      labels[[pd_field]] <- facet_data

    } else if (facet_mode == iSEE:::.facetByColSelectionsTitle) {
      target <- if (exists("col_selected", envir = envir, inherits = FALSE)) {
        "col_selected"
      } else {
        "list()"
      }
      facet_cmds[pd_field] <- sprintf(
        "plot.data$%s <- iSEE::multiSelectionToFactor(%s, colnames(se));",
        pd_field, target)
      labels[[pd_field]] <- "Column selection"
    }
  }

  iSEE:::.textEval(facet_cmds, envir)
  list(commands = facet_cmds, labels = labels)
})

############################################################
# Visual UI overrides
############################################################

# Return colData choices directly from se rather than relying on the
# ColumnDotPlot cache, so colour-by UI is always populated.
#' @importFrom iSEE .allowableColorByDataChoices .findAtomicFields
#' @importFrom SummarizedExperiment colData
setMethod(".allowableColorByDataChoices", "LinkedFeaturesAssayPlot", function(x, se) {
  .findAtomicFields(colData(se))
})

#' @importFrom iSEE .getCachedCommonInfo .getEncodedName .allowableColorByDataChoices
#'   .singleSelectionDimension .radioButtons.iSEE .conditionalOnRadio
#'   .sliderInput.iSEE
#' @importFrom colourpicker colourInput
#' @importFrom shiny hr tagList selectInput selectizeInput checkboxInput
setMethod(".defineVisualColorInterface", "LinkedFeaturesAssayPlot",
    function(x, se, select_info) {
  all_assays <- .getCachedCommonInfo(se, "LinkedFeaturesAssayPlot")$valid.assay.names

  plot_name      <- iSEE:::.getEncodedName(x)
  colorby_field  <- paste0(plot_name, "_", iSEE:::.colorByField)

  colorby          <- iSEE:::.getDotPlotColorConstants(x)
  mydim_single     <- iSEE:::.singleSelectionDimension(x)
  otherdim_single  <- setdiff(c("feature", "sample"), mydim_single)
  mydim_choices    <- select_info[[mydim_single]]
  otherdim_choices <- select_info[[otherdim_single]]

  covariates    <- iSEE:::.allowableColorByDataChoices(x, se)
  color_choices <- iSEE:::.colorByNothingTitle
  if (length(covariates)) {
    color_choices <- c(color_choices, iSEE:::.colorByColDataTitle)
  }

  tagList(
    hr(),
    iSEE:::.radioButtons.iSEE(x, iSEE:::.colorByField,
      label    = "Color by:",
      inline   = TRUE,
      choices  = color_choices,
      selected = slot(x, iSEE:::.colorByField)),
    iSEE:::.conditionalOnRadio(colorby_field, iSEE:::.colorByNothingTitle,
      colourpicker::colourInput(
        paste0(plot_name, "_", iSEE:::.colorByDefaultColor),
        label = NULL,
        value = slot(x, iSEE:::.colorByDefaultColor))),
    iSEE:::.conditionalOnRadio(colorby_field, colorby$metadata$title,
      selectInput(
        paste0(plot_name, "_", colorby$metadata$field),
        label    = NULL,
        choices  = c("altExp_feature_id", iSEE:::.allowableColorByDataChoices(x, se)),
        selected = x[[colorby$metadata$field]])),
    iSEE:::.conditionalOnRadio(colorby_field, colorby$name$title,
      selectizeInput(
        paste0(plot_name, "_", colorby$name$field),
        label = NULL, selected = NULL, choices = NULL, multiple = FALSE),
      selectInput(
        paste0(plot_name, "_", colorby$name$table),
        label    = NULL,
        choices  = mydim_choices,
        selected = iSEE:::.choose_link(x[[colorby$name$table]], mydim_choices)),
      colourpicker::colourInput(
        paste0(plot_name, "_", colorby$name$color),
        label = NULL,
        value = x[[colorby$name$color]]),
      checkboxInput(
        paste0(plot_name, "_", colorby$name$dynamic),
        label = sprintf("Use dynamic %s selection", mydim_single),
        value = x[[colorby$name$dynamic]])),
    iSEE:::.conditionalOnRadio(colorby_field, colorby$assay$title,
      selectizeInput(
        paste0(plot_name, "_", colorby$assay$field),
        label = NULL, choices = NULL, selected = NULL, multiple = FALSE),
      selectInput(
        paste0(plot_name, "_", colorby$assay$assay),
        label    = NULL,
        choices  = all_assays,
        selected = x[[colorby$assay$assay]]),
      selectInput(
        paste0(plot_name, "_", colorby$assay$table),
        label    = NULL,
        choices  = otherdim_choices,
        selected = iSEE:::.choose_link(x[[colorby$assay$table]], otherdim_choices)),
      checkboxInput(
        paste0(plot_name, "_", colorby$assay$dynamic),
        label = sprintf("Use dynamic %s selection", otherdim_single),
        value = x[[colorby$assay$dynamic]])),
    iSEE:::.sliderInput.iSEE(x, iSEE:::.selectTransAlpha,
      label = "Unselected point opacity:",
      min   = 0, max = 1,
      value = slot(x, iSEE:::.selectTransAlpha))
  )
})

# Shape interface: adds "altExp_feature_id" as an option in addition to
# the standard discrete colData covariates.
#' @importFrom iSEE .getEncodedName .radioButtons.iSEE .conditionalOnRadio
#' @importFrom SummarizedExperiment colData
#' @importFrom shiny hr tagList selectInput
setMethod(".defineVisualShapeInterface", "LinkedFeaturesAssayPlot", function(x, se) {
  cd_cls <- sapply(colData(se), class)
  discrete_covariates <- c(
    "altExp_feature_id",
    names(cd_cls)[cd_cls %in% c("factor", "character")])

  if (length(discrete_covariates)) {
    plot_name     <- iSEE:::.getEncodedName(x)
    shapeby_field <- paste0(plot_name, "_", iSEE:::.shapeByField)
    shapeby       <- iSEE:::.getDotPlotShapeConstants(x)

    tagList(
      hr(),
      iSEE:::.radioButtons.iSEE(x, iSEE:::.shapeByField,
        label    = "Shape by:",
        inline   = TRUE,
        choices  = c(iSEE:::.shapeByNothingTitle, shapeby$metadata$title),
        selected = slot(x, iSEE:::.shapeByField)),
      iSEE:::.conditionalOnRadio(shapeby_field, shapeby$metadata$title,
        selectInput(
          paste0(plot_name, "_", shapeby$metadata$field),
          label    = NULL,
          choices  = discrete_covariates,
          selected = x[[shapeby$metadata$field]]))
    )
  } else {
    NULL
  }
})

############################################################
# Panel tour
############################################################

#' @importFrom iSEE .getEncodedName .getPanelColor .addTourStep .dataParamBoxOpen
#'   .definePanelTour
setMethod(".definePanelTour", "LinkedFeaturesAssayPlot", function(x) {
  collated <- rbind(
    c(
      element = paste0("#", .getEncodedName(x)),
      intro   = sprintf(
        "The <font color=\"%s\">Linked features assay plot</font> panel visualises
assay values from one source — an alternative experiment, or the main
experiment — while using a <em>different</em> source for feature selection.
<br><br>
A feature is selected from the <em>Selection source</em> (typically via an
<em>AltExp row data table</em>, or any panel transmitting a row selection for
the main experiment).  The value in the <em>Lookup column</em> of that
source's <code>rowData</code> is used as a join key to retrieve matching rows
from the <em>Visualisation source</em> via the <em>Map column</em>.  All
matching features are displayed as separate traces per sample.
<br><br>
A typical use case: select a protein from the main experiment (or a
protein-level altExp) and display all peptides or precursors from a
peptide-level altExp that map to it.",
        .getPanelColor(x)
      )
    ),
    c(
      element = paste0("#", .getEncodedName(x)),
      intro   = "Each point represents one
<strong>(visualisation-altExp feature, sample)</strong> combination.  When
multiple rows match — for example several peptides mapping to the same protein
— all appear simultaneously, distinguishable by colour, shape, or a connecting
line."
    ),
    .addTourStep(x, .dataParamBoxOpen,
      "The <i>Data parameters</i> box contains all controls specific to this
panel.<br><br><strong>Action:</strong> click to expand.")
  )

  parent_tour <- callNextMethod()
  parent_tour <- parent_tour[
    !grepl("Feature assay plot", parent_tour$intro, fixed = TRUE), ]

  rbind(
    data.frame(element = collated[, 1], intro = collated[, 2],
               stringsAsFactors = FALSE),
    parent_tour
  )
})
