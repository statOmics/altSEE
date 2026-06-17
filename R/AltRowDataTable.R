############################################################
# AltRowDataTable
############################################################

# Extends RowTable (the same parent as RowDataTable) to display rowData from
# an alternative experiment or the main experiment (via the "(Main)"
# sentinel).  Row selections propagate the chosen experiment's feature names.

#' AltRowDataTable: row data table for alternative experiments
#'
#' An S4 class extending \code{\linkS4class{RowTable}} to display the
#' \code{\link[SummarizedExperiment]{rowData}} of an
#' \emph{alternative experiment}
#' (\code{\link[SingleCellExperiment]{altExp}}) stored inside a
#' \code{\linkS4class{SingleCellExperiment}}.
#' Each row of the table corresponds to a feature of the selected altExp;
#' clicking a row transmits that altExp feature name as a row selection to
#' other panels.
#'
#' @slot Experiment \code{character(1)}.  Name of the alternative experiment whose
#'   \code{rowData} is displayed, or the sentinel \code{"(Main)"}
#'   to display the main experiment's \code{rowData}.  Defaults to \code{NA},
#'   resolved at runtime to the first available alternative experiment.
#'
#' @section Inherited slots:
#' All slots from \code{\linkS4class{RowTable}} are inherited, including
#' \code{Selected}, \code{Search}, \code{SearchColumns}, and
#' \code{HiddenColumns}.
#'
#' @seealso
#' \code{\link{AltRowDataTable}} for the constructor.
#' \code{\linkS4class{RowTable}} for the parent class.
#' \code{\linkS4class{RowDataTable}} for the analogous main-experiment table.
#'
#' @name AltRowDataTable-class
#' @rdname AltRowDataTable-class
#' @exportClass AltRowDataTable
NULL

setClass("AltRowDataTable",
  contains = "RowTable",
  slots = c(Experiment = "character"))

#' @importFrom iSEE .emptyDefault
setMethod("initialize", "AltRowDataTable", function(.Object, ...) {
  args <- list(...)
  args <- .emptyDefault(args, "Experiment", NA_character_)
  do.call(callNextMethod, c(list(.Object), args))
})

#' Construct an AltRowDataTable panel
#'
#' Creates an instance of \code{\linkS4class{AltRowDataTable}} for use as
#' a panel in an iSEE application.
#'
#' @param ... Named arguments corresponding to slots of
#'   \code{\linkS4class{AltRowDataTable}} or its parent classes.
#'
#' @return An \code{\linkS4class{AltRowDataTable}} object.
#'
#' @seealso \code{\linkS4class{AltRowDataTable}} for slot details.
#' @export
AltRowDataTable <- function(...) new("AltRowDataTable", ...)

setMethod(".fullName",   "AltRowDataTable", function(x) "Alt row data table")
setMethod(".panelColor", "AltRowDataTable", function(x) "#AA5500")

#' @importFrom iSEE .getCachedCommonInfo .setCachedCommonInfo
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment rowData
setMethod(".cacheCommonInfo", "AltRowDataTable", function(x, se) {
  if (!is.null(.getCachedCommonInfo(se, "AltRowDataTable"))) return(se)
  se <- callNextMethod()

  all_names <- c(altExpNames(se), .selectionMainExpTitle)
  valid_by_ae <- lapply(all_names, function(ae_name) {
    ae_obj <- if (identical(ae_name, .selectionMainExpTitle)) se else altExp(se, ae_name)
    iSEE:::.findAtomicFields(rowData(ae_obj))
  })
  names(valid_by_ae) <- all_names

  .setCachedCommonInfo(se, "AltRowDataTable",
    valid.altExp.names            = all_names,
    valid.rowData.names.by.altExp = valid_by_ae)
})

#' @importFrom iSEE .replaceMissingWithFirst
#' @importFrom SingleCellExperiment altExp altExpNames
setMethod(".refineParameters", "AltRowDataTable", function(x, se) {
  # Save before callNextMethod() may validate Selected against rownames(se)
  # and reset an altExp feature name to "" or the first main-experiment row.
  saved_selected <- slot(x, "Selected")

  x <- callNextMethod()
  if (is.null(x)) return(NULL)

  x <- .replaceMissingWithFirst(x, "Experiment", c(altExpNames(se), .selectionMainExpTitle))

  ae_name <- slot(x, "Experiment")
  ae_rn <- if (identical(ae_name, .selectionMainExpTitle)) rownames(se) else rownames(altExp(se, ae_name))
  if (nzchar(saved_selected) && saved_selected %in% ae_rn) {
    slot(x, "Selected") <- saved_selected
  }

  x
})

#' @importFrom iSEE .getEncodedName .selectInput.iSEE
#' @importFrom SingleCellExperiment altExpNames
#' @importFrom shiny selectInput
setMethod(".defineDataInterface", "AltRowDataTable", function(x, se, select_info) {
  panel_name  <- .getEncodedName(x)
  all_altexps <- c(altExpNames(se), .selectionMainExpTitle)

  c(
    list(
      selectInput(
        paste0(panel_name, "_Experiment"),
        label    = "Experiment:",
        choices  = all_altexps,
        selected = iSEE:::.choose_link(slot(x, "Experiment"), all_altexps)
      )
    ),
    callNextMethod()   # hidden-columns selector from RowTable
  )
})

#' @importFrom iSEE .getEncodedName .createUnprotectedParameterObservers
#' @importFrom SingleCellExperiment altExpNames
setMethod(".createObservers", "AltRowDataTable",
    function(x, se, input, session, pObjects, rObjects) {
  callNextMethod()   # Table observers + RowTable row-propagation observer

  .createUnprotectedParameterObservers(
    .getEncodedName(x),
    fields   = "Experiment",
    input    = input,
    pObjects = pObjects,
    rObjects = rObjects)

  invisible(NULL)
})

#' @importFrom iSEE .getCachedCommonInfo .textEval
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment rowData
setMethod(".generateTable", "AltRowDataTable", function(x, envir) {
  ae_name <- slot(x, "Experiment")

  cmds <- c(
    if (identical(ae_name, .selectionMainExpTitle)) "ae <- se;" else
      sprintf("ae  <- altExp(se, %s);", deparse(ae_name)),
    "tab <- as.data.frame(rowData(ae));"
  )

  if (exists("row_selected", envir = envir, inherits = FALSE)) {
    cmds <- c(cmds,
      "tab <- tab[rownames(tab) %in% unlist(row_selected), , drop = FALSE]")
  }

  valid_names <- .getCachedCommonInfo(
    envir$se, "AltRowDataTable")$valid.rowData.names.by.altExp[[ae_name]]

  ae_obj <- if (identical(ae_name, .selectionMainExpTitle)) envir$se else altExp(envir$se, ae_name)
  if (!identical(colnames(rowData(ae_obj)), valid_names)) {
    cmds <- c(cmds, sprintf(
      "tab <- tab[, %s, drop = FALSE]",
      paste(deparse(valid_names), collapse = "\n     ")))
  }

  .textEval(cmds, envir)
  cmds
})

#' @importFrom iSEE .getEncodedName .getPanelColor .addTourStep .dataParamBoxOpen
setMethod(".definePanelTour", "AltRowDataTable", function(x) {
  collated <- rbind(
    c(
      element = paste0("#", .getEncodedName(x)),
      intro   = sprintf(
        "The <font color=\"%s\">Alt row data table</font> displays the
<code>rowData</code> of a selected alternative experiment — for example all
peptides and their annotations in an MS-based proteomics workflow, or all
antibody-derived tags in a CITE-seq experiment.
Clicking a row transmits that altExp feature as a row selection to other
panels.", .getPanelColor(x))
    ),
    .addTourStep(x, .dataParamBoxOpen,
      "The <i>Data parameters</i> box lets you choose which alternative
experiment to display and which columns to hide.<br><br>
<strong>Action:</strong> click to expand.")
  )

  rbind(
    data.frame(element = collated[, 1], intro = collated[, 2],
               stringsAsFactors = FALSE),
    callNextMethod()
  )
})
