############################################################
# AltExpRowDataTable
############################################################

# Extends RowTable (the same parent as RowDataTable) to display rowData from
# an alternative experiment instead of the main experiment.
# Row selections propagate altExp feature names; this is compatible with other
# altExp-aware panels but not with main-experiment-only panels.

#' AltExpRowDataTable: row data table for alternative experiments
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
#' @slot AltExp \code{character(1)}.  Name of the alternative experiment whose
#'   \code{rowData} is displayed.  Defaults to \code{NA}, resolved at runtime
#'   to the first entry of \code{altExpNames(se)}.
#'
#' @section Inherited slots:
#' All slots from \code{\linkS4class{RowTable}} are inherited, including
#' \code{Selected}, \code{Search}, \code{SearchColumns}, and
#' \code{HiddenColumns}.
#'
#' @seealso
#' \code{\link{AltExpRowDataTable}} for the constructor.
#' \code{\linkS4class{RowTable}} for the parent class.
#' \code{\linkS4class{RowDataTable}} for the analogous main-experiment table.
#'
#' @name AltExpRowDataTable-class
#' @rdname AltExpRowDataTable-class
#' @exportClass AltExpRowDataTable
NULL

setClass("AltExpRowDataTable",
  contains = "RowTable",
  slots = c(AltExp = "character"))

#' @importFrom iSEE .emptyDefault
setMethod("initialize", "AltExpRowDataTable", function(.Object, ...) {
  args <- list(...)
  args <- .emptyDefault(args, "AltExp", NA_character_)
  do.call(callNextMethod, c(list(.Object), args))
})

#' Construct an AltExpRowDataTable panel
#'
#' Creates an instance of \code{\linkS4class{AltExpRowDataTable}} for use as
#' a panel in an iSEE application.
#'
#' @param ... Named arguments corresponding to slots of
#'   \code{\linkS4class{AltExpRowDataTable}} or its parent classes.
#'
#' @return An \code{\linkS4class{AltExpRowDataTable}} object.
#'
#' @seealso \code{\linkS4class{AltExpRowDataTable}} for slot details.
#' @export
AltExpRowDataTable <- function(...) new("AltExpRowDataTable", ...)

setMethod(".fullName",   "AltExpRowDataTable", function(x) "AltExp row data table")
setMethod(".panelColor", "AltExpRowDataTable", function(x) "#AA5500")

#' @importFrom iSEE .getCachedCommonInfo .setCachedCommonInfo
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom SummarizedExperiment rowData
setMethod(".cacheCommonInfo", "AltExpRowDataTable", function(x, se) {
  if (!is.null(.getCachedCommonInfo(se, "AltExpRowDataTable"))) return(se)
  se <- callNextMethod()

  valid_by_ae <- lapply(altExpNames(se), function(ae_name) {
    iSEE:::.findAtomicFields(rowData(altExp(se, ae_name)))
  })
  names(valid_by_ae) <- altExpNames(se)

  .setCachedCommonInfo(se, "AltExpRowDataTable",
    valid.altExp.names            = altExpNames(se),
    valid.rowData.names.by.altExp = valid_by_ae)
})

#' @importFrom iSEE .replaceMissingWithFirst
#' @importFrom SingleCellExperiment altExp altExpNames
setMethod(".refineParameters", "AltExpRowDataTable", function(x, se) {
  # Save before callNextMethod() may validate Selected against rownames(se)
  # and reset an altExp feature name to "" or the first main-experiment row.
  saved_selected <- slot(x, "Selected")

  x <- callNextMethod()
  if (is.null(x)) return(NULL)

  x <- .replaceMissingWithFirst(x, "AltExp", altExpNames(se))

  ae_rn <- rownames(altExp(se, slot(x, "AltExp")))
  if (nzchar(saved_selected) && saved_selected %in% ae_rn) {
    slot(x, "Selected") <- saved_selected
  }

  x
})

#' @importFrom iSEE .getEncodedName .selectInput.iSEE
#' @importFrom SingleCellExperiment altExpNames
#' @importFrom shiny selectInput
setMethod(".defineDataInterface", "AltExpRowDataTable", function(x, se, select_info) {
  panel_name  <- .getEncodedName(x)
  all_altexps <- altExpNames(se)

  c(
    list(
      selectInput(
        paste0(panel_name, "_AltExp"),
        label    = "AltExp:",
        choices  = all_altexps,
        selected = iSEE:::.choose_link(slot(x, "AltExp"), all_altexps)
      )
    ),
    callNextMethod()   # hidden-columns selector from RowTable
  )
})

#' @importFrom iSEE .getEncodedName .createUnprotectedParameterObservers
#' @importFrom SingleCellExperiment altExpNames
setMethod(".createObservers", "AltExpRowDataTable",
    function(x, se, input, session, pObjects, rObjects) {
  callNextMethod()   # Table observers + RowTable row-propagation observer

  .createUnprotectedParameterObservers(
    .getEncodedName(x),
    fields   = "AltExp",
    input    = input,
    pObjects = pObjects,
    rObjects = rObjects)

  invisible(NULL)
})

#' @importFrom iSEE .getCachedCommonInfo .textEval
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment rowData
setMethod(".generateTable", "AltExpRowDataTable", function(x, envir) {
  ae_name <- slot(x, "AltExp")

  cmds <- c(
    sprintf("ae  <- altExp(se, %s);", deparse(ae_name)),
    "tab <- as.data.frame(rowData(ae));"
  )

  if (exists("row_selected", envir = envir, inherits = FALSE)) {
    cmds <- c(cmds,
      "tab <- tab[rownames(tab) %in% unlist(row_selected), , drop = FALSE]")
  }

  valid_names <- .getCachedCommonInfo(
    envir$se, "AltExpRowDataTable")$valid.rowData.names.by.altExp[[ae_name]]

  if (!identical(colnames(rowData(altExp(envir$se, ae_name))), valid_names)) {
    cmds <- c(cmds, sprintf(
      "tab <- tab[, %s, drop = FALSE]",
      paste(deparse(valid_names), collapse = "\n     ")))
  }

  .textEval(cmds, envir)
  cmds
})

#' @importFrom iSEE .getEncodedName .getPanelColor .addTourStep .dataParamBoxOpen
setMethod(".definePanelTour", "AltExpRowDataTable", function(x) {
  collated <- rbind(
    c(
      element = paste0("#", .getEncodedName(x)),
      intro   = sprintf(
        "The <font color=\"%s\">AltExp row data table</font> displays the
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
