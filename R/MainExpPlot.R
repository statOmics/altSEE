############################################################
# MainExpPlot — modified from FeatureAssayPlot
############################################################

# Extends FeatureAssayPlot from iSEE to add line-plot support for the main
# experiment assay.  Shares the .scatter_plot2 renderer with AltExpPlot

#' MainExpPlot: FeatureAssayPlot with optional line-plot support
#'
#' An S4 class that extends \code{\linkS4class{FeatureAssayPlot}} to add a
#' \code{PlotType} slot, allowing the user to overlay a \code{geom_line} layer
#' on top of the standard scatter plot.  This is useful when the x-axis
#' encodes an ordered variable such as time, pseudotime or an ordered sampleId 
#' e.g. according to treatement conditions and you want to connect measurements 
#' for the same feature across samples.
#'
#' @slot PlotType \code{character(1)}.  One of \code{"Auto"},
#'   \code{"Scatter"}, or \code{"Scatter + lines"}.  \code{"Auto"} delegates
#'   the renderer choice to iSEE's internal heuristic; \code{"Scatter"} always
#'   uses a dot plot; \code{"Scatter + lines"} adds a \code{geom_line} layer
#'   connecting all points (grouped by the selected feature via the \code{id}
#'   column in \code{plot.data}).
#'
#' @section Inherited slots:
#' All slots from \code{\linkS4class{FeatureAssayPlot}} are inherited,
#' including those controlling assay choice, x-axis, colour, shape, size,
#' faceting, and multi-selection behaviour.
#'
#' @seealso
#' \code{\link{MainExpPlot}} for the constructor.
#' \code{\linkS4class{FeatureAssayPlot}} for the parent class.
#' \code{\linkS4class{AltExpPlot}} for the analogous panel that plots data
#' from an alternative experiment.
#'
#' @name MainExpPlot-class
#' @rdname MainExpPlot-class
#' @exportClass MainExpPlot

setClass("MainExpPlot", contains="FeatureAssayPlot",
         slots=c(PlotType = "character"))

#' @importFrom iSEE .validStringError
setValidity2("MainExpPlot", function(object) {
  msg <- character(0)
  msg <- .validStringError(msg, object, "PlotType")
  if (length(msg)) return(msg)
  TRUE
})

#' @importFrom iSEE .emptyDefault
setMethod("initialize", "MainExpPlot", function(.Object, ...) {
  args <- list(...)
  args <- .emptyDefault(args, "PlotType", "Auto")
  do.call(callNextMethod, c(list(.Object), args))
})

#' Construct a MainExpPlot panel
#'
#' Creates an instance of the \code{\linkS4class{MainExpPlot}} class for use
#' as a panel in an iSEE application.
#'
#' @param ... Named arguments corresponding to slots of
#'   \code{\linkS4class{MainExpPlot}} or its parent classes.  Any slot not
#'   supplied receives a sensible default: \code{PlotType} defaults to
#'   \code{"Auto"}.
#'
#' @return A \code{\linkS4class{MainExpPlot}} object.
#'
#' @examples
#' \dontrun{
#' # Minimal panel using runtime defaults
#' MainExpPlot()
#'
#' # Pre-configure to show scatter plot with lines connecting points
#' MainExpPlot(PlotType = "Scatter + lines")
#' 
#' # Pre-configure YAxisFeatureSource, PlotType, XAxis and XAxisColumnData
#' MainExpPlot(YAxisFeatureSource = "RowDataTable1", 
#' PlotType = "Scatter + lines", 
#' XAxis = "Column data",
#' XAxisColumnData = "sampleId")
#' 
#' # Example of iSEE instance with mainExp and altExp
#' # linked via the "Proteins" in rowData columns of MainExp and AltExp 
#' data("proteinsPeptidesSce")
#' iSEE(
#' sceProteinsPeptides,
#' initial = list(
#' RowDataTable(),
#' MainExpPlot(YAxisFeatureSource = "RowDataTable1", 
#'             PlotType = "Scatter + lines", 
#'             XAxis = "Column data",
#'             XAxisColumnData = "sampleId"),
#' AltExpPlot(YAxisFeatureSource = "RowDataTable1", 
#'            AltAssay = "peptides_norm",
#'            PlotType = "Scatter + lines", 
#'            XAxis = "Column data", 
#'            XAxisColumnData = "sampleId",
#'            LookupColumn = "Proteins",
#'            MapColumn = "Proteins")
#' )
#' )
#' }
#'
#' @seealso \code{\linkS4class{MainExpPlot}} for a description of all slots.
#' @export
#' 

MainExpPlot <- function(...) {
  new("MainExpPlot", ...)
}

setMethod(".fullName", "MainExpPlot", function(x) "MainExp plot")

setMethod(".panelColor", "MainExpPlot", function(x) "#7BB854")


#' @importFrom iSEE .getEncodedName .getCachedCommonInfo .choose_link
#'   .selectizeInput.iSEE .radioButtons.iSEE .conditionalOnRadio
#'   .featAssayYAxisFeatName .featAssayYAxisRowTable .featAssayYAxisFeatDynamic
#'   .featAssayXAxis .featAssayXAxisNothingTitle .featAssayXAxisColDataTitle
#'   .featAssayXAxisColData .featAssayXAxisFeatNameTitle .featAssayXAxisFeatName
#'   .featAssayXAxisRowTable .featAssayXAxisFeatDynamic
#' @importFrom SummarizedExperiment assayNames
#' @importFrom shiny selectInput checkboxInput selectizeInput

setMethod(".defineDataInterface", "MainExpPlot", function(x, se, select_info) {
  panel_name <- .getEncodedName(x)
  .input_FUN <- function(field) { paste0(panel_name, "_", field) }

  tab_by_row <- select_info$single$feature

  all_assays <- assayNames(se)

  column_covariates <- .getCachedCommonInfo(se, "ColumnDotPlot")$valid.colData.names
  xaxis_choices <- c(iSEE:::.featAssayXAxisNothingTitle)
  if (length(column_covariates)) { 
    xaxis_choices <- c(xaxis_choices, iSEE:::.featAssayXAxisColDataTitle)
  }

  list(
    .selectizeInput.iSEE(
      x, iSEE:::.featAssayYAxisFeatName,
      label = "Feature in main assay:",
      choices = NULL,
      selected = NULL,
      multiple = FALSE),
    selectInput(
      .input_FUN(iSEE:::.featAssayYAxisRowTable), 
      label=NULL, 
      choices = tab_by_row,
      selected = iSEE:::.choose_link(slot(x, iSEE:::.featAssayYAxisRowTable), 
                                   tab_by_row)
      ),
    checkboxInput(
      .input_FUN(iSEE:::.featAssayYAxisFeatDynamic),
      label = "Use dynamic feature selection for the y-axis",
      value = slot(x, iSEE:::.featAssayYAxisFeatDynamic)
      ),
    selectInput(
      .input_FUN("Assay"), 
      label = "Assay:", 
      choices = all_assays,
      selected = iSEE:::.choose_link(slot(x, "Assay"), all_assays)
      ),
    selectInput(
      .input_FUN("PlotType"), 
      label = "Plot type:",
      choices = c("Auto","Scatter", "Scatter + lines"), 
      selected = slot(x,"PlotType")
      ),
    .radioButtons.iSEE(
      x, iSEE:::.featAssayXAxis,
      label="X-axis:",
      inline=TRUE,
      choices=xaxis_choices,
      selected=slot(x, iSEE:::.featAssayXAxis)
      ),
    .conditionalOnRadio(
      .input_FUN(iSEE:::.featAssayXAxis),
      iSEE:::.featAssayXAxisColDataTitle,
      selectInput(.input_FUN(iSEE:::.featAssayXAxisColData),
                  label = "X-axis column data:",
                  choices = column_covariates, 
                  selected = slot(x, iSEE:::.featAssayXAxisColData)
                  )
      ),
    .conditionalOnRadio(
      .input_FUN(iSEE:::.featAssayXAxis),
      iSEE:::.featAssayXAxisFeatNameTitle,
      selectizeInput(.input_FUN(iSEE:::.featAssayXAxisFeatName),
                     label = "X-axis feature:", 
                     choices = NULL, 
                     selected = NULL, 
                     multiple = FALSE
                     ),
      selectInput(.input_FUN(iSEE:::.featAssayXAxisRowTable), 
                  label = NULL,
                  choices = tab_by_row, 
                  selected = slot(x, iSEE:::.featAssayXAxisRowTable)
                  ),
      checkboxInput(.input_FUN(iSEE:::.featAssayXAxisFeatDynamic),
                    label = "Use dynamic feature selection for the x-axis",
                    value = slot(x, iSEE:::.featAssayXAxisFeatDynamic)
                    )
      )
    )
  })

#' @importFrom iSEE .getEncodedName .createUnprotectedParameterObservers
setMethod(".createObservers", 
          "MainExpPlot", 
          function(x, se, input, session, pObjects, rObjects) {
  callNextMethod()

  plot_name <- .getEncodedName(x)

  .createUnprotectedParameterObservers(
    plot_name,
    fields = c("PlotType"),
    input = input, 
    pObjects = pObjects, 
    rObjects = rObjects)

  invisible(NULL)
})


#' @importFrom iSEE .featAssayYAxisFeatName .featAssayAssay .featAssayXAxis
#'   .featAssayXAxisColDataTitle .featAssayXAxisFeatNameTitle
#'   .featAssayXAxisFeatName .featAssayXAxisSelectionsTitle .textEval
#' @importFrom SummarizedExperiment assay colData

setMethod(".generateDotPlotData", "MainExpPlot", function(x, envir) {
  data_cmds <- list()

  ## Setting up the y-axis:
  gene_selected_y <- slot(x, iSEE:::.featAssayYAxisFeatName)
  assay_choice <- slot(x, iSEE:::.featAssayAssay)
  plot_title <- gene_selected_y
  y_lab <- sprintf("%s (%s)", gene_selected_y, assay_choice)
  data_cmds[["y"]] <- sprintf(
    "plot.data <- data.frame(Y=assay(se, %s)[%s, ], row.names=colnames(se))",
    deparse(assay_choice), deparse(gene_selected_y)
  )
  data_cmds[["id"]] <-  sprintf("plot.data$id <- %s",deparse(gene_selected_y))

  ## Checking X axis choice:
  x_choice <- slot(x, iSEE:::.featAssayXAxis)

  if (x_choice == iSEE:::.featAssayXAxisColDataTitle) { 
    # colData column selected
    x_lab <- slot(x, iSEE:::.featAssayXAxisColData)
    plot_title <- paste(plot_title, "vs", x_lab)
    data_cmds[["x"]] <- sprintf(
      "plot.data$X <- colData(se)[, %s];", 
      deparse(x_lab)
      )
    } else 
      if (x_choice == iSEE:::.featAssayXAxisFeatNameTitle) { # gene selected
        gene_selected_x <- slot(x, iSEE:::.featAssayXAxisFeatName)
        plot_title <- paste(plot_title, "vs", gene_selected_x)
        x_lab <- sprintf("%s (%s)", gene_selected_x, assay_choice)
        data_cmds[["x"]] <- sprintf(
          "plot.data$X <- assay(se, %s)[%s, ];",
          deparse(assay_choice), 
          deparse(gene_selected_x)
          )
        } else if (x_choice == iSEE:::.featAssayXAxisSelectionsTitle) {
          x_lab <- "Column selection"
          plot_title <- paste(plot_title, "vs column selection")
          
          if (exists("col_selected", envir=envir, inherits=FALSE)) {
            target <- "col_selected"
            } else {
              target <- "list()"
              }
          data_cmds[["x"]] <- sprintf(
            "plot.data$X <- iSEE::multiSelectionToFactor(%s, colnames(se));",
            target
            )
          } else { # no x axis variable specified: show single violin
            x_lab <- ''
            data_cmds[["x"]] <- "plot.data$X <- factor(character(ncol(se)))"
            }
  data_cmds <- unlist(data_cmds)
  .textEval(data_cmds, envir)
  
  list(commands = data_cmds, 
       labels = list(title = plot_title, X = x_lab, Y = y_lab)
       )
})


############################################################
# Plot generation
############################################################


#' @importFrom iSEE .shapeByField .shapeByNothingTitle .buildAes
#'   .set_colorby_when_none .addFacets .addCustomLabelsCommands
#'   .addLabelCentersCommands .addMultiSelectionPlotCommands .textEval
#'   .square_plot .violin_plot .scatter_plot
#' @importFrom dplyr n_distinct
#' @importFrom ggplot2 geom_line

setMethod(".generateDotPlot", "MainExpPlot", function(x, labels, envir) {
  plot_data <- envir$plot.data

  is_subsetted <- exists("plot.data.all", envir=envir, inherits=FALSE)
  is_downsampled <- exists("plot.data.pre", envir=envir, inherits=FALSE)
  
  # PlotType "Auto" defers to whatever iSEE chose in plot.type;
  # "Scatter" and "Scatter + lines" both use scatter2 as the base renderer.
  add_lines <- grepl("lines", slot(x, "PlotType"), fixed = TRUE)
  plot_type  <- if (slot(x, "PlotType") == "Auto") {
    envir$plot.type
  } else {
    "scatter2"
  }

  args <- list(plot_data,
               param_choices = x,
               x_lab = labels$X,
               y_lab = labels$Y,
               color_lab = labels$ColorBy,
               shape_lab = labels$ShapeBy,
               size_lab = labels$SizeBy,
               title = labels$title,
               is_subsetted = is_subsetted,
               is_downsampled = is_downsampled)

  plot_cmds <- switch(plot_type,
                      square=do.call(iSEE:::.square_plot, args),
                      violin=do.call(iSEE:::.violin_plot, args),
                      violin_horizontal=do.call(iSEE:::.violin_plot, 
                                                c(args, list(horizontal=TRUE))),
                      scatter=do.call(iSEE:::.scatter_plot, args),
                      scatter2 = do.call(.scatter_plot2, args)
  )

  # Line layer — only when explicitly requested
  if (add_lines) {
    N <- length(plot_cmds)
    color_set <- !is.null(plot_data$ColorBy)
    aes_line  <- iSEE:::.buildAes(
      color = color_set, group = TRUE,
      alt   = c(color = iSEE:::.set_colorby_when_none(x))
    )
    # Map GroupBy --> id (the altExp feature identifier)
    aes_line <- gsub("GroupBy", "id", aes_line, fixed = TRUE)
    plot_cmds[[N]] <- paste(plot_cmds[[N]], "+")
    plot_cmds[["geom_line"]] <- sprintf(
      "geom_line(data = plot.data, mapping = %s)", 
      aes_line
    )
  }

  # Adding a faceting command, if applicable.
  facet_cmd <- iSEE:::.addFacets(x)
  if (length(facet_cmd)) {
    N <- length(plot_cmds)
    plot_cmds[[N]] <- paste(plot_cmds[[N]], "+")
    plot_cmds <- c(plot_cmds, facet_cmd)
  }

  plot_cmds <- iSEE:::.addCustomLabelsCommands(
    x, 
    commands = plot_cmds, 
    plot_type = plot_type)

  if (plot_type=="scatter") {
    plot_cmds <- iSEE:::.addLabelCentersCommands(x, commands=plot_cmds)
  }

  # Adding self-brushing boxes, if they exist.
  plot_cmds <- iSEE:::.addMultiSelectionPlotCommands(
    x,
    flip = (plot_type == "violin_horizontal"),
    envir = envir, 
    commands = plot_cmds)

  list(plot=iSEE:::.textEval(plot_cmds, envir), commands=plot_cmds)
})


############################################################
# Panel tour
############################################################

#' @importFrom iSEE .getEncodedName .getPanelColor .addTourStep .dataParamBoxOpen
#'   .definePanelTour
setMethod(".definePanelTour", "MainExpPlot", function(x) {
  
  collated <- rbind(
    c(
      element = paste0("#", .getEncodedName(x)),
      intro = sprintf(
        "The <font color=\"%s\">MainExp plot</font> panel displays assay values
from the <em>main experiment</em> of a <code>SummarizedExperiment</code> or
<code>SingleCellExperiment</code> object.  Each point represents one sample
(column), and the y-axis shows the assay value of the selected feature (row).
<br><br>
Unlike the standard Feature assay plot, this panel can optionally connect
measurements across samples using <strong>Scatter + lines</strong>, which is
useful for ordered sample layouts such as time courses, pseudotime trajectories,
or experiments where samples have a meaningful ordering.",
        .getPanelColor(x)
      )
    ),
    c(
      element = paste0("#", .getEncodedName(x)),
      intro = "The selected feature is taken from the main experiment.  Its
assay values are plotted directly on the y-axis.  Additional visual mappings
such as colour, shape, size, and faceting can be used to highlight sample-level
metadata or other biological annotations."
    ),
    c(
      element = paste0("#", .getEncodedName(x)),
      intro = "When <strong>Scatter + lines</strong> is selected, points
belonging to the same feature are connected with a line.  This is most useful
when the x-axis represents an ordered variable, for example sampling time,
developmental progression, treatment dose, or an ordered sample identifier."
    ),
    .addTourStep(
      x,
      .dataParamBoxOpen,
      "The <i>Data parameters</i> box contains the controls for choosing the
feature, assay, x-axis, and plotting mode.<br><br>
<strong>Action:</strong> click on this box to expand it and explore the
available options."
    )
  )
  
  parent_tour <- callNextMethod()
  
  parent_tour <- parent_tour[
    !grepl("Feature assay plot", parent_tour$intro),
  ]
  
  rbind(
    data.frame(
      element = collated[, 1],
      intro = collated[, 2],
      stringsAsFactors = FALSE
    ),
    parent_tour
  )
})