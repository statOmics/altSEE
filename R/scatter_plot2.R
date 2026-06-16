###########################################################################
# scatter_plot2 — base scatter plot renderer for LinkedFeaturesAssayPlot
###########################################################################

#' Base scatter plot renderer for LinkedFeaturesAssayPlot
#'
#' An internal helper that produces a \pkg{ggplot2} command list for a scatter
#' plot of long-format altExp data.  It is called by the
#' \code{.generateDotPlot} method for
#' \code{\linkS4class{LinkedFeaturesAssayPlot}} whenever \code{PlotType} is
#' \code{"Scatter"} or \code{"Scatter + lines"}.
#'
#' Unlike the standard \pkg{iSEE} scatter plot helpers, this function operates
#' on the \emph{long-format} \code{plot.data} produced by
#' \code{.generateDotPlotData} for \code{LinkedFeaturesAssayPlot}: the data
#' frame contains one row per (feature, sample) combination, with columns
#' \code{altExp_feature_id} (the visualised feature's identifier),
#' \code{sample} (column name in the main \code{SingleCellExperiment}),
#' \code{X}, and \code{Y}.
#'
#' If more than 1 000 distinct features are present — which typically
#' indicates a misconfigured \code{MapColumn} — the function returns a
#' placeholder plot with an explanatory title instead of attempting to render
#' thousands of overlapping lines.
#'
#' @param plot_data \code{data.frame}.  Long-format plot data as produced by
#'   \code{.generateDotPlotData} for \code{LinkedFeaturesAssayPlot}.  Must
#'   contain columns \code{id}, \code{sample}, \code{X}, \code{Y}, and
#'   optionally \code{ColorBy}, \code{ShapeBy}, \code{SizeBy}, \code{SelectBy}.
#' @param param_choices A \code{\linkS4class{LinkedFeaturesAssayPlot}} object
#'   carrying the current panel parameters (used to extract
#'   colour/shape/size/font settings).
#' @param x_lab,y_lab \code{character(1)}.  Axis labels.
#' @param color_lab,shape_lab,size_lab \code{character(1)} or \code{NULL}.
#'   Legend titles for colour, shape, and size aesthetics respectively.
#' @param title \code{character(1)}.  Plot title (typically the selected
#'   main-experiment feature name).
#' @param by_row \code{logical(1)}.  Unused; retained for signature
#'   compatibility with other iSEE scatter helpers.  Default \code{FALSE}.
#' @param is_subsetted \code{logical(1)}.  Whether a column selection subset
#'   is active.  Passed through to \code{\link[iSEE]{.create_points}}.
#'   Default \code{FALSE}.
#' @param is_downsampled \code{logical(1)}.  Whether downsampling has been
#'   applied.  Passed through to \code{\link[iSEE]{.create_points}}.
#'   Default \code{FALSE}.
#'
#' @return A named \code{character} vector of \pkg{ggplot2} commands (one
#'   element per layer / theme block) that, when evaluated with
#'   \code{plot.data} in scope, produces the final \code{ggplot} object stored
#'   in \code{dot.plot}.  If more than 1 000 features are detected, the vector
#'   contains a single \code{ggplot()} call that renders a warning title.
#'
#' @keywords internal
#' @importFrom iSEE .buildAes .buildLabs
#' @importFrom dplyr n_distinct
#' @importFrom ggplot2 ggplot theme_bw theme element_text
.scatter_plot2 <- function(plot_data, param_choices, x_lab, y_lab, color_lab,
                           shape_lab, size_lab, title,
                           by_row = FALSE, is_subsetted = FALSE,
                           is_downsampled = FALSE) {
  plot_cmds <- list()
  
  # Guard: too many features likely means the MapColumn is misconfigured.
  # Return a consistent named list so callers can always do result$commands.
  n_features <- dplyr::n_distinct(plot_data$id)
  if (n_features > 1000) {
    plot_cmds[["ggplot"]] <- paste0(
      "ggplot() + ggtitle('More than 1000 features selected in altExp — ",
      "is the Map column set correctly?')"
    )
    return(list(commands = unlist(plot_cmds)))
  }
  
  color_set <- !is.null(plot_data$ColorBy)
  shape_set <- slot(param_choices, 
                    iSEE:::.shapeByField) != iSEE:::.shapeByNothingTitle
  size_set  <- slot(param_choices, 
                    iSEE:::.sizeByField)  != iSEE:::.sizeByNothingTitle
  
  new_aes <- iSEE:::.buildAes(
    color = color_set, 
    shape = shape_set, 
    size = size_set,
    alt   = c(color = iSEE:::.set_colorby_when_none(param_choices))
  )
  
  plot_cmds[["ggplot"]]      <- "dot.plot <- ggplot() +"
  plot_cmds[["points"]]      <- iSEE:::.create_points(
    param_choices, !is.null(plot_data$SelectBy), new_aes, color_set, size_set
  )
  plot_cmds[["labs"]]        <- iSEE:::.buildLabs(
    x = x_lab, y = y_lab, color = color_lab,
    shape = shape_lab, size = size_lab, title = title
  )
  plot_cmds[["scale_color"]] <- iSEE:::.colorDotPlot(param_choices, plot_data$ColorBy)
  plot_cmds[["guides"]]      <- iSEE:::.create_guides_command(param_choices, plot_data$ColorBy)
  plot_cmds[["theme_base"]]  <- "theme_bw() +"
  
  font_size <- slot(param_choices, iSEE:::.plotFontSize)
  legend_pos <- tolower(slot(param_choices, iSEE:::.plotLegendPosition))
  
  # Shared theme elements
  common_theme <- sprintf(
    "legend.position='%s', legend.box='vertical',
     legend.text=element_text(size=%s),
     legend.title=element_text(size=%s),
     axis.title=element_text(size=%s),
     title=element_text(size=%s)",
    legend_pos,
    font_size * iSEE:::.plotFontSizeLegendTextDefault,
    font_size * iSEE:::.plotFontSizeLegendTitleDefault,
    font_size * iSEE:::.plotFontSizeAxisTitleDefault,
    font_size * iSEE:::.plotFontSizeTitleDefault
  )
  
  plot_cmds[["theme_custom"]] <- if (is.numeric(plot_data$X)) {
    sprintf(
      "theme(%s, axis.text=element_text(size=%s))",
      common_theme,
      font_size * iSEE:::.plotFontSizeAxisTextDefault
    )
  } else {
    sprintf(
      "theme(%s, axis.text.x=element_text(angle=45, size=%s, hjust=1))",
      common_theme,
      font_size * iSEE:::.plotFontSizeAxisTextDefault
    )
  }
  
  unlist(plot_cmds)
}