############################################################
# MainExpPlot — modified from FeatureAssayPlot
############################################################

# Start from native FeatureAssayPlot
# adds Lineplot
#' @export
setClass("MainExpPlot", contains="FeatureAssayPlot",
         slots=c(PlotType="character"))

setValidity2("MainExpPlot", function(object) {
  msg <- character(0)

  msg <- .validStringError(msg, object, "PlotType")

  if (length(msg)) {
    return(msg)
  }
  TRUE
})

setMethod("initialize", "MainExpPlot", function(.Object, ...) {
  args <- list(...)
  args <- .emptyDefault(args, "PlotType", "Auto")
  do.call(callNextMethod, c(list(.Object), args))
})

#' @export
MainExpPlot <- function(...) {
  new("MainExpPlot", ...)
}

#' @export
setMethod(".fullName", "MainExpPlot", function(x) "MainExp plot")

#' @export
setMethod(".panelColor", "MainExpPlot", function(x) "#7BB854")

#' @export
setMethod(".defineDataInterface", "MainExpPlot", function(x, se, select_info) {
  panel_name <- .getEncodedName(x)
  .input_FUN <- function(field) { paste0(panel_name, "_", field) }

  tab_by_row <- select_info$single$feature

  all_assays <- assayNames(se)

  column_covariates <- .getCachedCommonInfo(se, "ColumnDotPlot")$valid.colData.names
  xaxis_choices <- c(iSEE:::.featAssayXAxisNothingTitle)
  if (length(column_covariates)) { # As it is possible for this plot to be _feasible_ but for no column data to exist.
    xaxis_choices <- c(xaxis_choices, iSEE:::.featAssayXAxisColDataTitle)
  }
  #xaxis_choices <- c(xaxis_choices, iSEE:::.featAssayXAxisFeatNameTitle, iSEE:::.featAssayXAxisSelectionsTitle)

  list(
    .selectizeInput.iSEE(
      x, iSEE:::.featAssayYAxisFeatName,
      label="Feature in main assay:",
      choices=NULL,
      selected=NULL,
      multiple=FALSE),
    selectInput(
      .input_FUN(iSEE:::.featAssayYAxisRowTable), label=NULL, choices=tab_by_row,
      selected=iSEE:::.choose_link(slot(x, iSEE:::.featAssayYAxisRowTable), tab_by_row)),
    checkboxInput(
      .input_FUN(iSEE:::.featAssayYAxisFeatDynamic),
      label="Use dynamic feature selection for the y-axis",
      value=slot(x, iSEE:::.featAssayYAxisFeatDynamic)),
    selectInput(
      .input_FUN("Assay"), label="Assay:", choices=all_assays,
      selected=iSEE:::.choose_link(slot(x, "Assay"), all_assays)),
    selectInput(
      paste0(panel_name, "_", "PlotType"), label="Plot type:",
      choices=c("Auto","Scatter", "Scatter + lines"), selected=slot(x,"PlotType")),
    .radioButtons.iSEE(x, iSEE:::.featAssayXAxis,
                       label="X-axis:",
                       inline=TRUE,
                       choices=xaxis_choices,
                       selected=slot(x, iSEE:::.featAssayXAxis)),

    .conditionalOnRadio(.input_FUN(iSEE:::.featAssayXAxis),
                        iSEE:::.featAssayXAxisColDataTitle,
                        selectInput(.input_FUN(iSEE:::.featAssayXAxisColData),
                                    label="X-axis column data:",
                                    choices=column_covariates, selected=slot(x, iSEE:::.featAssayXAxisColData))),

    .conditionalOnRadio(.input_FUN(iSEE:::.featAssayXAxis),
                        iSEE:::.featAssayXAxisFeatNameTitle,
                        selectizeInput(.input_FUN(iSEE:::.featAssayXAxisFeatName),
                                       label="X-axis feature:", choices=NULL, selected=NULL, multiple=FALSE),
                        selectInput(.input_FUN(iSEE:::.featAssayXAxisRowTable), label=NULL,
                                    choices=tab_by_row, selected=slot(x, iSEE:::.featAssayXAxisRowTable)),
                        checkboxInput(.input_FUN(iSEE:::.featAssayXAxisFeatDynamic),
                                      label="Use dynamic feature selection for the x-axis",
                                      value=slot(x, iSEE:::.featAssayXAxisFeatDynamic))
    )
  )
})

#' @export
setMethod(".createObservers", "MainExpPlot", function(x, se, input, session, pObjects, rObjects) {
  callNextMethod()

  plot_name <- .getEncodedName(x)

  .createUnprotectedParameterObservers(
    plot_name,
    fields=c("PlotType"),
    input=input, pObjects=pObjects, rObjects=rObjects)

  invisible(NULL)
})


#' @export
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

  if (x_choice == iSEE:::.featAssayXAxisColDataTitle) { # colData column selected
    x_lab <- slot(x, iSEE:::.featAssayXAxisColData)
    plot_title <- paste(plot_title, "vs", x_lab)
    data_cmds[["x"]] <- sprintf("plot.data$X <- colData(se)[, %s];", deparse(x_lab))

  } else if (x_choice == iSEE:::.featAssayXAxisFeatNameTitle) { # gene selected
    gene_selected_x <- slot(x, iSEE:::.featAssayXAxisFeatName)
    plot_title <- paste(plot_title, "vs", gene_selected_x)
    x_lab <- sprintf("%s (%s)", gene_selected_x, assay_choice)
    data_cmds[["x"]] <- sprintf(
      "plot.data$X <- assay(se, %s)[%s, ];",
      deparse(assay_choice), deparse(gene_selected_x)
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

  list(commands=data_cmds, labels=list(title=plot_title, X=x_lab, Y=y_lab))
})

setMethod(".generateDotPlot", "MainExpPlot", function(x, labels, envir) {
  plot_data <- envir$plot.data

  is_subsetted <- exists("plot.data.all", envir=envir, inherits=FALSE)
  is_downsampled <- exists("plot.data.pre", envir=envir, inherits=FALSE)
  if(slot(x,"PlotType") == "Auto") plot_type <- envir$plot.type else
    plot_type = "scatter2"

  args <- list(plot_data,
               param_choices=x,
               x_lab=labels$X,
               y_lab=labels$Y,
               color_lab=labels$ColorBy,
               shape_lab=labels$ShapeBy,
               size_lab=labels$SizeBy,
               title=labels$title,
               is_subsetted=is_subsetted,
               is_downsampled=is_downsampled)

  plot_cmds <- switch(plot_type,
                      square=do.call(iSEE:::.square_plot, args),
                      violin=do.call(iSEE:::.violin_plot, args),
                      violin_horizontal=do.call(iSEE:::.violin_plot, c(args, list(horizontal=TRUE))),
                      scatter=do.call(iSEE:::.scatter_plot, args),
                      scatter2 = do.call(.scatter_plot2, args)
  )

  if (slot(x, iSEE:::.shapeByField) != iSEE:::.shapeByNothingTitle)
    if (dplyr::n_distinct(plot_data$ShapeBy) < 26){
      N <- length(plot_cmds)
      plot_cmds[[N]] <- paste(plot_cmds[[N]], "+")
      plot_cmds[[N+1]] <- "scale_shape_manual(values = seq_along(plot.data$ShapeBy))"
    }


  if(grepl("lines", slot(x,"PlotType"))) {
    N <- length(plot_cmds)
    plot_cmds[[N]] <- paste(plot_cmds[[N]], "+")
    color_set<- !is.null(plot_data$ColorBy)
    aes_line <- iSEE:::.buildAes(color = color_set, group = TRUE, alt = c(color = iSEE:::.set_colorby_when_none(x)))
    aes_line <- gsub(pattern="GroupBy", replacement = "id", x =aes_line)
    plot_cmds[["geom_line"]] <- sprintf("geom_line(data = plot.data, mapping = %s)", aes_line)
  }

  # Adding a faceting command, if applicable.
  facet_cmd <- iSEE:::.addFacets(x)
  if (length(facet_cmd)) {
    N <- length(plot_cmds)
    plot_cmds[[N]] <- paste(plot_cmds[[N]], "+")
    plot_cmds <- c(plot_cmds, facet_cmd)
  }

  plot_cmds <- iSEE:::.addCustomLabelsCommands(x, commands=plot_cmds, plot_type=plot_type)

  if (plot_type=="scatter") {
    plot_cmds <- iSEE:::.addLabelCentersCommands(x, commands=plot_cmds)
  }

  # Adding self-brushing boxes, if they exist.
  plot_cmds <- iSEE:::.addMultiSelectionPlotCommands(x,
                                                     flip=(plot_type == "violin_horizontal"),
                                                     envir=envir, commands=plot_cmds)

  list(plot=iSEE:::.textEval(plot_cmds, envir), commands=plot_cmds)
})

.scatter_plot2 <- function(plot_data, param_choices, x_lab, y_lab, color_lab,
                           shape_lab, size_lab, title, by_row = FALSE, is_subsetted = FALSE,
                           is_downsampled = FALSE) {

  plot_cmds <- list()
  if (dplyr::n_distinct(plot_data$id) > 1000)
  {
    commands <- "ggplot() + ggtitle('More than 1000 features selected of altExp, is mapping column correct?')"
    return(commands=commands)
  }


  plot_cmds[["ggplot"]] <- "dot.plot <- ggplot() +"
  color_set <- !is.null(plot_data$ColorBy)
  shape_set <- slot(param_choices, iSEE:::.shapeByField) != iSEE:::.shapeByNothingTitle
  size_set <- slot(param_choices, iSEE:::.sizeByField) != iSEE:::.sizeByNothingTitle
  new_aes <- iSEE:::.buildAes(color = color_set, shape = shape_set,
                              size = size_set, alt = c(color = iSEE:::.set_colorby_when_none(param_choices)))
  plot_cmds[["points"]] <- iSEE:::.create_points(param_choices, !is.null(plot_data$SelectBy),
                                                 new_aes, color_set, size_set)
  color_scale_cmd <- iSEE:::.colorDotPlot(param_choices, plot_data$ColorBy)
  guides_cmd <- iSEE:::.create_guides_command(param_choices, plot_data$ColorBy)
  plot_cmds[["labs"]] <- iSEE:::.buildLabs(x = x_lab, y = y_lab, color = color_lab,
                                           shape = shape_lab, size = size_lab, title = title)
  plot_cmds[["scale_color"]] <- color_scale_cmd
  plot_cmds[["guides"]] <- guides_cmd
  plot_cmds[["theme_base"]] <- "theme_bw() +"
  font_size <- slot(param_choices, iSEE:::.plotFontSize)

  if(is.numeric(plot_data$X))
    plot_cmds[["theme_custom"]] <- sprintf("theme(legend.position='%s',
                                           legend.box='vertical',
                                           legend.text=element_text(size=%s),
                                           legend.title=element_text(size=%s),
                                           axis.text=element_text(size=%s),
                                           axis.title=element_text(size=%s),
                                           title=element_text(size=%s))",
                                           tolower(slot(param_choices, iSEE:::.plotLegendPosition)),
                                           font_size *
                                             iSEE:::.plotFontSizeLegendTextDefault,
                                           font_size * iSEE:::.plotFontSizeLegendTitleDefault,
                                           font_size * iSEE:::.plotFontSizeAxisTextDefault,
                                           font_size * iSEE:::.plotFontSizeAxisTitleDefault,
                                           font_size * iSEE:::.plotFontSizeTitleDefault) else
                                             plot_cmds[["theme_custom"]] <- sprintf("theme(legend.position='%s',
                                           legend.box='vertical',
                                           legend.text=element_text(size=%s),
                                           legend.title=element_text(size=%s),
                                           axis.text.x=element_text(angle=45, size=%s, hjust=1),
                                           axis.title=element_text(size=%s),
                                           title=element_text(size=%s))",
                                                                                    tolower(slot(param_choices, iSEE:::.plotLegendPosition)),
                                                                                    font_size *
                                                                                      iSEE:::.plotFontSizeLegendTextDefault,
                                                                                    font_size * iSEE:::.plotFontSizeLegendTitleDefault,
                                                                                    font_size * iSEE:::.plotFontSizeAxisTextDefault,
                                                                                    font_size * iSEE:::.plotFontSizeAxisTitleDefault,
                                                                                    font_size * iSEE:::.plotFontSizeTitleDefault)

  unlist(plot_cmds)
}
