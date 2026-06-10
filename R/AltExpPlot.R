############################################################
# AltExpPlot
############################################################

# Start from native FeatureAssayPlot
# Lineplot is not available yet
# Everything works
# Has to be extended to include lineplot
#' @export

setClass("AltExpPlot", contains="FeatureAssayPlot",
         slots=c(AltExp="character",
                 AltAssay="character",
                 MapColumn="character",
                 PlotType="character"))

setValidity2("AltExpPlot", function(object) {
  msg <- character(0)

  #msg <- .validStringError(msg, object, "MapColumn")
  msg <- .validStringError(msg, object, "PlotType")

  if (length(msg)) {
    return(msg)
  }
  TRUE
})

setMethod("initialize", "AltExpPlot", function(.Object, ...) {
  args <- list(...)
  args <- .emptyDefault(args, "AltExp", NA_character_)
  args <- .emptyDefault(args, "AltAssay", NA_character_)
  args <- .emptyDefault(args, "MapColumn", NA_character_)
  args <- .emptyDefault(args, "PlotType", "Auto")
  args <- .emptyDefault(
    args,
    iSEE:::.shapeByField,
    iSEE:::.shapeByColDataTitle
  )

  args <- .emptyDefault(
    args,
    iSEE:::.shapeByColData,
    "id"
  )
  # args <- .emptyDefault(
  #   args,
  #   iSEE:::.visualParamChoice,
  #   c(
  #     iSEE:::.visualParamChoiceColorTitle,
  #     iSEE:::.visualParamChoiceShapeTitle
  #   )
  # )

  do.call(callNextMethod, c(list(.Object), args))
})

setMethod(".refineParameters", "AltExpPlot", function(x, se) {
  x <- callNextMethod()
  all_altexps <- altExpNames(se)
  all_assays <- unlist(lapply(altExpNames(se), function(nm) {
    assayNames(altExp(se, nm))
  }))
  if (length(all_assays)==0L) {
    warning(sprintf("no valid 'assays' for plotting '%s'", class(x)[1]))
    return(NULL)
  }
  x <- .replaceMissingWithFirst(x, "AltExp", all_altexps)
  x <- .replaceMissingWithFirst(x, "AltAssay", all_assays)

  if (is.null(slot(x, iSEE:::.shapeByColData)) ||
      !nzchar(slot(x, iSEE:::.shapeByColData)) ||
      slot(x, iSEE:::.shapeByColData) %in% colnames(colData(se))) {

    # Only override if user hasn’t explicitly chosen something else
    slot(x, iSEE:::.shapeByField) <- iSEE:::.shapeByColDataTitle
    slot(x, iSEE:::.shapeByColData) <- "id"
  }


  if (is.null(slot(x, iSEE:::.colorByColData)) ||
      !nzchar(slot(x, iSEE:::.colorByColData)) ||
      slot(x, iSEE:::.colorByColData) %in% colnames(colData(se))) {

    slot(x, iSEE:::.colorByField) <- iSEE:::.colorByColDataTitle
    slot(x, iSEE:::.colorByColData) <- "id"
  }

  x
})

#' @export
AltExpPlot <- function(...) {
  new("AltExpPlot", ...)
}

#' @export
setMethod(".fullName", "AltExpPlot", function(x) "AltExp plot")

#' @export
setMethod(".panelColor", "AltExpPlot", function(x) "#AA5500")

#' @export
setMethod(".defineDataInterface", "AltExpPlot", function(x, se, select_info) {
  panel_name <- .getEncodedName(x)
  .input_FUN <- function(field) { paste0(panel_name, "_", field) }

  all_altexps <- altExpNames(se)
  all_assays <- unlist(lapply(altExpNames(se), function(nm) {
    assayNames(altExp(se, nm))
  }))
  all_rowdata_colnames <- colnames(rowData(se))
  all_rowdata_colnames_altExp_sel <- colnames(rowData(altExp(se, slot(x,"AltExp"))))
  col_classes <- sapply(rowData(se), class)
  valid_mapcols <- all_rowdata_colnames[(col_classes %in% c("character", "factor")) & (all_rowdata_colnames %in%all_rowdata_colnames_altExp_sel)]

  tab_by_row <- select_info$single$feature

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
      .input_FUN("AltExp"), label="AltExp:", choices=all_altexps,
      selected=iSEE:::.choose_link(slot(x, "AltExp"), all_altexps)),
    selectInput(
      .input_FUN("AltAssay"), label="Alt assay:", choices=all_assays,
      selected=iSEE:::.choose_link(slot(x, "AltAssay"), all_assays)),
    selectInput(
      .input_FUN("MapColumn"), label="Map column:", choices=valid_mapcols,
      selected=iSEE:::.choose_link(slot(x, "MapColumn"), valid_mapcols)),
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
setMethod(".createObservers", "AltExpPlot", function(x, se, input, session, pObjects, rObjects) {
  callNextMethod()

  plot_name <- .getEncodedName(x)

  .createUnprotectedParameterObservers(
    plot_name,
    fields=c("AltExp", "AltAssay", "MapColumn", "PlotType"),
    input=input, pObjects=pObjects, rObjects=rObjects)

  invisible(NULL)
})

#' @export
setMethod(".generateDotPlotData", "AltExpPlot", function(x, envir) {
  data_cmds <- list()

  gene_selected_y <- slot(x, iSEE:::.featAssayYAxisFeatName)

  data_cmds[["y"]] <- c(
    sprintf("ae <- altExp(se, %s);", deparse(slot(x, "AltExp"))),
    sprintf("pg <- rowData(se)[[%s]][which(rownames(se) %%in%% %s)];",
            deparse(slot(x, "MapColumn")), deparse(gene_selected_y)),
    sprintf("alt <- ae[rowData(ae)[[%s]] %%in%% pg,];",
            deparse(slot(x, "MapColumn"))),
    sprintf(
      "plot.data <- data.frame(assay(alt, %s)) |> tibble::rownames_to_column('id') |> tidyr::pivot_longer(names_to='sample',values_to='Y',-'id');",
      deparse(slot(x, "AltAssay")))
    )

  x_choice <- slot(x, iSEE:::.featAssayXAxis)

  if (x_choice == iSEE:::.featAssayXAxisColDataTitle) { # colData column selected
    data_cmds[["x"]] <-  sprintf("plot.data$X <- rep(colData(se)[,%s], length(unique(plot.data$id)));", deparse(slot(x, iSEE:::.featAssayXAxisColData)))
  } else { # no x axis variable specified: show single violin
    data_cmds[["x"]] <- "plot.data$X <- as.factor('')"
  }

  data_cmds <- unlist(data_cmds)
  .textEval(data_cmds, envir)

  list(commands=data_cmds, labels=list(title= gene_selected_y, X="", Y=slot(x, "AltAssay")))
})


############################################################
# To control visual parameters
############################################################

# Copied from DotPlot and modified on one line:  .add_extra_aesthetic_columns_AltExp

#' @export
setMethod(".generateOutput", "AltExpPlot", function(x, se, all_memory, all_contents) {
  # Initialize an environment storing information for generating ggplot commands
  plot_env <- new.env()
  plot_env$se <- se
  plot_env$colormap <- iSEE:::.get_colormap(se)

  all_cmds <- list()
  all_labels <- list()

  # Doing this first so that .generateDotPlotData can respond to the selection.
  all_cmds$select <- iSEE:::.processMultiSelections(x, all_memory, all_contents, plot_env)

  xy_out <- iSEE:::.generateDotPlotData(x, plot_env)
  all_cmds$xy <- xy_out$commands
  all_labels <- c(all_labels, xy_out$labels)

  # Line below is the only change compared to the DotPlot method to ensure that the aesthetics dimensions are correct
  extra_out <- .add_extra_aesthetic_columns_AltExp(x, plot_env)
  all_cmds <- c(all_cmds, extra_out$commands)
  all_labels <- c(all_labels, extra_out$labels)

  select_out2 <- iSEE:::.add_selectby_column(x, plot_env)
  all_cmds <- c(all_cmds, select_out2)

  # We need to set up the plot type before downsampling,
  # to ensure the X/Y jitter is correctly computed.
  all_cmds$setup <-  iSEE:::.choose_plot_type(plot_env)

  # Also collect the plot coordinates before downsampling.
  panel_data <- plot_env$plot.data

  # Non-data-related fiddling to affect the visual display.
  # First, scrambling the plot.data to avoid biases.
  scramble_cmds <- c(
    "# Avoid visual biases from default ordering by shuffling the points",
    sprintf("set.seed(%i);", nrow(panel_data)), # Using a deterministically different seed to keep things exciting.
    "plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];"
  )
  iSEE:::.textEval(scramble_cmds, plot_env)
  all_cmds$shuffle <- scramble_cmds

  # Next, reordering by priority (this is stable so any ordering due to the
  # shuffling above is still preserved within each priority level).
  priority_out <- iSEE:::.prioritizeDotPlotData(x, plot_env)
  rescaled_res <- FALSE
  if (has_priority <- !is.null(priority_out)) {
    order_cmds <- "plot.data <- plot.data[order(.priority),,drop=FALSE];"
    iSEE:::.textEval(order_cmds, plot_env)
    all_cmds$priority <- c(priority_out$commands, order_cmds)
    rescaled_res <- priority_out$rescaled
  }
  # Finally, the big kahuna of downsampling.
  all_cmds$downsample <- iSEE:::.downsample_points(x, plot_env, priority=has_priority, rescaled=rescaled_res)

  plot_out <- iSEE:::.generateDotPlot(x, all_labels, plot_env)
  all_cmds$plot <- plot_out$commands

  list(commands=all_cmds, contents=panel_data, plot=plot_out$plot, varname="plot.data")

})

### Add function to override non exported functions
# Copied from .add_extra_aesthetic_columns in outputs_plot.R
# Include functions to override inside to adjust for dimensions
# .addDotPlotDataColor, adjusted to include feature id from altExp to color
# .addDotPlotDataShape, adjusted to include feature id from altExp for shape
# .addDotPlotDataSize
# .addDotPlotDataFacets


.add_extra_aesthetic_columns_AltExp <- function (x, envir)
{
  # We first include all .addDotPlotData functions that we have to override
  .addDotPlotDataColor <- function(x, envir) {
    color_choice <- slot(x, iSEE:::.colorByField)

    if (color_choice == iSEE:::.colorByColDataTitle) {
      covariate_name <- slot(x, iSEE:::.colorByColData)
      label <- covariate_name

      if (covariate_name != "id") {
        # Normal case: column exists in colData
        cmds <- sprintf(
          "plot.data$ColorBy <- rep(colData(se)[, %s], dplyr::n_distinct(plot.data$id));",
          deparse(covariate_name)
        )
      } else {
        # Special case: id is generated in plot.data
        cmds <- "plot.data$ColorBy <- plot.data$id"
      }

    } else if (color_choice == iSEE:::.colorByFeatNameTitle) {
      # Set the color to the selected gene
      chosen_gene <- slot(x, iSEE:::.colorByFeatName)
      assay_choice <- slot(x, iSEE:::.colorByFeatNameAssay)
      label <- sprintf("%s\n(%s)", chosen_gene, assay_choice)
      cmds <- sprintf(
        "plot.data$ColorBy <- rep(assay(se, %s)[%s, ], dplyr::n_distinct(plot.data$id));",
        deparse(assay_choice), deparse(chosen_gene)
      )

    } else if (color_choice == iSEE:::.colorBySampNameTitle) {
      chosen_sample <- slot(x, iSEE:::.colorBySampName)
      label <- chosen_sample
      cmds <- sprintf(
        "plot.data$ColorBy <- logical(nrow(plot.data));\nplot.data[%s, 'ColorBy'] <- TRUE;",
        deparse(chosen_sample)
      )

    } else if (color_choice == iSEE:::.colorByColSelectionsTitle) {
      label <- "Column selection"
      if (exists("col_selected", envir=envir, inherits=FALSE)) {
        target <- "col_selected"
      } else {
        target <- "list()"
      }
      cmds <- sprintf(
        "plot.data$ColorBy <- iSEE::multiSelectionToFactor(%s, colnames(se));",
        target
      )

    } else {
      return(NULL)
    }

    iSEE:::.textEval(cmds, envir)

    list(commands=cmds, labels=list(ColorBy=label))
  }

  .addDotPlotDataShape <- function(x, envir) {
    shape_choice <- slot(x, iSEE:::.shapeByField)

    if (shape_choice == iSEE:::.shapeByColDataTitle) {
      covariate_name <- slot(x, iSEE:::.shapeByColData)

      if (covariate_name == "id") {
        label <- "id"
        cmds <- "plot.data$ShapeBy <- plot.data$id"

      } else if (covariate_name %in% colnames(colData(envir$se))) {
        # normal iSEE behavior
        label <- covariate_name
        cmds <- sprintf(
          "plot.data$ShapeBy <- rep(colData(se)[, %s], dplyr::n_distinct(plot.data$id));",
          deparse(covariate_name)
        )

      } else {
        # fallback safety
        label <- "None"
        cmds <- "plot.data$ShapeBy <- NULL"
      }

    } else {
      # default fallback → use id
      label <- "id"
      cmds <- "plot.data$ShapeBy <- plot.data$id"
    }

    iSEE:::.textEval(cmds, envir)

    list(commands=cmds, labels=list(ShapeBy=label))
  }

  .addDotPlotDataSize <- function(x, envir) {
    size_choice <- slot(x, iSEE:::.sizeByField)

    if (size_choice == iSEE:::.sizeByColDataTitle) {
      covariate_name <- slot(x, iSEE:::.sizeByColData)
      label <- covariate_name
      cmds <- sprintf("plot.data$SizeBy <- rep(colData(se)[, %s], dplyr::n_distinct(plot.data$id));", deparse(covariate_name))

    } else {
      return(NULL)
    }

    iSEE:::.textEval(cmds, envir)

    list(commands=cmds, labels=list(SizeBy=label))
  }

  .addDotPlotDataFacets <-  function(x, envir) {
    facet_cmds <- NULL
    labels <- list()

    params <- list(
      list(iSEE:::.facetRow, "FacetRow", iSEE:::.facetRowByColData),
      list(iSEE:::.facetColumn, "FacetColumn", iSEE:::.facetColumnByColData)
    )

    for (f in seq_len(2)) {
      current <- params[[f]]
      param_field <- current[[1]]
      pd_field <- current[[2]]
      facet_mode <- slot(x, param_field)

      if (facet_mode == iSEE:::.facetByColDataTitle) {
        facet_data <- x[[current[[3]]]]
        facet_cmds[pd_field] <- sprintf("plot.data$%s <- rep(colData(se)[, %s], dplyr::n_distinct(plot.data$id));", pd_field, deparse(facet_data))
        labels[[pd_field]] <- facet_data

      } else if (facet_mode == iSEE:::.facetByColSelectionsTitle) {
      if (exists("col_selected", envir=envir, inherits=FALSE)) {
       target <- "col_selected"
      } else {
       target <- "list()"
      }
      facet_cmds[pd_field] <- sprintf("plot.data$%s <- iSEE::multiSelectionToFactor(%s, colnames(se));", pd_field, target)
      labels[[pd_field]] <- "Column selection"
      }
    }

    iSEE:::.textEval(facet_cmds, envir)

    list(commands=facet_cmds, labels=labels)
  }

  collected <- list()
  labels <- list()
  collected$coerce <- iSEE:::.coerce_dataframe_columns(envir,
                                                              fields = c("X", "Y"),
                                                              df = "plot.data",
                                                              max_levels = iSEE:::.get_factor_maxlevels())
  # Add commands adding optional columns to plot.data
  out_color <- .addDotPlotDataColor(x, envir)
  collected$color <- out_color$commands
  labels <- c(labels, out_color$labels)
  if (!is.null(envir$plot.data$ColorBy)) {
    collected$color <- c(collected$color, iSEE:::.coerce_dataframe_columns(envir,
                                                                    fields = "ColorBy", df = "plot.data", max_levels = iSEE:::.get_color_maxlevels()))
  }
  out_shape <- .addDotPlotDataShape(x, envir)
  collected$shape <- out_shape$commands
  labels <- c(labels, out_shape$labels)

  out_size <- .addDotPlotDataSize(x, envir)
  collected$size <- out_size$commands
  labels <- c(labels, out_size$labels)

  out_facets <- .addDotPlotDataFacets(x, envir)
  collected$facets <- out_facets$commands
  labels <- c(labels, out_facets$labels)

  list(commands = collected, labels = labels)
}

### Redefine .defineVisual interfaces
### Reduce color options as compared to DotPlot class
setMethod(".defineVisualColorInterface", "AltExpPlot", function(x, se, select_info) {
  all_assays <- iSEE:::.getCachedCommonInfo(se, "AltExpPlot")$valid.assay.names

  plot_name <- iSEE:::.getEncodedName(x)
  colorby_field <- paste0(plot_name, "_", iSEE:::.colorByField)

  colorby <- iSEE:::.getDotPlotColorConstants(x)
  mydim_single <- iSEE:::.singleSelectionDimension(x)
  otherdim_single <- setdiff(c("feature", "sample"), mydim_single)
  mydim_choices <- select_info[[mydim_single]]
  otherdim_choices <- select_info[[otherdim_single]]

  #color_choices <-iSEE:::.defineDotPlotColorChoices(x, se)


  # Actually creating the UI.
  # custom color_choices deviating from DotPlot class
  covariates <- iSEE:::.allowableColorByDataChoices(x, se)
  color_choices <- c(iSEE:::.colorByNothingTitle)
  if (length(covariates)) {
    color_choices <- c(color_choices, iSEE:::.colorByColDataTitle)
  }

  tagList(
    hr(),
    iSEE:::.radioButtons.iSEE(x, iSEE:::.colorByField,
                       label="Color by:",
                       inline=TRUE,
                       choices=color_choices,
                       selected=slot(x, iSEE:::.colorByField)
    ),
    iSEE:::.conditionalOnRadio(
      colorby_field, iSEE:::.colorByNothingTitle,
      colourpicker::colourInput(
        paste0(plot_name, "_", iSEE:::.colorByDefaultColor), label=NULL,
        value=slot(x, iSEE:::.colorByDefaultColor))
    ),
    iSEE:::.conditionalOnRadio(
      colorby_field, colorby$metadata$title,
      selectInput(
        paste0(plot_name, "_", colorby$metadata$field), label=NULL,
        choices=c("id",iSEE:::.allowableColorByDataChoices(x, se)),
        selected=x[[colorby$metadata$field]])
    ),
    iSEE:::.conditionalOnRadio(colorby_field, colorby$name$title,
                        selectizeInput(paste0(plot_name, "_", colorby$name$field),
                                       label=NULL, selected=NULL, choices=NULL, multiple=FALSE),
                        selectInput(
                          paste0(plot_name, "_", colorby$name$table), label=NULL, choices=mydim_choices,
                          selected=iSEE:::.choose_link(x[[colorby$name$table]], mydim_choices)),
                        colourpicker::colourInput(paste0(plot_name, "_", colorby$name$color), label=NULL,
                                    value=x[[colorby$name$color]]),
                        checkboxInput(
                          paste0(plot_name, "_", colorby$name$dynamic),
                          label=sprintf("Use dynamic %s selection", mydim_single),
                          value=x[[colorby$name$dynamic]]),
    ),
    iSEE:::.conditionalOnRadio(colorby_field, colorby$assay$title,
                        selectizeInput(paste0(plot_name, "_", colorby$assay$field),
                                       label=NULL, choices=NULL, selected=NULL, multiple=FALSE),
                        selectInput(
                          paste0(plot_name, "_", colorby$assay$assay), label=NULL,
                          choices=all_assays, selected=x[[colorby$assay$assay]]),
                        selectInput(
                          paste0(plot_name, "_", colorby$assay$table), label=NULL, choices=otherdim_choices,
                          selected=iSEE:::.choose_link(x[[colorby$assay$table]], otherdim_choices)),
                        checkboxInput(
                          paste0(plot_name, "_", colorby$assay$dynamic),
                          label=sprintf("Use dynamic %s selection", otherdim_single),
                          value=x[[colorby$assay$dynamic]])
    ),
    iSEE:::.sliderInput.iSEE(x, iSEE:::.selectTransAlpha,
                      label="Unselected point opacity:",
                      min=0, max=1, value=slot(x, iSEE:::.selectTransAlpha)
    )
  )
})


#Custom Shape deviating from DotPlot
#' @export
setMethod(".defineVisualShapeInterface", "AltExpPlot", function (x, se)
  {
    discrete_covariates <- c("id",iSEE:::.getDiscreteMetadataChoices(x, se))
    if (length(discrete_covariates)) {
      plot_name <- iSEE:::.getEncodedName(x)
      shapeby_field <- paste0(plot_name, "_", iSEE:::.shapeByField)
      shapeby <- iSEE:::.getDotPlotShapeConstants(x)


      tagList(hr(), iSEE:::.radioButtons.iSEE(x, iSEE:::.shapeByField, label = "Shape by:",
                                       inline = TRUE, choices = c(iSEE:::.shapeByNothingTitle,
                                                                  shapeby$metadata$title), selected = slot(x, iSEE:::.shapeByField)),
              iSEE:::.conditionalOnRadio(shapeby_field, shapeby$metadata$title,
                                  selectInput(paste0(plot_name, "_", shapeby$metadata$field),
                                              label = NULL, choices = discrete_covariates,
                                              selected = x[[shapeby$metadata$field]])))
    }
    else {
      NULL
    }
  })

#' @export
setMethod(".generateDotPlot", "AltExpPlot", function(x, labels, envir) {
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
