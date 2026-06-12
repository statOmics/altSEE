############################################################
# AltExpPlot
############################################################

# Extends FeatureAssayPlot from iSEE to support AltExp visualisation.
# Requires a fork of iSEE from statomics/iSEE that exports:
#   .addDotPlotDataColor, .addDotPlotDataShape,
#   .addDotPlotDataSize,  .addDotPlotDataFacets
# Currently also added Remotes: statomics/iSEE in DESCRIPTION

#' AltExpPlot: visualise alternative experiment data in iSEE
#'
#' An S4 class that extends \code{\linkS4class{FeatureAssayPlot}} to plot
#' assay values from an \emph{alternative experiment}
#' (\code{\link[SingleCellExperiment]{altExp}}) stored inside a
#' \code{\linkS4class{SingleCellExperiment}}.  Features in the altExp are
#' linked to features in the main experiment through annotation
#' columns in the altExp (\code{MapColumn}) and mainExp (\code{LookupColumn}), 
#' allowing the visualisation of e.g. protein abundances alongside peptide or 
#' precursor abundances, RNA-seq alongside CITE-seq or ATAC-seq; in the same 
#' iSEE session.
#'
#' @slot AltExp \code{character(1)}.  Name of the alternative experiment to
#'   visualise (must be present in \code{altExpNames(se)}).
#' @slot AltAssay \code{character(1)}.  Name of the assay within \code{AltExp}
#'   whose values are shown on the y-axis.
#' @slot LookupColumn \code{character(1)}.  Name of a \code{character} or
#'   \code{factor} column in \code{rowData} of the \emph{main} experiment.
#'   The value stored for the selected feature in this column is used as the
#'   lookup key — for example \code{"Proteins"} when the main experiment holds
#'   protein-level data identified by a column named "Proteins".
#' @slot MapColumn \code{character(1)}.  Name of a \code{character} or
#'   \code{factor} column in \code{rowData} of the selected \code{AltExp}.
#'   AltExp rows whose value in this column matches the key retrieved via
#'   \code{LookupColumn} are shown — for example \code{"Proteins"} in a
#'   peptide-level altExp that has the protein name stored in a column named 
#'   "Proteins".
#' @slot PlotType \code{character(1)}.  One of \code{"Auto"},
#'   \code{"Scatter"}, or \code{"Scatter + lines"}.  \code{"Auto"} delegates
#'   the choice to iSEE's internal heuristic; \code{"Scatter"} always uses a
#'   dot plot; \code{"Scatter + lines"} adds a \code{geom_line} layer that
#'   connects points belonging to the same altExp feature (useful for ordered
#'   x-axes such as time points or sampleIds that include condition and repeat).
#'
#' @section Inherited slots:
#' All slots from \code{\linkS4class{FeatureAssayPlot}} are inherited,
#' including those controlling x-axis choice, colour, shape, size, faceting,
#' and multi-selection behaviour.
#'
#' @seealso
#' \code{\link{AltExpPlot}} for the constructor.
#' \code{\linkS4class{FeatureAssayPlot}} for the parent class.
#' \code{\link[SingleCellExperiment]{altExp}} for the underlying data accessor.
#'
#' @name AltExpPlot-class
#' @rdname AltExpPlot-class
#' @exportClass AltExpPlot
#' 

setClass("AltExpPlot", contains = "FeatureAssayPlot",
         slots = c(
           AltExp        = "character",
           AltAssay      = "character",
           LookupColumn  = "character",
           MapColumn     = "character",
           PlotType      = "character"
         ))

#' @importFrom iSEE .validStringError
setValidity2("AltExpPlot", function(object) {
  msg <- character(0)
  msg <- .validStringError(msg, object, "PlotType")
  if (length(msg)) return(msg)
  TRUE
})

#' @importFrom iSEE .emptyDefault .shapeByField .shapeByColDataTitle .shapeByColData
setMethod("initialize", "AltExpPlot", function(.Object, ...) {
  args <- list(...)
  args <- .emptyDefault(args, "AltExp",       NA_character_)
  args <- .emptyDefault(args, "AltAssay",     NA_character_)
  args <- .emptyDefault(args, "LookupColumn", NA_character_)
  args <- .emptyDefault(args, "MapColumn",    NA_character_)
  args <- .emptyDefault(args, "PlotType",     "Auto")
  args <- .emptyDefault(args, 
                        iSEE:::.shapeByField,   
                        iSEE:::.shapeByColDataTitle)
  args <- .emptyDefault(args, iSEE:::.shapeByColData, "altExp_feature_id")
  do.call(callNextMethod, c(list(.Object), args))
})


#' @importFrom iSEE .replaceMissingWithFirst .colorByColDataTitle .colorByField
#'   .colorByColData .shapeByField .shapeByColDataTitle .shapeByColData
#' @importFrom SingleCellExperiment altExpNames altExp assayNames
#' @importFrom SummarizedExperiment colData
setMethod(".refineParameters", "AltExpPlot", function(x, se) {
  x <- callNextMethod()
  
  all_altexps <- altExpNames(se)
  
  # Validate selected AltExp first so assay lookup is correct
  x <- .replaceMissingWithFirst(x, "AltExp", all_altexps)
  current_altexp <- slot(x, "AltExp")
  
  all_assays <- assayNames(altExp(se, current_altexp))
  if (length(all_assays) == 0L) {
    warning(sprintf("no valid 'assays' for plotting '%s'", class(x)[1]))
    return(NULL)
  }
  
  # Only replace AltAssay if it is not valid for the current AltExp
  if (!slot(x, "AltAssay") %in% all_assays) {
    x <- .replaceMissingWithFirst(x, "AltAssay", all_assays)
  }
  
  if (is.null(slot(x, iSEE:::.shapeByColData)) ||
      !nzchar(slot(x, iSEE:::.shapeByColData)) ||
      slot(x, iSEE:::.shapeByColData) %in% colnames(colData(se))) {
    slot(x, iSEE:::.shapeByField)   <- iSEE:::.shapeByColDataTitle
    slot(x, iSEE:::.shapeByColData) <- "altExp_feature_id"
  }
  
  if (is.null(slot(x, iSEE:::.colorByColData)) ||
      !nzchar(slot(x, iSEE:::.colorByColData)) ||
      slot(x, iSEE:::.colorByColData) %in% colnames(colData(se))) {
    slot(x, iSEE:::.colorByField)   <- iSEE:::.colorByColDataTitle
    slot(x, iSEE:::.colorByColData) <- "altExp_feature_id"
  }
  
  x
})


#' Construct an AltExpPlot panel
#'
#' Creates an instance of the \code{\linkS4class{AltExpPlot}} class for use
#' as a panel in an iSEE application.
#'
#' @param ... Named arguments corresponding to slots of
#'   \code{\linkS4class{AltExpPlot}} or its parent classes.  Any slot not
#'   supplied receives a sensible default: \code{AltExp} and \code{AltAssay}
#'   default to \code{NA} (resolved at runtime from the
#'   \code{SingleCellExperiment}); \code{MapColumn} defaults to \code{NA};
#'   \code{PlotType} defaults to \code{"Auto"}.
#'
#' @return An \code{\linkS4class{AltExpPlot}} object.
#'
#' @examples
#' \dontrun{
#' # Minimal panel using runtime defaults
#' AltExpPlot()
#'
#' # Pre-configure YAxisFeatureSource, AltAssay, PlotType, XAxis and 
#' # XAxisColumnData
#' # linked via the "Proteins" in rowData columns of MainExp and AltExp, 
#' # with lines connecting points
#' AltExpPlot(YAxisFeatureSource = "RowDataTable1", 
#'            AltAssay = "peptides_norm",
#'            PlotType = "Scatter + lines", 
#'            XAxis = "Column data", 
#'            XAxisColumnData = "sampleId",
#'            LookupColumn = "Proteins",
#'            MapColumn = "Proteins"
#'            )
#'            
#' # Example of iSEE instance with main and altExp
#' data("sceProteinsPeptides")
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
#' @seealso \code{\linkS4class{AltExpPlot}} for a description of all slots.
#' @export

AltExpPlot <- function(...) {new("AltExpPlot", ...)}

setMethod(".fullName", "AltExpPlot", function(x) "AltExp plot")

setMethod(".panelColor", "AltExpPlot", function(x) "#AA5500")

#' @importFrom iSEE .getEncodedName .getCachedCommonInfo .choose_link
#'   .selectizeInput.iSEE .radioButtons.iSEE .conditionalOnRadio
#'   .featAssayYAxisFeatName .featAssayYAxisRowTable .featAssayYAxisFeatDynamic
#'   .featAssayXAxis .featAssayXAxisColDataTitle .featAssayXAxisColData
#'   .featAssayXAxisFeatNameTitle .featAssayXAxisFeatName
#'   .featAssayXAxisRowTable .featAssayXAxisFeatDynamic
#'   .featAssayXAxisNothingTitle .addSpecificTour
#' @importFrom SingleCellExperiment altExpNames altExp assayNames
#' @importFrom SummarizedExperiment rowData
#' @importFrom shiny selectInput checkboxInput selectizeInput
setMethod(".defineDataInterface", "AltExpPlot", function(x, se, select_info) {
  panel_name <- .getEncodedName(x)
  .input_FUN <- function(field) paste0(panel_name, "_", field)
  
  all_altexps  <- altExpNames(se)
  current_ae   <- slot(x, "AltExp")
  
  # Only show assays that belong to the currently selected AltExp
  all_assays <- if (!is.na(current_ae) && current_ae %in% all_altexps) {
    assayNames(altExp(se, current_ae))
  } else {
    unlist(lapply(all_altexps, function(nm) assayNames(altExp(se, nm))))
  }
  
  # LookupColumn: rowname sentinel + character/factor columns in the main 
  # experiment rowData.
  # The value for the selected feature in this column is the join key.
  main_rd_cls   <- sapply(rowData(se), class)
  valid_lookcols <- c("(rowname)", 
                      names(main_rd_cls)[main_rd_cls %in% 
                                           c("character", "factor")])
  
  # MapColumn: character/factor columns in the altExp rowData.
  # AltExp rows whose value matches the lookup key are shown.
  ae_rd_cls     <- sapply(rowData(altExp(se, current_ae)), class)
  valid_mapcols <- names(ae_rd_cls)[ae_rd_cls %in% c("character", "factor")]
  
  tab_by_row        <- select_info$single$feature
  column_covariates <- .getCachedCommonInfo(se, "ColumnDotPlot")$valid.colData.names
  
  xaxis_choices <- iSEE:::.featAssayXAxisNothingTitle
  if (length(column_covariates)) {
    xaxis_choices <- c(xaxis_choices, iSEE:::.featAssayXAxisColDataTitle)
  }
  
  # Tour steps for the y-axis feature selector and its linked controls
  .addSpecificTour(class(x)[1], 
                   iSEE:::.featAssayYAxisFeatName, 
                   function(plot_name) {
    data.frame(rbind(
      c(
        element = paste0("#", plot_name, "_", iSEE:::.featAssayYAxisFeatName,
                         " + .selectize-control"),
        intro = "Here, we choose the feature from the <em>main</em> experiment
whose row annotation will be used to look up matching rows in the selected
alternative experiment.  The assay values of those altExp rows — not those of
the main-experiment feature itself — are what appears on the y-axis."
      ),
      c(
        element = paste0("#", plot_name, "_", iSEE:::.featAssayYAxisRowTable,
                         " + .selectize-control"),
        intro = "The feature selection can also be driven automatically by a
selection made in another panel, for example a <em>Row Data Table</em>.
Choose the transmitting panel here so that clicking a row in that table
immediately updates this plot."
      ),
      c(
        element = paste0("#", plot_name, "_", iSEE:::.featAssayYAxisFeatDynamic),
        intro = "If this box is checked, the panel will react to a feature
selection from <em>any</em> row-based panel in the application rather than
the one panel named above.  This is useful when you want to browse many
features without reconfiguring the link each time."
      )
    ), stringsAsFactors = FALSE)
  })
  
  # Tour steps for the AltExp-specific data controls
  .addSpecificTour(class(x)[1], "AltExp", function(plot_name) {
    data.frame(rbind(
      c(
        element = paste0("#", plot_name, "_AltExp + .selectize-control"),
        intro = "Select which alternative experiment to visualise.  Alternative
experiments are stored alongside the main <code>SingleCellExperiment</code>
and typically hold complementary modalities such as precursor, PSM and/or peptide 
abundances in MS-based proteomics that were summarised in protein abundances 
stored in the mainExp; CITE-seq protein counts,
ATAC-seq peaks; or CRISPR guide capture data."),
      c(
        element = paste0("#", plot_name, "_AltAssay + .selectize-control"),
        intro = "Choose which assay within the selected alternative experiment
to display on the y-axis — for example <code>peptide_norm</code> for
log2-normalised values or <code>peptide_log</code> for log2 transformed and non 
normalised values, etc.  Only assays belonging to the currently selected 
        alternative experiment are listed."
      ),
      c(
        element = paste0("#", plot_name, "_LookupColumn + .selectize-control"),
        intro = "Choose a <code>character</code> or <code>factor</code> column
from the <em>main experiment</em> <code>rowData</code>.  The value stored for
the selected feature in this column is used as the join key — for example
<code>\"Proteins\"</code> when the main experiment holds protein-level data
identified by a protein name column."
      ),
      c(
        element = paste0("#", plot_name, "_MapColumn + .selectize-control"),
        intro = "Choose a <code>character</code> or <code>factor</code> column
from the <em>altExp</em> <code>rowData</code>.  AltExp rows whose value in
this column matches the key read from the Lookup column are shown — for
example <code>\"Proteins\"</code> in a peptide-level altExp."
      ),
      c(
        element = paste0("#", plot_name, "_PlotType"),
        intro = "Controls how the altExp measurements are rendered.
<strong>Auto</strong> lets iSEE choose based on the x-axis type.
<strong>Scatter</strong> always produces a dot plot.
<strong>Scatter + lines</strong> additionally connects the dots for each
altExp feature with a line, which is most informative when the x-axis encodes
an ordered variable such as pseudotime or sampleIds ordered according to 
        treaments."
      )
    ), stringsAsFactors = FALSE)
  })
  
  # Tour steps for the x-axis controls
  .addSpecificTour(class(x)[1], iSEE:::.featAssayXAxis, function(plot_name) {
    data.frame(rbind(
      c(
        element = paste0("#", plot_name, "_", iSEE:::.featAssayXAxis),
        intro = "Choose what to show on the x-axis.  Leaving this as
<strong>None</strong> produces a single column of points (useful for a simple
strip chart or violin).  Selecting <strong>Column data</strong> lets you
stratify samples by a sample-level covariate."
      ),
      c(
        element = paste0("#", plot_name, "_", iSEE:::.featAssayXAxis),
        intro = "If you select <strong>Column data</strong>..."
      ),
      c(
        element = paste0("#", plot_name, "_", iSEE:::.featAssayXAxisColData,
                         " + .selectize-control"),
        intro = "... you can pick any valid <code>colData</code> field to
place on the x-axis.  Numeric fields produce a scatter plot; factor or
character fields produce grouped strips or violins.  Combined with
<strong>Scatter + lines</strong>, trajectories connecting each altExp feature 
        across samples."
      )
    ), stringsAsFactors = FALSE)
  })
  
  list(
    .selectizeInput.iSEE(
      x, iSEE:::.featAssayYAxisFeatName,
      label    = "Feature in main assay:",
      choices  = NULL,
      selected = NULL,
      multiple = FALSE
    ),
    selectInput(
      .input_FUN(iSEE:::.featAssayYAxisRowTable),
      label    = NULL,
      choices  = tab_by_row,
      selected = iSEE:::.choose_link(slot(x, iSEE:::.featAssayYAxisRowTable), 
                                     tab_by_row)
    ),
    checkboxInput(
      .input_FUN(iSEE:::.featAssayYAxisFeatDynamic),
      label = "Use dynamic feature selection for the y-axis",
      value = slot(x, iSEE:::.featAssayYAxisFeatDynamic)
    ),
    selectInput(
      .input_FUN("AltExp"),
      label    = "AltExp:",
      choices  = all_altexps,
      selected = iSEE:::.choose_link(slot(x, "AltExp"), all_altexps)
    ),
    selectInput(
      .input_FUN("AltAssay"),
      label    = "Alt assay:",
      choices  = all_assays,
      selected = iSEE:::.choose_link(slot(x, "AltAssay"), all_assays)
    ),
    selectInput(
      .input_FUN("LookupColumn"),
      label    = "Lookup column (main assay):",
      choices  = valid_lookcols,
      selected = iSEE:::.choose_link(slot(x, "LookupColumn"), valid_lookcols)
    ),
    selectInput(
      .input_FUN("MapColumn"),
      label    = "Map column (altExp):",
      choices  = valid_mapcols,
      selected = iSEE:::.choose_link(slot(x, "MapColumn"), valid_mapcols)
    ),
    selectInput(
      .input_FUN("PlotType"),
      label    = "Plot type:",
      choices  = c("Auto", "Scatter", "Scatter + lines"),
      selected = slot(x, "PlotType")
    ),
    .radioButtons.iSEE(
      x, iSEE:::.featAssayXAxis,
      label    = "X-axis:",
      inline   = TRUE,
      choices  = xaxis_choices,
      selected = slot(x, iSEE:::.featAssayXAxis)
    ),
    .conditionalOnRadio(
      .input_FUN(iSEE:::.featAssayXAxis),
      iSEE:::.featAssayXAxisColDataTitle,
      selectInput(
        .input_FUN(iSEE:::.featAssayXAxisColData),
        label    = "X-axis column data:",
        choices  = column_covariates,
        selected = slot(x, iSEE:::.featAssayXAxisColData)
      )
    ),
    .conditionalOnRadio(
      .input_FUN(iSEE:::.featAssayXAxis),
      iSEE:::.featAssayXAxisFeatNameTitle,
      selectizeInput(
        .input_FUN(iSEE:::.featAssayXAxisFeatName),
        label    = "X-axis feature:",
        choices  = NULL,
        selected = NULL,
        multiple = FALSE
      ),
      selectInput(
        .input_FUN(iSEE:::.featAssayXAxisRowTable),
        label    = NULL,
        choices  = tab_by_row,
        selected = slot(x, iSEE:::.featAssayXAxisRowTable)
      ),
      checkboxInput(
        .input_FUN(iSEE:::.featAssayXAxisFeatDynamic),
        label = "Use dynamic feature selection for the x-axis",
        value = slot(x, iSEE:::.featAssayXAxisFeatDynamic)
      )
    )
  )
})

#' @importFrom iSEE .getEncodedName .createUnprotectedParameterObservers
setMethod(".createObservers", 
          "AltExpPlot", 
          function(x, se, input, session, pObjects, rObjects) {
            callNextMethod()
            
            plot_name <- .getEncodedName(x)
            
            .createUnprotectedParameterObservers(
              plot_name,
              fields=c("AltExp", 
                       "AltAssay", 
                       "LookupColumn", 
                       "MapColumn", 
                       "PlotType"),
              input=input,
              pObjects=pObjects,
              rObjects=rObjects)
            
            invisible(NULL)
          })


#' @importFrom iSEE .featAssayYAxisFeatName .featAssayXAxis
#'   .featAssayXAxisColDataTitle .featAssayXAxisColData .textEval
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment rowData assay colData
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
setMethod(".generateDotPlotData", "AltExpPlot", function(x, envir) {
  data_cmds <- list()
  
  gene_selected_y <- slot(x, iSEE:::.featAssayYAxisFeatName)
  lookup_col      <- slot(x, "LookupColumn")
  map_col         <- slot(x, "MapColumn")
  alt_exp_name    <- slot(x, "AltExp")
  alt_assay_name  <- slot(x, "AltAssay")
  
  pg_cmd <- if (identical(lookup_col, "(rowname)")) {
    sprintf("pg <- %s;", deparse(gene_selected_y))
  } else {
    sprintf("pg <- rowData(se)[%s, %s];", deparse(gene_selected_y), deparse(lookup_col))
  }

  data_cmds[["y"]] <- c(
    # Retrieve the altExp
    sprintf("ae <- altExp(se, %s);", deparse(alt_exp_name)),

    # Join key: the feature's rowname when "(rowname)" is selected, 
    # otherwise a rowData value
    pg_cmd,
    
    # Guard: stop early with an informative message
    "if (length(pg) == 0 || all(is.na(pg))) stop('No value found in the Lookup column for the selected feature. Check the Lookup column setting.');",
    
    # Subset altExp rows whose MapColumn value matches the lookup key
    sprintf("alt <- ae[rowData(ae)[[%s]] %%in%% pg, ];", deparse(map_col)),
    
    # Pivot to long format; 'altExp_feature_id' is the altExp feature identifier
    sprintf(paste0(
      "plot.data <- as.data.frame.matrix(assay(alt, %s)) |>",
      " tibble::rownames_to_column('altExp_feature_id') |>",
      " tidyr::pivot_longer(names_to = 'sample', values_to = 'Y', -'altExp_feature_id');"
    ), deparse(alt_assay_name))
  )
  
  x_choice <- slot(x, iSEE:::.featAssayXAxis)
  
  if (x_choice == iSEE:::.featAssayXAxisColDataTitle) {
    # Index into colData by sample name so row order always matches,
    # regardless of how pivot_longer arranged the rows.
    data_cmds[["x"]] <- sprintf(
      "plot.data$X <- colData(se)[plot.data$sample, %s];",
      deparse(slot(x, iSEE:::.featAssayXAxisColData))
    )
  } else {
    data_cmds[["x"]] <- "plot.data$X <- as.factor('');"
  }
  
  data_cmds <- unlist(data_cmds)
  .textEval(data_cmds, envir)

  pg_value <- paste(get("pg", envir = envir), collapse = ", ")
  plot_title <- if (pg_value == gene_selected_y) pg_value else sprintf("%s (%s)", pg_value, gene_selected_y)

  list(
    commands = data_cmds,
    labels   = list(title = plot_title, X = "", Y = alt_assay_name)
  )
})

############################################################
# Adjust .addDotPlotDataXXX methods for multiple features
############################################################

# index colData directly by plot.data$sample


#' @importFrom iSEE .colorByField .colorByColDataTitle .colorByFeatNameTitle
#'   .colorByFeatName .colorByFeatNameAssay .colorBySampNameTitle
#'   .colorBySampName .colorByColSelectionsTitle .textEval
#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr n_distinct
setMethod(".addDotPlotDataColor", "AltExpPlot", function(x, envir) {
  color_choice <- slot(x, iSEE:::.colorByField)
  
  if (color_choice == iSEE:::.colorByColDataTitle) {
    covariate_name <- slot(x, iSEE:::.colorByColData)
    label <- covariate_name
    
    cmds <- if (covariate_name == "altExp_feature_id") {
      # 'altExp_feature_id' lives in plot.data itself, not in colData
      "plot.data$ColorBy <- plot.data$altExp_feature_id;"
    } else {
      sprintf(
        "plot.data$ColorBy <- colData(se)[plot.data$sample, %s];",
        deparse(covariate_name)
      )
    }
    
  } else if (color_choice == iSEE:::.colorByFeatNameTitle) {
    chosen_gene  <- slot(x, iSEE:::.colorByFeatName)
    assay_choice <- slot(x, iSEE:::.colorByFeatNameAssay)
    label <- sprintf("%s\n(%s)", chosen_gene, assay_choice)
    cmds  <- sprintf(
      "plot.data$ColorBy <- rep(assay(se, %s)[%s, ], dplyr::n_distinct(plot.data$altExp_feature_id));",
      deparse(assay_choice), deparse(chosen_gene)
    )
    
  } else if (color_choice == iSEE:::.colorBySampNameTitle) {
    chosen_sample <- slot(x, iSEE:::.colorBySampName)
    label <- chosen_sample
    cmds  <- sprintf(
      "plot.data$ColorBy <- logical(nrow(plot.data));\nplot.data[plot.data$sample == %s, 'ColorBy'] <- TRUE;",
      deparse(chosen_sample)
    )
    
  } else if (color_choice == iSEE:::.colorByColSelectionsTitle) {
    label  <- "Column selection"
    target <- if (exists("col_selected", envir = envir, inherits = FALSE)) {
      "col_selected"
    } else {
      "list()"
    }
    cmds <- sprintf(
      "plot.data$ColorBy <- iSEE::multiSelectionToFactor(%s, colnames(se));",
      target
    )
    
  } else {
    return(NULL)
  }
  
  iSEE:::.textEval(cmds, envir)
  list(commands = cmds, labels = list(ColorBy = label))
})


#' @importFrom iSEE .shapeByField .shapeByColDataTitle .shapeByColData .textEval
#' @importFrom SummarizedExperiment colData
setMethod(".addDotPlotDataShape", "AltExpPlot", function(x, envir) {
  shape_choice   <- slot(x, iSEE:::.shapeByField)
  covariate_name <- slot(x, iSEE:::.shapeByColData)
  
  if (shape_choice == iSEE:::.shapeByColDataTitle) {
    
    if (covariate_name == "altExp_feature_id") {
      label <- "altExp_feature_id"
      cmds  <- "plot.data$ShapeBy <- plot.data$altExp_feature_id;"
      
    } else if (covariate_name %in% colnames(colData(envir$se))) {
      label <- covariate_name
      cmds  <- sprintf(
        "plot.data$ShapeBy <- colData(se)[plot.data$sample, %s];",
        deparse(covariate_name)
      )
      
    } else {
      # Covariate no longer present — fall back silently
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


#' @importFrom iSEE .sizeByField .sizeByColDataTitle .sizeByColData .textEval
#' @importFrom SummarizedExperiment colData
setMethod(".addDotPlotDataSize", "AltExpPlot", function(x, envir) {
  size_choice <- slot(x, iSEE:::.sizeByField)
  
  if (size_choice == iSEE:::.sizeByColDataTitle) {
    covariate_name <- slot(x, iSEE:::.sizeByColData)
    label <- covariate_name
    cmds  <- sprintf(
      "plot.data$SizeBy <- colData(se)[plot.data$sample, %s];",
      deparse(covariate_name)
    )
  } else {
    return(NULL)
  }
  
  iSEE:::.textEval(cmds, envir)
  list(commands = cmds, labels = list(SizeBy = label))
})

#' @importFrom iSEE .facetRow .facetColumn .facetRowByColData
#'   .facetColumnByColData .facetByColDataTitle .facetByColSelectionsTitle
#'   .textEval
#' @importFrom SummarizedExperiment colData
setMethod(".addDotPlotDataFacets", "AltExpPlot", function(x, envir) {
  facet_cmds <- NULL
  labels      <- list()
  
  params <- list(
    list(iSEE:::.facetRow,    "FacetRow",    iSEE:::.facetRowByColData),
    list(iSEE:::.facetColumn, "FacetColumn", iSEE:::.facetColumnByColData)
  )
  
  for (f in seq_len(2)) {
    current    <- params[[f]]
    param_field <- current[[1]]
    pd_field    <- current[[2]]
    facet_mode  <- slot(x, param_field)
    
    if (facet_mode == iSEE:::.facetByColDataTitle) {
      facet_data <- x[[current[[3]]]]
      facet_cmds[pd_field] <- sprintf(
        "plot.data$%s <- colData(se)[plot.data$sample, %s];",
        pd_field, deparse(facet_data)
      )
      labels[[pd_field]] <- facet_data
      
    } else if (facet_mode == iSEE:::.facetByColSelectionsTitle) {
      target <- if (exists("col_selected", envir = envir, inherits = FALSE)) {
        "col_selected"
      } else {
        "list()"
      }
      facet_cmds[pd_field] <- sprintf(
        "plot.data$%s <- iSEE::multiSelectionToFactor(%s, colnames(se));",
        pd_field, target
      )
      labels[[pd_field]] <- "Column selection"
    }
  }
  
  iSEE:::.textEval(facet_cmds, envir)
  list(commands = facet_cmds, labels = labels)
})



############################################################
# Visual UI overrides
############################################################

### Redefine .defineVisual interfaces
### Reduce color options as compared to DotPlot class


#' @importFrom iSEE .getCachedCommonInfo .getEncodedName .colorByField
#'   .colorByNothingTitle .colorByColDataTitle .allowableColorByDataChoices
#'   .getDotPlotColorConstants .singleSelectionDimension .radioButtons.iSEE
#'   .conditionalOnRadio .colorByDefaultColor .choose_link .selectTransAlpha
#'   .sliderInput.iSEE
#' @importFrom colourpicker colourInput
#' @importFrom shiny hr tagList selectInput selectizeInput checkboxInput
setMethod(".defineVisualColorInterface", 
          "AltExpPlot", 
          function(x, se, select_info) {
            all_assays <- iSEE:::.getCachedCommonInfo(se, "AltExpPlot")$valid.assay.names
            
            plot_name    <- iSEE:::.getEncodedName(x)
            colorby_field <- paste0(plot_name, "_", iSEE:::.colorByField)
            
            colorby        <- iSEE:::.getDotPlotColorConstants(x)
            mydim_single   <- iSEE:::.singleSelectionDimension(x)
            otherdim_single <- setdiff(c("feature", "sample"), mydim_single)
            mydim_choices  <- select_info[[mydim_single]]
            otherdim_choices <- select_info[[otherdim_single]]
            
            covariates    <- iSEE:::.allowableColorByDataChoices(x, se)
            color_choices <- iSEE:::.colorByNothingTitle
            if (length(covariates)) {
              color_choices <- c(color_choices, iSEE:::.colorByColDataTitle)
            }
            
            tagList(
              hr(),
              iSEE:::.radioButtons.iSEE(
                x, iSEE:::.colorByField,
                label    = "Color by:",
                inline   = TRUE,
                choices  = color_choices,
                selected = slot(x, iSEE:::.colorByField)
              ),
              iSEE:::.conditionalOnRadio(
                colorby_field, iSEE:::.colorByNothingTitle,
                colourpicker::colourInput(
                  paste0(plot_name, "_", iSEE:::.colorByDefaultColor),
                  label = NULL,
                  value = slot(x, iSEE:::.colorByDefaultColor)
                )
              ),
              iSEE:::.conditionalOnRadio(
                colorby_field, colorby$metadata$title,
                selectInput(
                  paste0(plot_name, "_", colorby$metadata$field),
                  label    = NULL,
                  choices  = c("altExp_feature_id", iSEE:::.allowableColorByDataChoices(x, se)),
                  selected = x[[colorby$metadata$field]]
                )
              ),
              iSEE:::.conditionalOnRadio(
                colorby_field, colorby$name$title,
                selectizeInput(
                  paste0(plot_name, "_", colorby$name$field),
                  label = NULL, selected = NULL, choices = NULL, multiple = FALSE
                ),
                selectInput(
                  paste0(plot_name, "_", colorby$name$table),
                  label    = NULL,
                  choices  = mydim_choices,
                  selected = iSEE:::.choose_link(x[[colorby$name$table]], mydim_choices)
                ),
                colourpicker::colourInput(
                  paste0(plot_name, "_", colorby$name$color),
                  label = NULL,
                  value = x[[colorby$name$color]]
                ),
                checkboxInput(
                  paste0(plot_name, "_", colorby$name$dynamic),
                  label = sprintf("Use dynamic %s selection", mydim_single),
                  value = x[[colorby$name$dynamic]]
                )
              ),
              iSEE:::.conditionalOnRadio(
                colorby_field, colorby$assay$title,
                selectizeInput(
                  paste0(plot_name, "_", colorby$assay$field),
                  label = NULL, choices = NULL, selected = NULL, multiple = FALSE
                ),
                selectInput(
                  paste0(plot_name, "_", colorby$assay$assay),
                  label    = NULL,
                  choices  = all_assays,
                  selected = x[[colorby$assay$assay]]
                ),
                selectInput(
                  paste0(plot_name, "_", colorby$assay$table),
                  label    = NULL,
                  choices  = otherdim_choices,
                  selected = iSEE:::.choose_link(x[[colorby$assay$table]], otherdim_choices)
                ),
                checkboxInput(
                  paste0(plot_name, "_", colorby$assay$dynamic),
                  label = sprintf("Use dynamic %s selection", otherdim_single),
                  value = x[[colorby$assay$dynamic]]
                )
              ),
              iSEE:::.sliderInput.iSEE(
                x, iSEE:::.selectTransAlpha,
                label = "Unselected point opacity:",
                min   = 0, max = 1,
                value = slot(x, iSEE:::.selectTransAlpha)
              )
            )
          })



# Custom Shape interface: adds "altExp_feature_id" as a shape-by option in 
# addition to the standard colData discrete covariates provided by the DotPlot 
# parent class.
#' @importFrom iSEE .getEncodedName .shapeByField .shapeByNothingTitle
#'   .getDiscreteMetadataChoices .getDotPlotShapeConstants .radioButtons.iSEE
#'   .conditionalOnRadio
#' @importFrom shiny hr tagList selectInput

setMethod(".defineVisualShapeInterface", "AltExpPlot", function(x, se) {
  discrete_covariates <- c("altExp_feature_id", 
                           iSEE:::.getDiscreteMetadataChoices(x, se))
  
  if (length(discrete_covariates)) {
    plot_name    <- iSEE:::.getEncodedName(x)
    shapeby_field <- paste0(plot_name, "_", iSEE:::.shapeByField)
    shapeby       <- iSEE:::.getDotPlotShapeConstants(x)
    
    tagList(
      hr(),
      iSEE:::.radioButtons.iSEE(
        x, iSEE:::.shapeByField,
        label    = "Shape by:",
        inline   = TRUE,
        choices  = c(iSEE:::.shapeByNothingTitle, shapeby$metadata$title),
        selected = slot(x, iSEE:::.shapeByField)
      ),
      iSEE:::.conditionalOnRadio(
        shapeby_field, shapeby$metadata$title,
        selectInput(
          paste0(plot_name, "_", shapeby$metadata$field),
          label    = NULL,
          choices  = discrete_covariates,
          selected = x[[shapeby$metadata$field]]
        )
      )
    )
  } else {
    NULL
  }
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

setMethod(".generateDotPlot", "AltExpPlot", function(x, labels, envir) {
  plot_data <- envir$plot.data
  
  is_subsetted  <- exists("plot.data.all", envir = envir, inherits = FALSE)
  is_downsampled <- exists("plot.data.pre", envir = envir, inherits = FALSE)
  
  # PlotType "Auto" defers to whatever iSEE chose in plot.type;
  # "Scatter" and "Scatter + lines" both use scatter2 as the base renderer.
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
                      square = do.call(iSEE:::.square_plot,  args),
                      violin = do.call(iSEE:::.violin_plot,  args),
                      violin_horizontal = do.call(iSEE:::.violin_plot, 
                                                  c(args, 
                                                    list(horizontal = TRUE)
                                                  )
                      ),
                      scatter = do.call(iSEE:::.scatter_plot, args),
                      scatter2 = do.call(.scatter_plot2, args)
  )
  
  # Shape scale: use n_distinct to get the right number of shapes
  if (slot(x, iSEE:::.shapeByField) != iSEE:::.shapeByNothingTitle) {
    n_shapes <- dplyr::n_distinct(plot_data$ShapeBy)
    if (n_shapes < 26) {
      N <- length(plot_cmds)
      plot_cmds[[N]]     <- paste(plot_cmds[[N]], "+")
      plot_cmds[[N + 1]] <- sprintf(
        "scale_shape_manual(values = seq_len(%d))", n_shapes
      )
    }
  }
  
  # Line layer — only when explicitly requested
  if (add_lines) {
    N <- length(plot_cmds)
    color_set <- !is.null(plot_data$ColorBy)
    aes_line  <- iSEE:::.buildAes(
      color = color_set, group = TRUE,
      alt   = c(color = iSEE:::.set_colorby_when_none(x))
    )
    # Map GroupBy --> altExp_feature_id (the altExp feature identifier)
    aes_line <- gsub("GroupBy", "altExp_feature_id", aes_line, fixed = TRUE)
    plot_cmds[[N]] <- paste(plot_cmds[[N]], "+")
    plot_cmds[["geom_line"]] <- sprintf(
      "geom_line(data = plot.data, mapping = %s)", 
      aes_line
    )
  }
  
  # Faceting
  facet_cmd <- iSEE:::.addFacets(x)
  if (length(facet_cmd)) {
    N <- length(plot_cmds)
    plot_cmds[[N]] <- paste(plot_cmds[[N]], "+")
    plot_cmds <- c(plot_cmds, facet_cmd)
  }
  
  plot_cmds <- iSEE:::.addCustomLabelsCommands(
    x, 
    commands = plot_cmds,
    plot_type = plot_type
  )
  
  if (plot_type == "scatter") {
    plot_cmds <- iSEE:::.addLabelCentersCommands(x, commands = plot_cmds)
  }
  
  plot_cmds <- iSEE:::.addMultiSelectionPlotCommands(
    x,
    flip   = (plot_type == "violin_horizontal"),
    envir  = envir,
    commands = plot_cmds
  )

  if (dplyr::n_distinct(plot_data$altExp_feature_id) > 10) {
    N <- length(plot_cmds)
    plot_cmds[[N]] <- paste(plot_cmds[[N]], "+")
    plot_cmds[["no_legend"]] <- "theme(legend.position = 'none')"
  }

  list(plot = iSEE:::.textEval(plot_cmds, envir), commands = plot_cmds)
})


############################################################
# Panel tour
############################################################

#' @importFrom iSEE .getEncodedName .getPanelColor .addTourStep .dataParamBoxOpen
#'   .definePanelTour
setMethod(".definePanelTour", "AltExpPlot", function(x) {
  collated <- rbind(
    c(
      element = paste0("#", .getEncodedName(x)),
      intro   = sprintf(
        "The <font color=\"%s\">AltExp plot</font> panel displays assay values
from an <em>alternative experiment</em> stored in a
<code>SingleCellExperiment</code> object — for example precursor, PSM and/or 
peptide log2-intensities for MS-based protemics; CITE-seq protein
measurements; ATAC-seq peaks; or CRISPR guide capture data.
<br><br>
A feature from the <em>main</em> experiment is selected, and the panel uses a
row annotation columns (<em>Map column</em>) from the AltExp and the MainExp
(<em> Lookup column<em>) to retrieve the corresponding row(s) in the alternative 
experiment.  All matching altExp features are shown as separate points (or 
lines) per sample, making it easy to compare related measurements from different 
modalities side by side.",
        .getPanelColor(x)
      )
    ),
    c(
      element = paste0("#", .getEncodedName(x)),
      intro   = "Each point in the panel represents one
<strong>(altExp feature, sample)</strong> combination.  When multiple altExp
rows match the selected main-experiment feature — for example several
precursors mapping to the same protein, or several isoforms or several 
antibody-derived tags mapping to the same gene symbol —
all of them appear simultaneously, distinguished by colour, shape, or a
connecting line depending on your visual settings."
    ),
    .addTourStep(x, .dataParamBoxOpen,
                 "The <i>Data parameters</i> box contains all controls specific 
                 to this panel.<br><br><strong>Action:</strong> click on this 
                 box to expand it and explore the available options."
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
