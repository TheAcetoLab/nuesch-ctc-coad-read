# Functions to apply to a Seurat objects



#' generate a ComplexHeatmap object with feature expresion from a Seurat object
#' gene set collection list. Returns a similarity matrix
##' @title Feature expression heatmap using a Seurat object as input
##' @description Draws a heatmap of single cell feature expression using CompleHeatmap functions.
##' @param object Seurat object
##' @param features A vector of features to plot, defaults to VariableFeatures(object = object).
##' @param features.group A vector of feature groups to use for row splitting in the heatmap.
##' @param cells  A vector of cells to plot, defaults to all cells in the object.
##' @param slot A character string indicating the slot to use for the heatmap, either "data" or "scale.data".
##' @param assay A character string indicating the assay to use for the heatmap, defaults to "SCT".
##' @param group.by A character string indicating the metadata column to use for column splitting in the heatmap.
##' @param group.col A named vector of colors to use for the groups in the heatmap.
##' @param disp.max A numeric value indicating the maximum value for the color scale in the heatmap.
##' @param brewer.col A character string indicating the Brewer color palette to use for the heatmap, defaults to "RdBu".
##' @param ... Additional parameters passed to the ComplexHeatmap::Heatmap function.
##' @template roxygen-template
##' @details
##' This function is a wrapper to rgenerate a heatmap of feature expression from a Seurat object.
##' input and multiple parameters.
##' @return A Heatmap-class object.
##'
##'
heatmap_seurat <- function(
    obj,
    features = NULL,
    features.group = NULL,
    cells = NULL,
    slot = "data",
    assay = "SCT",
    group.by = "ident",
    group.col = NULL,
    disp.max = NULL,
    brewer.col = "RdBu",
    group.labels.gp =  gpar(col = "black", fontsize = 8),
    ...
){
  require(ComplexHeatmap)
  require(circlize)
  library(RColorBrewer)

  # Subset Seurat object
  if(!is.null(cells)){
    obj <- subset(obj, cells = cells)
  }

  if(is.null(features)){
    features <- VariableFeatures(object = obj)
  }
  obj <- subset(obj, features = features)

  # Combine features annotation
  if(!is.null(features.group))
    names(features.group) <- features


  # Extract matrix of expression
  if(slot == 'data'){
    mat <- obj[[assay]]@data
    zmat <- t(apply(mat, 1, scale, center = TRUE, scale = TRUE))
    colnames(zmat) <- colnames(mat)
  } else if(slot == 'scale.data') {
    zmat <- obj[[assay]]@scale.data[features,]
  } else {
    stop("Slot must be either 'data' or 'scale.data'")
  }

  # Remove missing
  zmat[is.na(zmat)] <- 0

  # Reorder features annotation
  if(!is.null(features.group))
    features.group <- features.group[rownames(zmat)]

  # Column Annotation
  if(!is.null(group.by)) {
    column_annotation <- obj@meta.data[,group.by]
    if(is.factor(column_annotation)){
      column_annotation <- as.character(column_annotation)
    }
    use_annot_levels <- unique(column_annotation) %>% sort
    # if(is.factor(column_annotation)){
    #   use_annot_levels <- levels(column_annotation)
    # } else {
    #   column_annotation <- as.factor(column_annotation)
    #   use_annot_levels <- unique(column_annotation)
    # }

    if(!is.null(group.col)){

      column_ha <-  HeatmapAnnotation(
        condition = anno_block(
          gp = gpar(fill = group.col),
          labels = use_annot_levels,
          labels_gp = group.labels.gp
        )
      )
    } else {
      column_ha <-  HeatmapAnnotation(
        condition = anno_block(
          labels = use_annot_levels,
          labels_gp = group.labels.gp
        )
      )
    }

  } else {
    column_ha <- NULL
  }

  # Heatmap color
  zmax <- zmat %>% abs %>% max(na.rm = TRUE)
  if(!is.null(disp.max))
    zmax <- disp.max

  col_fun <-  colorRamp2(
    seq(-zmax, zmax, length.out = 11),
    rev(brewer.pal(11, brewer.col))
  )

  # Generate Heatmap
  ht <- Heatmap(
    zmat,
    name= 'Scaled\nexpression',
    col = col_fun,

    show_column_names = FALSE,
    show_column_dend = FALSE,
    top_annotation = column_ha,
    column_split = column_annotation,
    cluster_column_slices = FALSE,
    column_title = NULL,


    show_row_dend = FALSE,
    row_split = features.group,
    cluster_row_slices = FALSE,
    row_title = NULL,
    ...
  )

  return(ht)
}
