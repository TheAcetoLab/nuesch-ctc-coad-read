# Functions adapted from semla

#' @description Extract barcodes from the borders of
#' a selected defined region. This function was extracted from semla::RadialDistance.default
#'
#' @section Scenario:
#' The region of interest could for example be an isolated tumor in the tissue
#' section surrounded by stroma. If we are interested in expressional changes
#' from the tumor core and outwards, we can use radial distances to model such
#' changes. Below are a few examples for how radial distances can be used to in
#' address certain question in this scenario:
#'
#'
#' @section Algorithm:
#' The border spots of the selected region is identified based on the
#' spatial network of nearest neighbors identified with \code{\link{GetSpatialNetwork}}.
#'
#' @param spots A character vector with spot IDs present \code{object}. These spots typically
#' represent one particular tissue structure identified either by data-driven clustering
#' or by the tissue histology.

ExtractRegionNeighbors <- function(
    object,
    spots,
    verbose = verbose,
    remove_singletons = TRUE,
    ...
) {

  require(cli)
  # Set global variables to NULL
  barcode <- x <- y <- sampleID <- from  <- to <- centroids <- diff_x <- diff_y <- intervals <- r_dist <- nn <- NULL

  # Check object class
  if (!any(class(object) %in% c("data.frame", "matrix", "tbl")))
    abort(glue("Invalid class '{class(object)}'."))
  if (ncol(object) != 4)
    abort(glue("Invalid number of columns '{ncol(object)}'. Expected 4."))
  if (!all(colnames(object) == c("barcode", "x", "y", "sampleID")))
    abort("Required columns are: 'barcode', 'x', 'y', and 'sampleID'")
  if (!all(
    object |> summarize(
      check_barcode = is.character(barcode),
      check_x = is.numeric(x),
      check_y = is.numeric(y),
      check_sample = is.numeric(sampleID)
    ) |>
    unlist()
  )) {
    abort(glue("Invalid column class(es)."))
  }

  # Check that spots are in object
  stopifnot(
    inherits(spots, what = "character"),
    length(spots) > 0,
    all(spots %in% object$barcode)
  )

  # Get spatial network
  spatnet <- GetSpatialNetwork(object, ...)
  # spatnet <- GetSpatialNetwork(object, maxDist = use_max_distance)
  if (length(spatnet) > 1) {
    abort(glue("Default method can only handle 1 tissue section at the time, got {length(spatnet)}"))
  }
  spatnet_region <- spatnet[[1]] |>
    filter(from %in% spots, to %in% spots) |>
    group_by(from) |>
    mutate(nn = n())

  # Check for singletons
  if (remove_singletons) {
    spatnet_region2 <- spatnet_region |> filter(nn > 1)
    spots_filtered <- spatnet_region2$from |> unlist() |> unique()
    if (verbose) cli_alert_info("Removing {length(spots) - length(spots_filtered)} spots with 0 neighbors.")
    spots <- spots_filtered
  }

  if(length(spots) == 0) {
    if (verbose) cli_alert_warning("Found no spots. Returning NA values")
    # res <- rep(NA_real_, nrow(object)) |>
    #   setNames(nm = object$barcode)
    res <- list(
      border_spots = NULL,
      inside_spots = NULL,
      outside_spots = object$barcode
    )
    return(res)
  }

  if (verbose) cli_alert_info("Extracting border spots from a region with {length(spots)} spots")

  # Find border
  border_spots <- RegionNeighbors(spatnet, spots = spots, outer_border = FALSE,
                                  ...)
  # border_spots <- RegionNeighbors(spatnet, spots = spots, outer_border = FALSE)
  inside_spots <- setdiff(spots, border_spots)

  # Get indices for spot groups
  border_spots_indices <- match(border_spots, object$barcode)
  outside_spots_indices <- -match(spots, object$barcode)
  inside_spots_indices <- match(inside_spots, object$barcode)
  if (verbose) {
    cli_alert("  Detected {length(border_spots_indices)} spots on borders")
    cli_alert("  Detected {length(inside_spots_indices)} spots inside borders")
    cli_alert("  Detected {nrow(object) - length(spots)} spots outside borders")
  }

  # Return border barcodes
  res <- list(
    border_spots = border_spots,
    inside_spots = inside_spots,
    outside_spots = setdiff(object$barcode, spots)
  )
  return(res)
}



#' @description Calculates the 3D radial distances to all spots from the borders of
#' a selected defined region.
#'
#' @section Scenario:
#' The region of interest could for example be an isolated tumor in the tissue
#' section surrounded by stroma. If we are interested in expressional changes
#' from the tumor core and outwards, we can use radial distances to model such
#' changes. Below are a few examples for how radial distances can be used to in
#' address certain question in this scenario:
#'
#' @section Algorithm:
#' First, the border spots of the selected region is identified based on the
#' spatial network of nearest neighbors identified with \code{\link{GetSpatialNetwork}}.
#' For each spot outside of this border, the distance is calculated to its nearest border spot.
#' Spots located inside the selected region will have negative distances and
#' spots located outside of the selected region will have positive distances.
#'
#' @section Search interval:
#' The microenvironment of the region of interest might be extremely heterogeneous
#' depending on the direction from its center. For this reason, it can be useful to narrow
#' down the search area by defining a smaller angle interval with \code{angles}. Alternatively,
#' you can split the radial distances into an even number of slices with \code{angles_nbreaks}.
#' When using a predefined search interval, the region of interest (e.g. manual annotation)
#' should not contain multiple spatially disconnected regions. Angles are calculated from
#' center of the region of interest so it only makes sense to investigate one region at the time.
#' You can use \code{\link{DisconnectRegions}} to split a categorical variable that contains
#' multiple spatially disconnected regions.
#'
#' @param object A \code{Seurat} object created with \code{semla},
#' and the results are returned to the \code{meta.data} slot. The metadata \code{meta.data} must
#' contain \code{sampleID}, \code{barcode}, and the 3D coordinates named \code{x}, \code{y},
#' and \code{z}.
#' @param remove_singletons Logical specifying if 'singletons' should be excluded. Spatially disconnected
#' regions are not allowed when a "search interval" is defined, but single spots without neighbors are
#' not detected as disconnected components. Single spots will most likely not interfere when calculating
#' the centroid of the region of interest and can therefore be kept.
#' @param column_name A character specifying the name of a column in your meta data that contains
#'  categorical data, e.g. clusters or manual selections
#' @param selected_groups A character vector to select specific groups in \code{column_name} with.
#' All groups are selected by default, but the common use case is to select a region of interest.
#' @param column_suffix Suffix to column names returned in the Seurat object.
#' @param verbose Print messages
#' @param ... Arguments passed to other methods (\code{ExtractRegionNeighbors} and \code{GetSpatialNetwork}.
#' Of special interest are variables \code{maxDist} and \code{k})

RadialDistance3D.Seurat <- function (
    object,
    column_name,
    selected_groups = "yes",
    column_suffix = NULL,
    verbose = TRUE,
    convert_to_microns = FALSE,
    ...
) {
  require(cli)
  require(dbscan)
  require(purrr)

  # Set global variables to NULL
  barcode <- pxl_col_in_fullres <- pxl_row_in_fullres <- sampleID <- angle <- r_dist <- NULL

  # validate input
  semla:::.check_seurat_object(object)
  semla:::.validate_column_name(object, column_name)
  selected_groups <- semla:::.validate_selected_labels(object, selected_groups, column_name)

  # Set new column name suffix
  lbl_suffix <- ifelse(is.null(column_suffix), "", column_suffix)

  # Select spots
  spots_list <- semla:::.get_spots_list(object, selected_groups, column_name)

  # Get 2D coordinates
  coords_list <- semla:::.get_coords_list(object)

  # Extract border spots for each sample
  group_borders <- foreach(lbl = names(spots_list)) %do% {
    if (verbose) cli_alert_info("Calculating radial distances for group '{lbl}'")
    sample_borders <- foreach(nm = names(coords_list)) %do% {
      if (verbose) cli_alert_info("Running calculations for sample {nm}")
      if (is.null(spots_list[[lbl]][[nm]])) {
        if (verbose) cli_alert_warning("Found no spots for groups '{lbl}' in section {nm}. Returning NA values for section {nm}")
        res <- list(
          border_spots = NULL,
          inside_spots = NULL,
          outside_spots = coords_list[[nm]]$barcode
        )
      } else {
        res <- ExtractRegionNeighbors(
          coords_list[[nm]],
          spots = spots_list[[lbl]][[nm]],
          verbose = TRUE,
          ...
        )
        # res <- ExtractRegionNeighbors(
        #   coords_list[[nm]],
        #   spots = spots_list[[lbl]][[nm]],
        #   verbose = TRUE,
        #   maxDist = 33.80579,
        #   remove_singletons = TRUE
        # )
      }
      return(res)
    }

    # Combine list sample_borders
    list(
      border_spots = purrr::map(sample_borders, 'border_spots') %>% unlist,
      inside_spots = purrr::map(sample_borders, 'inside_spots') %>% unlist,
      outside_spots = purrr::map(sample_borders, 'outside_spots') %>%  unlist
    )
  }

  names(group_borders) <- names(spots_list)

  # Calculate radial distances
  res <- lapply(names(group_borders), function(lbl) {

    if (!requireNamespace("dbscan")) {
      abort(glue("Package {cli::col_br_magenta('dbscan')} is required. Please install it with: \n",
                 "install.packages('dbscan')"))
    }

    use_group_borders <- group_borders[[lbl]]
    border_spots <- use_group_borders$border_spots
    outside_spots <- use_group_borders$outside_spots
    inside_spots <- use_group_borders$inside_spots

    knn_spatial_outside <- kNN(x = object@meta.data[border_spots, c("x", "y", "z")] |> as.matrix(),
                                       query = object@meta.data[outside_spots, c("x", "y", "z")] |> as.matrix(),
                                       k = 1)
    knn_spatial_inside <- kNN(x = object@meta.data[border_spots, c("x", "y", "z")] |> as.matrix(),
                                      query = object@meta.data[inside_spots, c("x", "y", "z")] |> as.matrix(),
                                      k = 1)
    if (verbose) cli_alert_success("Returning radial distances for groups '{lbl}'")

    # Get radial distances
    radial_dists <- setNames(numeric(length = ncol(object)), nm = object$barcode)
    radial_dists_outside <- setNames(knn_spatial_outside$dist[, 1, drop = TRUE], nm = outside_spots)
    radial_dists_inside <- setNames(knn_spatial_inside$dist[, 1, drop = TRUE], nm = inside_spots)
    radial_dists[names(radial_dists_outside)] <- radial_dists_outside
    radial_dists[names(radial_dists_inside)] <- -radial_dists_inside

    # Convert to microns
    if (convert_to_microns) {
      if (!requireNamespace("dbscan")) {
        abort(glue("Package {cli::col_br_magenta('dbscan')} is required. Please install it with: \n",
                   "install.packages('dbscan')"))
      }
      center_to_center_pixel_distances <- sapply(coords_list, function(xy) {
        kNN(x = xy |> select(x, y), k = 1)$dist |> min()
      })
      radial_dists <- radial_dists/(mean(center_to_center_pixel_distances)/100)
    }

    radial_dists <- tibble(barcode = names(radial_dists), radial_dists) |>
      setNames(nm = c("barcode", paste0("r_dist_", lbl, lbl_suffix)))

    # Return radial distances
    return(radial_dists)

  })
  res <- Reduce(\(x, y) left_join(x, y, by = "barcode"), res)

  # Return data to Seurat object
  object@meta.data <- object@meta.data |>
    rownames_to_column("rowname") |>
    left_join(res, by = "barcode") |>
    column_to_rownames("rowname")

  return(object)
}



#' @description Runs run DisconnectRegions.default for each sample separately,
#' allowing adjust the scaling factor of maxDist and other parameteres for each
#' sample. DisconnectRegions.default already does this, but using a fixed scaling
#' factor of 1.2 and cannot be adjusted.

#' @section Method:
#' Takes a \code{tibble} and set of spot IDs and returns a data.frame with the
#' input spot IDs, hotspot IDs, and sampleID.

DisconnectRegions.sample <- function(use_coords, use_spots, scaling_factor = 1.2, ...) {
  sample_levels <- use_coords$sampleID %>% unique %>% sort
  si <- sample_levels[4]
  DisconnectRegions_combined <- foreach(si = sample_levels, .errorhandling = "pass", .combine = rbind) %do% {
    use_coords_sample <- use_coords %>% filter(sampleID == si)
    use_spots_sample <- intersect(use_spots, use_coords_sample$barcode)

    # Get spatial network for the sample
    spatnet_sample <- GetSpatialNetwork(use_coords_sample)

    # Get minimum distance from the spatial network
    min_distance_sample <- foreach(i = spatnet_sample) %do% {
      min(i$distance)
    } %>% unlist %>% min

    # Define max distance for DisconnectRegions
    use_max_distance_sample <- min_distance_sample * scaling_factor

    # Run DisconnectRegions with the defined max distance for the sample
    # safe_fun <- possibly(DisconnectRegions, otherwise = NULL)
    # safe_fun(
    #   use_coords_sample,
    #   use_spots_sample,
    #   maxDist = use_max_distance_sample,
    #   verbose = FALSE)

    res <- try(
      DisconnectRegions(
        use_coords_sample,
        use_spots_sample,
        maxDist = use_max_distance_sample,
        ...),
      silent = TRUE)
    if(class(res) != 'try-error') {
      res <- res %>%
        data.frame(check.names = FALSE) %>%
        set_names('hotspot') %>%
        rownames_to_column('sample_spot') %>%
        mutate(sampleID = si)
      return(res)
    } else {
      NULL
    }
  }
  return(DisconnectRegions_combined)
}

