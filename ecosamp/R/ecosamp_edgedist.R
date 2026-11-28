#' EcoSamp distance-to-habitat-edge function
#'
#' This function is a sub-function of ecosamp() function, which doesn't work on its own.
#' It loads the habitat type map produced by the main function and calculates,
#' for every pixel, the distance to the nearest pixel of a different habitat type
#' (including pixels outside the habitat map, i.e. NA pixels).
#' 
#' For each habitat type, the calculation is restricted to the bounding box of that habitat,
#' grown by one pixel. The nearest pixel of a different habitat is always adjacent to
#' the habitat itself, so this gives exactly the same distances as calculating over the whole map,
#' but is much faster for habitats that cover a small part of the map.
#' Habitat types with value 0 (excluded from sample design) are skipped, and keep a distance of 0.
#' If a habitat type has no pixels of a different type around it, its distance is set to Inf.
#' 
#' @param map_habitat RasterLayer, the habitat type map produced by the ecosamp() function.
#' 
#' @return A RasterLayer with the same extent and resolution as map_habitat,
#' containing the distance (in metres) from each pixel to the nearest habitat edge.
#' Pixels outside the habitat map, and pixels of habitat type 0, have a value of 0.
#' @export
ecosamp_edgedist <- function(map_habitat){
  # Work with a SpatRaster and a plain vector of pixel values (much faster than RasterLayer indexing)
  map_rast <- terra::rast(map_habitat)
  n_row <- terra::nrow(map_rast)
  n_col <- terra::ncol(map_rast)
  values_habitat <- terra::values(map_rast, mat = FALSE)
  # Output distances, 0 by default (pixels outside the map and excluded habitats)
  dist_edge <- numeric(length(values_habitat))
  
  # Loop through each type of habitat in the map, skipping excluded habitats (value 0)
  habitat_types <- sort(unique(values_habitat[!is.na(values_habitat) & values_habitat != 0]))
  for (i in habitat_types) {
    # Pixels of the current habitat, and their row and column numbers
    cells_own <- which(values_habitat == i)
    rows_own <- (cells_own - 1) %/% n_col + 1
    cols_own <- (cells_own - 1) %% n_col + 1
    # Bounding box of the current habitat, grown by one pixel (clamped to the map)
    row_min <- max(1, min(rows_own) - 1)
    row_max <- min(n_row, max(rows_own) + 1)
    col_min <- max(1, min(cols_own) - 1)
    col_max <- min(n_col, max(cols_own) + 1)
    # Cell numbers (in the full map) of every pixel in the window, in row-major order
    cells_window <- as.vector(t(outer((row_min:row_max - 1) * n_col, col_min:col_max, "+")))
    values_window <- values_habitat[cells_window]
    is_own <- !is.na(values_window) & values_window == i
    
    # If there is no other habitat in the window, there is no edge
    if (all(is_own)) {
      dist_edge[cells_window[is_own]] <- Inf
      next
    }
    # Pixels of the current habitat become NA, all other pixels ("foreign habitats") become 1,
    #   so terra::distance() gives each NA pixel its distance to the nearest foreign pixel
    map_window <- map_rast[row_min:row_max, col_min:col_max, drop = FALSE]
    terra::values(map_window) <- ifelse(is_own, NA, 1)
    dist_window <- terra::values(terra::distance(map_window, unit = "m"), mat = FALSE)
    # Store distances of the pixels of the current habitat
    dist_edge[cells_window[is_own]] <- dist_window[is_own]
  }
  
  # Return as a RasterLayer, matching the other maps used in ecosamp()
  map_edge <- map_habitat
  raster::values(map_edge) <- dist_edge
  return(map_edge)
}
