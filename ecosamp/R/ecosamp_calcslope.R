#' EcoSamp slope angle calculation function
#'
#' This function is a sub-function of ecosamp() function, which doesn't work on its own.
#' It loads an input elevation map with a RasterLayer or SpatRaster object class, 
#' and calculates regional slope angles using elevation provided.
#' 
#' Using the inputs, it converts the map into a SpatRaster object with regional slope angles.
#' 
#' It then outputs the transformed SpatRaster map.
#' 
#' @param map RasterLayer or SpatRaster object, the elevation map needing transformation, specified in the ecosamp() function.
#' 
#' @return A SpatRaster object with regional slope angles calculated using elevation map provided.
#' @export
ecosamp_calcslope <- function(map){
  # Check class and transform elevation map input
  if (is(map, "RasterLayer")) {
    map <- terra::rast(map)
  } else if (is(map, "SpatRaster")) {
    map <- map
  } else {
    stop("map_elevation must be a RasterLayer or SpatRaster object.")
  }
  
  # Calculate slope angle using transformed elevation map
  map_slope <- terra::terrain(map, v="slope",neighbors=8, unit='degrees')
  
  # Return slope map to the main function
  return(map_slope)

}