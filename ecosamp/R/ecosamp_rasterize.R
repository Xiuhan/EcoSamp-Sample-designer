#' EcoSamp rasterizing function
#'
#' This function is a sub-function of ecosamp() function, which doesn't work on its own.
#' It loads an input habitat/treatment map with a SF object class, to be converted into RasterLayer,
#' and desired resolution and reference for text-numeric transformation.
#' 
#' Using the inputs, it converts the map into a RasterLayer object.
#' When the map is in geographic CRS, it automatically transforms it into local UTM CRS.
#' If sd_crs was not specified in ecosamp(), this CRS will be used as the target CRS.
#' 
#' It then outputs the transformed RasterLayer map.
#' 
#' @param map SP object, the target map needing transformation, specified in the ecosamp() function.
#' @param index Data frame, a dataframe for reference of text-numeric transformation of SP maps.
#' It has 2 columns: "Field" contains all the unique values in the input map, 
#' and "Index" contains a unique numeric value for each of them. 
#' Rows with value of 0 in Index will be excluded from sample design.
#' @param map_res Numeric, resolution (in metres) of sample map processing, specified in the ecosamp() function.
#' 5(m) by default. 
#' 
#' @return A RasterLayer transformed version of input map.
#' @export
ecosamp_rasterize <- function(map_res = 5, map,index){
  # Rename column names of input map for standarsation
  colnames(map) <- c("Field","geometry")
  # Merge map and numerical index for conversion into Raster object
  map <- merge(map,index,by="Field")
  
  # Check whether CRS is geographic (lon/lat)
  if (sf::st_is_longlat(map)) {
    # Determine appropriate UTM zone from the centre of the map
    centre <- sf::st_coordinates(sf::st_centroid(sf::st_union(map)))
    lon <- centre[1]
    lat <- centre[2]
    utm_zone <- floor((lon + 180) / 6) + 1
    # Compute EPSG code:
    epsg <- if (lat >= 0) {
      32600 + utm_zone
    } else {
      32700 + utm_zone
    }
    # Transform to UTM
    map <- sf::st_transform(map, epsg)
    message("Input was in geographic CRS. Transformed to local UTM CRS.")
  }
  
  # Buffer polygon by 1000 m
  map_buff <- sf::st_buffer(map, dist = 1000)
  # Get extent of buffered polygon
  bb <- sf::st_bbox(map)
  # Get extent of polygon input
  ext <- raster::extent(bb["xmin"],bb["xmax"],bb["ymin"],bb["ymax"])
  # Create an empty background for the sample area (resolution = 5m)
  base <- raster::raster(ext,res = map_res,crs = raster::crs(map))
  # Convert SF objects into raster objects using the background as template
  map_out <- fasterize::fasterize(map, base,field = "Index",fun = "max")
  #plot(map_out)
  return(map_out)
}