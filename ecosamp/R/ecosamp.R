#' EcoSamp sample designing function
#'
#' This function loads the habitat type map, total number of sample points,
#' and optional inputs including treatment, roads and waters maps, 
#' in addition to optional specification parameters (available below).
#' The function uses the input to determine valid sample regions in your map, 
#' and generates sample points matching your specification parameters.
#'
#' @param n_total Integer, total number of sample points to generate.
#' @param map_habitat RasterLayer, SpatRaster or SF object, habitat type map of sample region.
#' 
#' @param index_habitat (optional) Data frame. If map_habitat is a SF object, 
#' a data frame with 2 columns is required to convert it into a RasterLayer: 
#' a "Field" column containing all unique types of habitats in map_habitat, 
#' and an "Index" column containing a unique numerical value for each of them.
#' Types of habitats assigned Index = 0 will not be considered for sample generation.
#' If left as NULL, an unique positive number will be assigned to each habitat type present in the map.
#' NULL by default.
#' 
#' @param map_treatmt (optional) RasterLayer, SpatRaster or SF object, treatment type map of sample region.
#' @param index_treatmt (optional) Data frame. If map_treatmt is a SF object, 
#' a data frame with 2 columns is required to convert it into a RasterLayer: 
#' a "Field" column containing all unique types of habitats in map_treatmt, 
#' and an "Index" column containing a unique numerical value for each of them.
#' Types of treatments assigned Index = 0 will not be considered for sample generation.
#' If left as NULL, an unique positive number will be assigned to each treatment type present in the map.
#' NULL by default.
#' 
#' @param map_roads (optional) SF object, map of roads in sample region.
#' @param map_waters (optional) SF object, map of waters/hydrography in sample region.
#' @param map_elevation (optional) RasterLayer or SpatRaster, elevation map of sample region.
#' @param map_slope (optional) RasterLayer or SpatRaster, slope angle map of sample region. 
#' Can be calculated from map_elevation if not provided.
#' 
#' @param map_res (optional) Numeric, resolution (in metres) of sample map processing.
#' @param sd_crs (optional) String, UTM CRS code of desired output. 
#' If left as NULL, a UTM CRS matching the region of map_habitat will be automatically assigned.
#' @param latlon_crs (optional) String, an alternative CRS code of desired output. 
#' If left as NULL, WGS 84 (EPSG:4326) (latitude/longitude in degrees) will be used.
#' @param point_number_min (optional) Integer, minimum number of sample points in each bin of habitat and/or treatment. 0 by default.
#' @param point_dist_min (optional) Numeric, minimum distance between sample points (metres). 0 by default.
#' @param edge_dist_min (optional) Numeric, minimum distance from points to habitat edge (metres). 0 by default.
#' @param edge_dist_max (optional) Numeric, maximum distance from points to habitat edge (metres). 
#' If left as 0, no maximum distance will be considered. 0 by default.
#' @param water_dist_min (optional) Numeric, minimum distance from points to waters (metres). 
#' Requires valid input for param map_waters. 0 by default.
#' @param water_dist_max (optional) Numeric, maximum distance from points to waters (metres). 
#' Requires valid input for param map_waters. If left as 0, no maximum distance will be considered. 0 by default.
#' @param road_dist_min (optional) Numeric, minimum distance from points to roads (metres). 
#' Requires valid input for param map_roads. 0 by default.
#' @param road_dist_max (optional) Numeric, maximum distance from points to roads (metres). 
#' Requires valid input for param map_roads. If left as 0, no maximum distance will be considered. 0 by default.
#' @param elevation_min (optional) Numeric, minimum elevation of points (metres).
#' Requires valid input for param map_elevation. If left as 0, it will not be considered. 0 by default.
#' @param elevation_max (optional) Numeric, maximum elevation of points (metres).
#' Requires valid input for param map_elevation. If left as 0, it will not be considered. 0 by default.
#' @param slope_min (optional) Numeric, minimum slope angle to allocate sample points (degrees).
#' Requires valid input for param map_slope. If left as 0, it will not be considered. 0 by default.
#' @param slope_max (optional) Numeric, maximum slope angle to allocate sample points (degrees).
#' Requires valid input for param map_slope. If left as 0, it will not be considered. 0 by default.
#' 
#' @param max_distance_sampling (optional) Logical. If TRUE, sample points will be allocated at the eligible pixel furthest away from existing points;
#' if FALSE, sample points will be randomly allocated within eligible regions.
#' FALSE by default.
#' @param area_weighting (optional) Logical. If TRUE, sample points will be generated in ascending order of area of eligible bins,
#' this reduces likelihood of bins with small area being covered by sample points from neighbouring bins with large area;
#' if FALSE, sample points will be randomly selected from eligible bins.
#' TRUE by default.
#' 
#' @param plot_results (optional) Logical. If TRUE, a map of output sample points will be plotted on back ground of habitat type map. 
#' FALSE by default.
#' 
#' @return A data frame containing sample point IDs, coordinates in specified CRSs, 
#' type of habitat (and treatment, distance to edge, distance to road, distance to water, 
#' and slope angle if specified) of generated sample points.
#' @importFrom dplyr %>% n
#' @importFrom methods is
#' @export
ecosamp <- function(# Required inputs
  n_total,
  map_habitat,
  
  # Optional inputs
  sd_crs = NULL,
  point_number_min = 0,
  point_dist_min = 0,
  edge_dist_min = 0, edge_dist_max = 0,
  water_dist_min = 0, water_dist_max = 0,
  road_dist_min = 0, road_dist_max = 0,
  elevation_min = 0, elevation_max = 0,
  slope_min = 0, slope_max = 0, map_res = 5,
  
  latlon_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
  map_treatmt = NULL,
  map_roads = NULL,map_waters = NULL,map_elevation = NULL, map_slope = NULL,
  index_habitat = NULL, index_treatmt = NULL,
  max_distance_sampling = FALSE, area_weighting = TRUE,
  plot_results = FALSE){
  ###############
  ### 2.4.1. If habitat/treatment maps are polygons, transform them into rasters
  # If an index is not provided, create a default one with all types of habitats
  if (is(map_habitat, "sf")) {
    if (is.null(index_habitat)){
      index_habitat <- data.frame(Field=unique(map_habitat$Habitat),
                                  Index=1:length(unique(map_habitat$Habitat)))
    }
    map_habitat <- ecosamp_rasterize(map = map_habitat, index = index_habitat, 
                                     map_res = map_res)
  } else if (is(map_habitat, "SpatRaster")){
    map_habitat <- raster::raster(map_habitat)
  }
  # Do the same for treatment map (if exists)
  if (exists("map_treatmt") & !is.null(map_treatmt) & is(map_treatmt, "sf")) {
    # If an index is not provided, create a default one with all types of treatments
    if (is.null(index_treatmt)){
      index_treatmt <- data.frame(Field=unique(map_treatmt$Treatment),
                                  Index=1:length(unique(map_treatmt$Treatment)))
    }
    map_treatmt <- ecosamp_rasterize(map = map_treatmt, index = index_treatmt, 
                                     map_res = map_res)
  } else if (is(map_treatmt, "SpatRaster")){
    map_treatmt <- raster::raster(map_treatmt)
  }
  
  # Report CRS system and units used
  if (exists("sd_crs") & is.null(sd_crs)){
    sd_crs <- raster::crs(map_habitat)
    message ("No CRS provided, using regional UTM CRS by default.")
  } else if (exists("sd_crs") & !is.null(sd_crs)){
    if ((sf::st_crs(sd_crs)$IsGeographic==TRUE)){
      stop("Input CRS is in unit of degrees, please use CRS with metric unit.")
    } else if (is.na(sf::st_crs(sd_crs))){
      stop("Input CRS is not a valid CRS code.")
    }
  } 
  print(paste("Current CRS:",sd_crs))
  print(paste("CRS unit:",sf::st_crs(sd_crs)$units_gdal))
  
  ##########
  ## 2.4.2. Check if inputs are correct formats, and standardise their CRS
  print ("Processing habitat type map...")
  map_habitat <- ecosamp_transcrs(map = map_habitat, sd_crs = sd_crs, 
                                  map_res = map_res)
  if (exists("map_treatmt") & !is.null(map_treatmt)){
    print ("Optional input: treatment map detected. Processing...")
    map_treatmt <- ecosamp_transcrs(map = map_treatmt, sd_crs = sd_crs, 
                                    map_res = map_res)
  }
  # If roads map exists and is an SF object, then reproject its CRS
    # If it exists and is not an SF object, return an error message
  if (exists("map_roads") & !is.null(map_roads) & (road_dist_min>0 | road_dist_max>0)){
    print ("Optional input: roads map detected. Processing...")
    if (is(map_roads, "sf")) {
      map_roads <- ecosamp_transcrs(map = map_roads, sd_crs = sd_crs, 
                                    map_res = map_res)
    } else  {
      stop("Input roads map is not SF object")
    }
  } else if ((road_dist_min>0 | road_dist_max>0) & is.null(map_roads)) {
    stop("Minimum/maximum distance to roads specified but no roads map provided. Please provide the roads map.")
  }
  # If water map exists and is an SF object, then reproject its CRS
    # If it exists and is not an SF object, return an error message
  if (exists("map_waters") & !is.null(map_waters) & (water_dist_min>0 | water_dist_max>0)){
    print ("Optional input: water map detected. Processing...")
    if (is(map_waters, "sf")) {
      map_waters <- ecosamp_transcrs(map = map_waters, sd_crs = sd_crs, 
                                     map_res = map_res)
    } else  {
      stop("Input waters map is not SF object")
    }
  } else if ((water_dist_min>0 | water_dist_max>0) & is.null(map_waters)) {
    stop("Minimum/maximum distance to waters specified but no waters map provided. Please provide the waters map.")
  }  
  
  # If slope map exists and is a RasterLayer/SpatRaster object, then reproject its CRS.
  # If it exists and is not a RasterLayer/SpatRaster object, return an error message.
  if (exists("map_slope") & !is.null(map_slope) & (slope_max>0 | slope_min>0)){
    print ("Optional input: slope map detected. Processing...")
    if (is(map_slope, "SpatRaster") | is(map_slope, "RasterLayer")) {
      # Create a temporary habitat map for cropping map_slope
      map_ref <- ecosamp_transcrs(map = map_habitat, sd_crs = raster::crs(map_slope), 
                                  map_res = raster::res(map_slope))
      map_slope <- raster::crop(map_slope,terra::ext(map_ref))
      # Then standardise CRS and class object of map_slope
      map_slope <- ecosamp_transcrs(map = map_slope, sd_crs = sd_crs, 
                                    map_res = map_res)
    } else {
      stop("Input slope map is not RasterLayer or SpatRaster object")
    }
  } else if ((slope_max>0 | slope_min>0) & is.null(map_slope) & !is.null(map_elevation)) {
    # If slope_max/slope_min specified without map_slope, but map_elevation provided, 
    # calculate map_slope from map_elevation.
    print ("Optional input: calculating slope angle from elevation map...")
    map_slope <- ecosamp_calcslope(map = map_elevation)
    # Create a temporary habitat map for cropping map_slope
    map_ref <- ecosamp_transcrs(map = map_habitat, sd_crs = raster::crs(map_slope), 
                                map_res = raster::res(map_slope))
    map_slope <- raster::crop(map_slope,terra::ext(map_ref))
    # Then standardise CRS and class object of map_slope
    map_slope <- ecosamp_transcrs(map = map_slope, sd_crs = sd_crs, 
                                  map_res = map_res)

  } else if ((slope_max>0 | slope_min>0) & is.null(map_slope) & is.null(map_elevation)) {
    # If slope_max/slope_min specified without map_slope nor map_elevation,
    # stop function and return an error message.
    stop("Maximum slope angle specified but no slope/elevation map provided. Please provide a slope/elevation map.")
  }
  
  # If elevation map exists and is a RasterLayer/SpatRaster object, then reproject its CRS
    # If it exists and is not a RasterLayer/SpatRaster object, return an error message
  if (exists("map_elevation") & !is.null(map_elevation) & (elevation_max>0 | elevation_min>0)){
    print ("Optional input: elevation map detected. Processing...")
    if (is(map_elevation, "SpatRaster") | is(map_elevation, "RasterLayer")) {
      # Create a temporary habitat map for cropping map_elevation
      map_ref <- ecosamp_transcrs(map = map_habitat, sd_crs = raster::crs(map_elevation), map_res = raster::res(map_elevation))
      map_elevation <- raster::crop(map_elevation,terra::ext(map_ref))
      # Then standardise CRS and class object of map_elevation
      map_elevation <- ecosamp_transcrs(map = map_elevation, sd_crs = sd_crs, 
                                    map_res = map_res)
    } else {
      stop("Input elevation map is not RasterLayer or SpatRaster object")
    }
  } else if ((elevation_max>0 | elevation_min>0) & is.null(map_elevation)) {
    # If elevation_max/elevation_min specified without map_elevation provided,
    # stop function and return an error message.
    stop("Minimum/maximum elevation specified but no elevation map provided. Please provide a elevation map.")
  }

  # Create a RasterStack for later extraction of values
  if (exists("map_treatmt") & !is.null(map_treatmt)){
    map_refmap <- raster::stack(map_habitat,map_treatmt)
    # Rename columns of reference map stack
    names(map_refmap) <- c("Habitat","Treatment")
  } else {
    map_refmap <- raster::stack(map_habitat)
    names(map_refmap) <- c("Habitat")
  }
  
  ###############
  ## 2.4.3. Calculate distances to input habitat features
  # Calculate distance to roads 
  # Calculate the distance from each pixel on the map to the nearest road vector
  # map_habitat needs to be converted to SpatRaster object first, because "distance doesn't work
  # for RasterLayer object
  # Only perform this section if a min/max distance to roads is specified
  if (exists("map_roads") & !is.null(map_roads) & (road_dist_min>0 | road_dist_max>0)){
    print ("Calculating distance to roads...")   
    dist_rd <- terra::distance(terra::rast(map_habitat), map_roads, unit="m", rasterize=TRUE)
    # Convert distance map to raster object
    dist_rd <- raster::raster(dist_rd)
  }
  # Calculate distance to waters 
  # Only perform this section if a min/max distance to waters is specified
  if (exists("map_waters") & !is.null(map_waters) & (water_dist_min>0 | water_dist_max>0)){
    print ("Calculating distance to waters...") 
    dist_wt <- terra::distance(terra::rast(map_habitat), map_waters, unit="m", rasterize=TRUE)
    # Convert distance map to raster object
    dist_wt <- raster::raster(dist_wt)
  }
  
  ###############
  ## 2.4.4. Calculate distances to habitat edge
  # Only run this section if a min/max distance to edge is specified
  if (edge_dist_min > 0 | edge_dist_max > 0){  
    print ("Calculating distance to habitat edge...") 
    dist_ed <- ecosamp_edgedist(map_habitat)
  }
  
  ###############
  ## 2.4.5. Selection of eligible sample region
  # Work with plain vectors of pixel values (one value per pixel of the habitat map), which is
  #   much faster and lighter on memory than repeatedly copying RasterLayer objects.
  #   The road, water and edge distance maps are on the same grid as the habitat map.
  values_habitat <- raster::values(map_habitat)
  # Pixels excluded from the sample design: outside the map, or habitat assigned a value of 0
  excluded <- is.na(values_habitat) | values_habitat == 0
  
  # Exclude pixels too close to (or too far from) the road
  if (exists("map_roads") & !is.null(map_roads)){
    if (road_dist_min > 0){
      excluded[which(raster::values(dist_rd) < road_dist_min)] <- TRUE
    }
    if (road_dist_max > 0){
      excluded[which(raster::values(dist_rd) > road_dist_max)] <- TRUE
    }
  }
  
  # Exclude pixels too close to (or too far from) the water
  if (exists("map_waters") & !is.null(map_waters)){
    if (water_dist_min > 0){
      excluded[which(raster::values(dist_wt) < water_dist_min)] <- TRUE
    }
    if (water_dist_max > 0){
      excluded[which(raster::values(dist_wt) > water_dist_max)] <- TRUE
    }
  }
  
  # Exclude pixels too close to (or too far from) the habitat edge
  if (edge_dist_min > 0){
    excluded[which(raster::values(dist_ed) < edge_dist_min)] <- TRUE
  }
  if (edge_dist_max > 0){
    excluded[which(raster::values(dist_ed) > edge_dist_max)] <- TRUE
  }
  
  # The elevation and slope maps are not necessarily on the same grid as the habitat map,
  #   so these rules are applied with RasterLayer indexing, as in previous versions
  if (elevation_min > 0 | elevation_max > 0 | 
      (exists("map_slope") & !is.null(map_slope) & (slope_min > 0 | slope_max > 0))){
    map_mask <- map_habitat
    raster::values(map_mask) <- 1
    # Pixels with too low (or too high) elevation
    if (elevation_min > 0){
      map_mask[map_elevation < elevation_min] <- 0
    }
    if (elevation_max > 0){
      map_mask[map_elevation > elevation_max] <- 0
    }
    # Pixels with slopes too gentle (or too steep)
    if (exists("map_slope") & !is.null(map_slope)){
      if (slope_min > 0){
        map_mask[map_slope < slope_min] <- 0
      }
      if (slope_max > 0){
        map_mask[map_slope > slope_max] <- 0
      }
    }
    excluded[raster::values(map_mask) == 0] <- TRUE
    rm(map_mask)
  }
  
  # Exclude pixels of treatments assigned a value of 0
  if (exists("map_treatmt") & !is.null(map_treatmt)){
    values_treatmt <- raster::values(map_treatmt)
    excluded[which(values_treatmt == 0)] <- TRUE
  }
  
  # Eligible pixels, in increasing order of cell number, with their coordinates
  cand_cell <- which(!excluded)
  cand_xy <- raster::xyFromCell(map_habitat, cand_cell)
  cand_x <- cand_xy[, 1]
  cand_y <- cand_xy[, 2]
  rm(excluded, cand_xy)
  
  ###############
  ## 2.4.6. Acquire all bins (combinations of habitat*treatment) on map
  # Unique combinations of habitat*treatment among eligible pixels, in order of first appearance
  if (exists("map_treatmt") & !is.null(map_treatmt)){
    landscapes <- data.frame(Treatment = values_treatmt[cand_cell],
                             Habitat = values_habitat[cand_cell])
  } else {
    landscapes <- data.frame(Habitat = values_habitat[cand_cell])
  }
  # Bin of each eligible pixel (NA if it has no treatment)
  cand_key <- do.call(paste, landscapes)
  landscapes <- unique(landscapes)
  # Remove NA values
  landscapes <- stats::na.omit(landscapes)
  cand_bin <- match(cand_key, do.call(paste, landscapes))
  rm(cand_key)
  
  ###############
  ## 2.4.7. Generation of sample points
  ## 2.4.7.1. Creating an empty data frame for storing sample points and summarising
  if (exists("map_treatmt") & !is.null(map_treatmt)){
    # Storage: contains columns for XY axis, treatment and habitat type of sample point
    sd_points <- data.frame(x=numeric(),y=numeric(),treatment=numeric(),habitat=numeric())
    # Summarising: count the number of points for each combination
    sd_pointcount <- data.frame(habitat=landscapes$Habitat,treatment=landscapes$Treatment,points=0)
  } else {
    sd_points <- data.frame(x=numeric(), y=numeric(), habitat=numeric())
    sd_pointcount <- data.frame(habitat=landscapes$Habitat, points=0)
  }
  ## 2.4.7.2. Loop to generate sample points
  # Create a copy of combinations, with a bin number and the number of eligible pixels in each bin
  sd_landscape <- landscapes
  sd_landscape$bin <- seq_len(nrow(landscapes))
  sd_landscape$n_pixels <- tabulate(cand_bin, nbins = nrow(landscapes))
  
  # If no minimum number of points per bin is required, skip the bin-filling
  #   stage, so all points are sampled from the whole eligible region
  if (point_number_min <= 0){
    sd_landscape <- sd_landscape[0, ]
  }
  
  # For max-distance sampling: distance from each eligible pixel to its nearest existing point
  cand_mindist <- rep(Inf, length(cand_cell))
  # Name of the value column returned by raster::sampleRandom() in previous versions
  layer_name <- names(map_habitat)
  
  # Select one of the eligible pixels in idx (positions in cand_cell). This makes the same single
  #   random draw over the same pixels, in the same order, as raster::sampleRandom(size = 1) and
  #   ecosamp_maxdist_select(), so results are identical to previous versions for the same seed.
  select_cell <- function(idx){
    if (max_distance_sampling == FALSE){
      if (length(idx) == 1) {
        return(cand_cell[idx])
      }
      return(cand_cell[idx[ceiling(stats::runif(1) * length(idx))]])
    }
    if (nrow(sd_points) == 0) {
      return(sample(cand_cell[idx], size = 1))
    }
    return(cand_cell[idx[which.max(cand_mindist[idx])]])
  }
  
  # Record a new sample point at pixel cell_new
  #   (value: the value of the selected pixel, reported in column "layer" by previous versions)
  add_point <- function(cell_new, value, treatment_first){
    xy_new <- raster::xyFromCell(map_habitat, cell_new)
    site <- data.frame(x = xy_new[1], y = xy_new[2])
    if (max_distance_sampling == FALSE){
      site[[layer_name]] <- value
    }
    # Extract the treatment and habitat type of the point
    if (exists("map_treatmt") & !is.null(map_treatmt) & treatment_first){
      site$treatment <- values_treatmt[cell_new]
    }
    site$habitat <- values_habitat[cell_new]
    if (exists("map_treatmt") & !is.null(map_treatmt) & !treatment_first){
      site$treatment <- values_treatmt[cell_new]
    }
    # Distance from each remaining eligible pixel to the new point
    if (point_dist_min > 0 | max_distance_sampling == TRUE){
      dist_new <- sqrt((cand_x - xy_new[1])^2 + (cand_y - xy_new[2])^2)
      cand_mindist <<- pmin(cand_mindist, dist_new)
      # Remove pixels too close to the point (prevent points being generated too close)
      if (point_dist_min > 0){
        keep <- !(dist_new <= point_dist_min)
        cand_cell <<- cand_cell[keep]
        cand_x <<- cand_x[keep]
        cand_y <<- cand_y[keep]
        cand_bin <<- cand_bin[keep]
        cand_mindist <<- cand_mindist[keep]
      }
    }
    return(site)
  }
  
  # A variable to count number of generated points
  count <- 0
  
  print("Generating sample points...")
  repeat{  # Loop start
    if (length(cand_cell) == 0) {
      message("No further valid region on sample map. Sample point generation stopping early.")
      break
    }
    # Count of point number +1
    count <- count+1
    # A variable to select combinations that have no reached the minimum number of points
    sample_id <- NA
    
    # Check if there is any remaining in the list of combinations not reaching minimum point number
    if (nrow(sd_landscape) >= 1){
      if (area_weighting == TRUE){
        # If area_weighting is TRUE, the combination with smallest remaining pixel will be selected
        sample_id <- sd_landscape[which.min(sd_landscape$n_pixels), ]
      } else {
        # If area_weighting is FALSE, randomly select a row in the list and return to sample_id,
        #   if not, sample_id will return NA
        sample_id <- sd_landscape[sample(nrow(sd_landscape), size = 1, replace = FALSE), ]
      }
    }
    
    # If sample_id is not NA, sample from the pixels with the same treatment*habitat combination
    #   as sample_id
    if(!is.na(unlist(sample_id)[1])){
      # Remaining eligible pixels with the same treatment*habitat combination as sample_id
      idx_bin <- which(cand_bin == sample_id$bin)
      # If there's no further valid region for this bin on the map, 
      # skip this bin and remove it from the list
      if (length(idx_bin) == 0) {
        message("Not enough area for this bin, minimum number of points cannot be reached:")
        message(paste("Habitat:",index_habitat$Field[index_habitat$Index==sample_id$Habitat]))
        if (exists("map_treatmt") & !is.null(map_treatmt)){
          message(paste("Treatment:",index_treatmt$Field[index_treatmt$Index==sample_id$Treatment]))
        }
        count <- count - 1
        # Only remove this habitat*treatment combination, not other treatments of the same habitat
        sd_landscape <- sd_landscape[sd_landscape$bin != sample_id$bin, ]
      } else {
        # Select a pixel (randomly, or furthest from existing points) and store the new point
        site <- add_point(select_cell(idx_bin), value = 1, treatment_first = TRUE)
        sd_points <- rbind(sd_points,site)
        
        # Update n_pixels column in sd_landscape
        sd_landscape$n_pixels <- tabulate(cand_bin, nbins = nrow(landscapes))[sd_landscape$bin]
        
        # Add 1 to the count of this combination
        if (exists("map_treatmt") & !is.null(map_treatmt)){
          sd_pointcount$points[sd_pointcount$habitat==site$habitat & sd_pointcount$treatment==site$treatment] <- 
            sd_pointcount$points[sd_pointcount$habitat==site$habitat & sd_pointcount$treatment==site$treatment]+1
          # If minimum number is reached for this combination, remove it from the list
          if (sd_pointcount$points[sd_pointcount$habitat==site$habitat & sd_pointcount$treatment==site$treatment]>=point_number_min){
            sd_landscape <- sd_landscape[!(sd_landscape$Treatment==site$treatment & sd_landscape$Habitat==site$habitat),]
          }
        } else {
          sd_pointcount$points[sd_pointcount$habitat==site$habitat] <- 
            sd_pointcount$points[sd_pointcount$habitat==site$habitat]+1
          # If minimum number is reached for this combination, remove it from the list
          if (sd_pointcount$points[sd_pointcount$habitat==site$habitat]>=point_number_min){
            sd_landscape <- sd_landscape[!(sd_landscape$Habitat==site$habitat),]
          }      
        }
      }
      # If sample_id is NA, sample from all the remaining pixels on the map
    } else{
      # Select a pixel (randomly, or furthest from existing points) and store the new point
      cell_new <- select_cell(seq_along(cand_cell))
      site <- add_point(cell_new, value = values_habitat[cell_new], treatment_first = FALSE)
      sd_points <- rbind(sd_points,site)
      
      # Add 1 to the count of this combination
      if (exists("map_treatmt") & !is.null(map_treatmt)){
        sd_pointcount$points[sd_pointcount$habitat==site$habitat & sd_pointcount$treatment==site$treatment] <- 
          sd_pointcount$points[sd_pointcount$habitat==site$habitat & sd_pointcount$treatment==site$treatment]+1
      } else {
        sd_pointcount$points[sd_pointcount$habitat==site$habitat] <- 
          sd_pointcount$points[sd_pointcount$habitat==site$habitat]+1
      }
    }
    
    # Break the loop when enough sample points have been generated
    if (count >= n_total){
      break
    }
  } # Loop end

  # Report number of sample points generated
  print(paste("Complete! Generated",nrow(sd_points),"sample points."))  
  
  # Remove unnecessary columns in output
  sd_points <- sd_points %>% dplyr::select(-dplyr::any_of("layer"))
  # Add an ID label column
  sd_points$ID <- 1:nrow(sd_points)
  
  ###############
  ## 2.4.8. Report generation outputs
  # Summarise number of points in each habitat*treatment combination
  
  if (exists("map_treatmt") & !is.null(map_treatmt)){
    sd_pointcount <- sd_points %>% dplyr::group_by(habitat,treatment) %>% dplyr::summarise(num=n())
    colnames(sd_pointcount)[length(colnames(sd_pointcount))] <- "Number of points"
    # Convert numeric habitat / treatment index back to text and rename column
    if (exists("index_habitat") & !is.null(index_habitat)){
      sd_pointcount <- merge(sd_pointcount, index_habitat, by.x="habitat",by.y="Index")
      colnames(sd_pointcount)[length(colnames(sd_pointcount))] <- "Habitat"
    }
    if (exists("index_treatmt") & !is.null(index_treatmt)){
      sd_pointcount <- merge(sd_pointcount, index_treatmt, by.x="treatment",by.y="Index")
      colnames(sd_pointcount)[length(colnames(sd_pointcount))] <- "Treatment"
    }
  } else {
    sd_pointcount <- sd_points %>% dplyr::group_by(habitat) %>% dplyr::summarise(num=n())
    colnames(sd_pointcount)[length(colnames(sd_pointcount))] <- "Number of points"
    # Convert numeric habitat / treatment index back to text and rename column
    if (exists("index_habitat") & !is.null(index_habitat)){
      sd_pointcount <- merge(sd_pointcount, index_habitat, by.x="habitat",by.y="Index")
      colnames(sd_pointcount)[length(colnames(sd_pointcount))] <- "Habitat"
    }
  }
  # Report number of points in each type of habitat and/or treatment
  print("Below are the details:")
  print(sd_pointcount %>% dplyr::select(dplyr::any_of("Habitat"),
                                        dplyr::any_of("habitat"),
                                        dplyr::any_of("Treatment"),
                                        dplyr::any_of("treatment"),
                                        dplyr::any_of("Number of points")))

  # Create a copy for assigning spatial feature
  sd_points_sp <- sd_points
  # Convert sd_points_sp into spatial data frame
  sp::coordinates(sd_points_sp) <- ~x+y
  raster::crs(sd_points_sp) <- sd_crs
  
  # Convert numeric habitat and treatment features back to text
  if (exists("index_habitat") & !is.null(index_habitat)){
    sd_points <- sp::merge(sd_points, index_habitat, by.x="habitat",by.y="Index")
    colnames(sd_points)[length(colnames(sd_points))] <- "Habitat"
  }
  if (exists("map_treatmt") & !is.null(map_treatmt) & exists("index_treatmt") & !is.null(index_treatmt)){
    sd_points <- sp::merge(sd_points, index_treatmt, by.x="treatment",by.y="Index")
    colnames(sd_points)[length(colnames(sd_points))] <- "Treatment"
  }
  
  # Extract and report distances to habitat edge, roads and waters, if provided
  if (exists("dist_ed")) {
    sd_points$Dist_edge <- raster::extract(dist_ed, sd_points_sp)
  }
  #print("Minimum/maximum distance to habitat edge: ",min(sd_points$Dist_edge),"m/ ",
  #      max(sd_points$Dist_edge), "m")
  if (exists("dist_rd")){
    sd_points$Dist_road <- raster::extract(dist_rd, sd_points_sp)
    #  print("Minimum/maximum distance to road: ",min(sd_points$Dist_road),"m/ ",
    #        max(sd_points$Dist_road), "m")
  }
  if (exists("dist_wt")){
    sd_points$Dist_water <- raster::extract(dist_wt, sd_points_sp)
    #  print("Minimum/maximum distance to water: ",min(sd_points$Dist_water),"m/ ",
    #        max(sd_points$Dist_water), "m")
  }
  # Extract slope angle, if provided
  if (exists("map_slope") & !is.null(map_slope)){
    sd_points$Slope <- raster::extract(map_slope, sd_points_sp)
    #  print("Minimum/maximum slope angle: ",min(sd_points$Slope),"°/ ",
    #        max(sd_points$Slope), "°")
  }
  
  ###############
  ## 2.4.9. Assign spatial property to data frame and export for output
  # Convert coordinates into latitude-longitude CRS (in degrees)
  sd_points_latlon <- sp::spTransform(sd_points_sp,latlon_crs)  
  
  # Create a data frame for exporting
  points_shape <- sd_points
  # Include columns for latitude and longitude of sample points
  points_shape$Longitude <- sp::coordinates(sd_points_latlon)[,1]
  points_shape$Latitude <- sp::coordinates(sd_points_latlon)[,2]
  # Create coordinates columns
  points_shape$coords <- paste(points_shape$x,points_shape$y,sep=", ")  
  points_shape$latlon_coords <- paste(points_shape$Longitude,points_shape$Latitude,sep=", ")
  
  # If plot_results is TRUE, then plot the habitat map and sample points
  if (plot_results == TRUE) {
    print("Plotting sample points...")
    raster::plot(map_habitat)
    raster::plot(sd_points_sp,add=T)
  }
  
  return(points_shape)
}