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
#' @param map_slope (optional) RasterLayer or SpatRaster, slope angle map of sample region.
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
#' @param slope_max (optional) Numeric, maximum slope angle to allocate sample points (degrees.)
#' Requires valid input for param map_slope. If left as 0, no slope angle will be considered. 0 by default.
#' 
#' @param max_distance (optional) Logical. If TRUE, sample points will be allocated at the eligible pixel furthest away from existing points;
#' if FALSE, sample points will be randomly allocated within eligible regions.
#' FALSE by default.
#' @param plot_results (optional) Logical. If TRUE, a map of output sample points will be plotted on back ground of habitat type map. 
#' FALSE by default.
#' 
#' @return A data frame containing sample point IDs, coordinates in specified CRSs, 
#' type of habitat (and treatment, distance to edge, distance to road, distance to water, 
#' and slope angle if specified) of generated sample points.
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
  slope_max = 0, map_res = 5,
  
  latlon_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
  map_treatmt = NULL,
  map_roads = NULL,map_waters = NULL,map_slope = NULL,
  index_habitat = NULL,index_treatmt = NULL,
  max_distance = FALSE,
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
    if (is.null(index_habitat)){
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
  # If slope map exists and is an SF object, then reproject its CRS
  # If it exists and is not an SF object, return an error message
  if (exists("map_slope") & !is.null(map_slope) & (slope_max>0)){
    print ("Optional input: slope map detected. Processing...")
    if (is(map_slope, "SpatRaster") | is(map_slope, "RasterLayer")) {
      map_slope <- ecosamp_transcrs(map = map_slope, sd_crs = sd_crs, 
                                    map_res = map_res, ref_map = map_habitat)
    } else {
      stop("Input slope map is not RasterLayer or SpatRaster object")
    }
  } else if ((slope_max>0) & is.null(map_slope)) {
    stop("Maximum slope angle specified but no slope map provided. Please provide the slope map.")
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
    dist_rd <- terra::distance(terra::rast(map_habitat), map_roads, unit="m", rasterize=TRUE, haversine=TRUE)
    # Convert distance map to raster object
    dist_rd <- raster::raster(dist_rd)
  }
  # Calculate distance to waters 
  # Only perform this section if a min/max distance to waters is specified
  if (exists("map_waters") & !is.null(map_waters) & (water_dist_min>0 | water_dist_max>0)){
    print ("Calculating distance to waters...") 
    dist_wt <- terra::distance(terra::rast(map_habitat), map_waters, unit="m", rasterize=TRUE, haversine=TRUE)
    # Convert distance map to raster object
    dist_wt <- raster::raster(dist_wt)
  }
  
  ###############
  ## 2.4.4. Calculate distances to habitat edge
  # Only run this section if a min/max distance to edge is specified
  if (edge_dist_min > 0 | edge_dist_max > 0){  
    # Create an empty RasterLayer object of the region for storing data
    sites_edge <- map_habitat
    values(sites_edge) <- 0
    
    # Loop through each type of habitat in the map to calculate distance to edge
    print ("Calculating distance to habitat edge...") 
    for (i in raster::unique(map_habitat)) {
      # Temporal storage
      edge_diff <- map_habitat
      # Assign a value to out-of-boundary pixels in the map
      edge_diff[is.na(edge_diff)] <- 99999
      # Create an SF object for all other habitats ("foreign habitats") except the current one in the loop
      #   Empty pixels with the same type of habitat
      edge_diff[edge_diff==i] <- NA
      #   Assign an uniform value to all "foreign habitats"
      edge_diff[edge_diff!=i] <- 1
      #   Convert the "foreign habitats" into a SpatRaster object, then into SF polygon object
      edge_diff <- terra::rast(edge_diff)
      edge_diff <- terra::as.polygons(edge_diff)
      # Calculate the distance from each pixel on the map to the nearest "foreign habitat" SF
      dist_tmp <- terra::distance(terra::rast(map_habitat), edge_diff, unit="m", rasterize=TRUE, haversine=TRUE)
      # Convert into a RasterLayer object
      dist_tmp <- raster::raster(dist_tmp)
      # Only keep distance values of pixels with habitat in the current loop
      dist_tmp[map_habitat!=i] <- 0
      
      # Sum distances of current loop and previous loops
      store_tmp <- raster::stack(sites_edge,dist_tmp)
      sites_edge <- raster::calc(store_tmp,sum)
    }
    # Copy final output
    dist_ed <- sites_edge
  }
  
  ###############
  ## 2.4.5. Selection of eligible sample region
  # Make a copy of habitat map
  map_temp_habitat <- map_habitat
  
  # Assign 0 to pixels too close to the road
  if (exists("map_roads") & !is.null(map_roads)){
    if (road_dist_min > 0){
      map_temp_habitat[dist_rd < road_dist_min] <- 0
    }
    # If road_dist_max > 0, assign 0 to pixels too far from the road
    if (road_dist_max > 0){
      map_temp_habitat[dist_rd > road_dist_max] <- 0
    }
  }
  
  # Assign 0 to pixels too close to the water
  if (exists("map_waters") & !is.null(map_waters)){
    if (water_dist_min > 0){
      map_temp_habitat[dist_wt < water_dist_min] <- 0
    }
    # If water_dist_max > 0, assign 0 to pixels too far from the water
    if (water_dist_max > 0){
      map_temp_habitat[dist_wt > water_dist_max] <- 0
    }
  }
  
  # Assign 0 to pixels too close to the habitat edge
  if (edge_dist_min > 0){
    map_temp_habitat[dist_ed < edge_dist_min] <- 0
  }
  # If edge_dist_max > 0, assign 0 to pixels too far from the habitat edge
  if (edge_dist_max > 0){
    map_temp_habitat[dist_ed > edge_dist_max] <- 0
  }
  
  # If slope_max > 0, assign 0 to pixels with slopes too steep
  if (exists("map_slope") & !is.null(map_slope)){
    if (slope_max > 0){
      map_temp_habitat[map_slope > slope_max] <- 0
    }
  }
  
  # Remove regions assigned value of 0 from valid sample map
  map_temp_habitat[map_temp_habitat==0] <- NA
  if (exists("map_treatmt") & !is.null(map_treatmt)){
    map_temp_treatmt <- map_treatmt
    map_temp_habitat[map_temp_treatmt==0] <- NA
    map_temp_treatmt[is.na(map_temp_habitat)] <- NA
  }
  
  ###############
  ## 2.4.6. Acquire all bins (combinations of habitat*treatment) on map
  # Convert raster maps into data frames to get values
  map_temp_habitat_data <- raster::as.data.frame(map_temp_habitat)
  # Rename columns
  colnames(map_temp_habitat_data) <- "Habitat"
  # Repeat for treatment map if it exists
  if (exists("map_treatmt") & !is.null(map_treatmt)){
    map_temp_treatmt_data <-  raster::as.data.frame(map_temp_treatmt)
    colnames(map_temp_treatmt_data) <-  "Treatment"
    # Combine the treatment and habitat data frames
    landscapes <- cbind(map_temp_treatmt_data,map_temp_habitat_data)
  } else {
    # If treatment map doesn't exist, then just use habitat map
    landscapes <- map_temp_habitat_data
  }
  # Get all unique combinations of habitat*treatment
  landscapes <- unique(landscapes)
  # Remove NA values
  landscapes <- stats::na.omit(landscapes)
  
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
  # Create a copy of combinations
  sd_landscape <- landscapes
  # A variable to count number of generated points
  count <- 0
  
  print("Generating sample points...")
  repeat{  # Loop start
    if (!any(!is.na(raster::values(map_temp_habitat)))) {
      message("No further valid region on sample map. Sample point generation stopping early.")
      break
    }
    # Count of point number +1
    count <- count+1
    # A variable to select combinations that have no reached the minimum number of points
    sample_id <- NA
    
    #browser()
    
    # Check if there is any remaining in the list of combinations not reaching minimum point number
    if (nrow(sd_landscape) >= 1){
      # If so, randomly select a row in the list and return to sample_id,
      #   if not, sample_id will return NA
      sample_id <- sd_landscape[sample(nrow(sd_landscape), size = 1, replace = FALSE), ]
      
    }
    
    # If sample_id is not NA, sample from the pixels with the same treatment*habitat combination
    #   as sample_id
    #sample_check <- ifelse((exists("map_treatmt") & !is.null(map_treatmt)),sample_id[1,],sample_id[1])
    if(!is.na(unlist(sample_id)[1])){
      # Create a copy of the base map
      sites_sample <- map_temp_habitat
      # Only keep pixels with the same treatment*habitat combination as sample_id
      #   discard other pixels
      if (exists("map_treatmt") & !is.null(map_treatmt)){
        sites_sample[map_temp_treatmt==sample_id$Treatment & map_temp_habitat==sample_id$Habitat] <- 1
        sites_sample[map_temp_treatmt!=sample_id$Treatment | map_temp_habitat!=sample_id$Habitat] <- NA
      } else {
        sites_sample[map_temp_habitat==sample_id] <- 1
        sites_sample[map_temp_habitat!=sample_id] <- NA
      }
      # If there's no further valid region for this bin on the map, 
      # skip this bin and remove it from the list
      if (!any(!is.na(raster::values(sites_sample)))) {
        message("Not enough area for this bin, minimum number of points cannot be reached:")
        message(paste("Habitat:",index_habitat$Field[index_habitat$Index==sample_id$Habitat]))
        if (exists("map_treatmt") & !is.null(map_treatmt)){
          message(paste("Treatment:",index_treatmt$Field[index_treatmt$Index==sample_id$Treatment]))
        }
        count <- count - 1
        sd_landscape <- sd_landscape[!(sd_landscape$Habitat==sample_id$Habitat),]
      } else {
        # If maxdistance is FALSE, randomly select a point in the remaining pixels
        if (max_distance == FALSE){
          samp <- raster::sampleRandom(x=sites_sample, size = 1, na.rm = TRUE, xy = TRUE)
        } else  {
          samp <- ecosamp_maxdist_select(sites_sample, sd_points)
        }
        
        # Convert output into a data frame and create a backup
        samp <- as.data.frame(samp)
        site <- samp
        
        # Convert samp into a spatial data frame
        sp::coordinates(samp) <- ~ x + y
        # Extract the treatment and habitat type of the point
        if (exists("map_treatmt") & !is.null(map_treatmt)){
          site$treatment <- raster::extract(map_refmap$Treatment,samp)
        }
        site$habitat <- raster::extract(map_refmap$Habitat,samp)
        
        # Combine and store new point with previous points
        sd_points <- rbind(sd_points,site)
        
        # Remove pixels too close to the point generated (prevent points being generated too close)
        #   Calculate distance from each pixel to the new point
        d1 <- raster::distanceFromPoints(map_temp_habitat, samp) 
        #   Remove pixels too close to the point
        map_temp_habitat[d1 <= point_dist_min] <- 0 
        map_temp_habitat[map_temp_habitat==0] <- NA   
        if (exists("map_treatmt") & !is.null(map_treatmt)){
          map_temp_treatmt[map_temp_habitat==0] <- NA
        }
        
        # Add 1 to the count of this combination
        if (exists("map_treatmt") & !is.null(map_treatmt)){
          sd_pointcount$points[sd_pointcount$habitat==site$habitat & sd_pointcount$treatment==site$treatment] <- 
            sd_pointcount$points[sd_pointcount$habitat==site$habitat & sd_pointcount$treatment==site$treatment]+1
          # If minimum number is reached for this combination, remove it from the list
          if (sd_pointcount$points[sd_pointcount$habitat==site$habitat & sd_pointcount$treatment==site$treatment]==point_number_min){
            sd_landscape <- sd_landscape[!(sd_landscape$Treatment==site$treatment & sd_landscape$Habitat==site$habitat),]
          }
        } else {
          sd_pointcount$points[sd_pointcount$habitat==site$habitat] <- 
            sd_pointcount$points[sd_pointcount$habitat==site$habitat]+1
          # If minimum number is reached for this combination, remove it from the list
          if (sd_pointcount$points[sd_pointcount$habitat==site$habitat]==point_number_min){
            sd_landscape <- sd_landscape[!(sd_landscape$Habitat==site$habitat),]
          }      
        }
      }
      # If sample_id is NA, sample from all the remaining pixels on the map
    } else{
      # Create a copy of the base map
      sites_sample <- map_temp_habitat
      # Randomly select a point in the remaining pixels
      # If maxdistance is FALSE, randomly select a point in the remaining pixels
      if (max_distance == FALSE){
        samp <- raster::sampleRandom(x=sites_sample, size = 1, na.rm = TRUE, xy = TRUE)
      } else  {
        samp <- ecosamp_maxdist_select(sites_sample, sd_points)
      }
      
      # Convert output into a data frame and create a backup
      samp <- as.data.frame(samp)
      site <- samp
      
      # Convert samp into a spatial data frame      
      sp::coordinates(samp) <- ~ x + y
      # Extract the habitat (and treatment) type of the point
      site$habitat <- raster::extract(map_refmap$Habitat,samp)
      if (exists("map_treatmt") & !is.null(map_treatmt)){
        site$treatment <- raster::extract(map_refmap$Treatment,samp)
      }
      
      # Combine and store new point with previous points
      sd_points <- rbind(sd_points,site)
      
      # Remove pixels too close to the point generated
      #   Calculate distance from each pixel to the new point
      d1 <- raster::distanceFromPoints(map_temp_habitat, samp) 
      #   Remove pixels too close to the point
      map_temp_habitat[d1<=point_dist_min] <- 0 
      map_temp_habitat[map_temp_habitat==0] <- NA
      if (exists("map_treatmt") & !is.null(map_treatmt)){
        map_temp_treatmt[map_temp_habitat==0] <-  NA
      }
      
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
  
  # Remove unnecessary columns in output
  sd_points <- sd_points %>% dplyr::select(-dplyr::any_of("layer"))
  # Add an ID label column
  sd_points$ID <- 1:nrow(sd_points)
  
  ###############
  ## 2.4.8. Report generation outputs
  # Summarise number of points in each habitat*treatment combination
  
  if (exists("map_treatmt") & !is.null(map_treatmt)){
    sd_pointcount <- sd_points %>% dplyr::group_by(habitat,treatment) %>% dplyr::summarise(num=n())
    # Convert numeric habitat / treatment index back to text and rename column
    sd_pointcount <- merge(sd_pointcount, index_habitat, by.x="habitat",by.y="Index")
    colnames(sd_pointcount)[length(colnames(sd_pointcount))-1] <- "Number of points"
    colnames(sd_pointcount)[length(colnames(sd_pointcount))] <- "Habitat"
    sd_pointcount <- merge(sd_pointcount, index_treatmt, by.x="treatment",by.y="Index")
    colnames(sd_pointcount)[length(colnames(sd_pointcount))] <- "Treatment"
  } else {
    sd_pointcount <- sd_points %>% dplyr::group_by(habitat) %>% dplyr::summarise(num=n())
    # Convert numeric habitat / treatment index back to text and rename column
    sd_pointcount <- merge(sd_pointcount, index_habitat, by.x="habitat",by.y="Index")
    colnames(sd_pointcount)[length(colnames(sd_pointcount))-1] <- "Number of points"
    colnames(sd_pointcount)[length(colnames(sd_pointcount))] <- "Habitat"
  }
  # Report number of sample points generated
  print(paste("Complete! Generated",nrow(sd_points),"sample points."))
  print("Below are the details")
  
  if (exists("map_treatmt") & !is.null(map_treatmt)){
    print(sd_pointcount %>% dplyr::select("Habitat","Treatment","Number of points"))
  } else {
    print(sd_pointcount %>% dplyr::select("Habitat","Number of points"))
  }
  # Create a copy for assigning spatial feature
  sd_points_sp <- sd_points
  # Convert sd_points_sp into spatial data frame
  sp::coordinates(sd_points_sp) <- ~x+y
  crs(sd_points_sp) <- sd_crs
  
  # Convert numeric habitat and treatment features back to text
  sd_points <- sp::merge(sd_points, index_habitat, by.x="habitat",by.y="Index")
  colnames(sd_points)[length(colnames(sd_points))] <- "Habitat"
  if (exists("map_treatmt") & !is.null(map_treatmt)){
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