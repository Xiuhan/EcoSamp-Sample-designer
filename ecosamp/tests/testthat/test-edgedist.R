# The original distance-to-edge code from ecosamp() (before it was moved into
# ecosamp_edgedist()), kept here as a reference to check the new version against.
edge_distance_reference <- function(map_habitat) {
  sites_edge <- map_habitat
  raster::values(sites_edge) <- 0
  for (i in raster::unique(map_habitat)) {
    edge_diff <- map_habitat
    edge_diff[is.na(edge_diff)] <- 99999
    edge_diff[edge_diff == i] <- NA
    edge_diff[edge_diff != i] <- 1
    edge_diff <- terra::as.polygons(terra::rast(edge_diff))
    dist_tmp <- terra::distance(terra::rast(map_habitat), edge_diff,
                                unit = "m", rasterize = TRUE)
    dist_tmp <- raster::raster(dist_tmp)
    dist_tmp[map_habitat != i] <- 0
    sites_edge <- raster::calc(raster::stack(sites_edge, dist_tmp), sum)
  }
  sites_edge
}

# A 1000 m x 600 m map at 10 m resolution with four habitat types, an excluded
# habitat (value 0, like a road) and a hole outside the map (NA).
mixed_habitat_raster <- function() {
  r <- raster::raster(xmn = 0, xmx = 1000, ymn = 0, ymx = 600, res = 10,
                      crs = "EPSG:32723")
  xy <- raster::coordinates(r)
  v <- ifelse(xy[, 1] < 400, 1, ifelse(xy[, 1] < 700, 2, 3))
  v[xy[, 1] > 450 & xy[, 1] < 550 & xy[, 2] > 200 & xy[, 2] < 300] <- 4  # small patch
  v[abs(xy[, 2] - 400) < 10] <- 0                                         # road
  v[xy[, 1] > 850 & xy[, 2] > 450] <- NA                                  # outside the map
  raster::values(r) <- v
  r
}

test_that("ecosamp_edgedist() matches the original edge-distance code", {
  r <- mixed_habitat_raster()
  new <- raster::values(ecosamp_edgedist(r))
  old <- raster::values(edge_distance_reference(r))
  v <- raster::values(r)
  sampled <- !is.na(v) & v != 0       # pixels that can hold sample points
  expect_identical(new[sampled], old[sampled])
  expect_true(all(new[!sampled] == 0))
})

test_that("ecosamp_edgedist() gives the expected distances", {
  # Forest for x < 1500 m, Pasture beyond (10 m pixels): the Forest pixel
  # centred at x = 1495 is 10 m from the nearest Pasture pixel centre (1505),
  # and the one at x = 1405 is 100 m away.
  r <- raster::raster(xmn = 0, xmx = 2000, ymn = 0, ymx = 1000, res = 10,
                      crs = "EPSG:32723")
  raster::values(r) <- ifelse(raster::coordinates(r)[, 1] < 1500, 1, 2)
  d <- ecosamp_edgedist(r)
  at <- function(x, y) raster::extract(d, cbind(x, y))
  expect_equal(at(1495, 505), 10)
  expect_equal(at(1405, 505), 100)
  expect_equal(at(1505, 505), 10)
})

test_that("a habitat with no edge gets an infinite distance", {
  r <- raster::raster(xmn = 0, xmx = 100, ymn = 0, ymx = 100, res = 10,
                      crs = "EPSG:32723")
  raster::values(r) <- 1
  expect_true(all(is.infinite(raster::values(ecosamp_edgedist(r)))))
})

test_that("ecosamp() respects edge_dist_min", {
  set.seed(1)
  pts <- run_quiet(n_total = 6, map_habitat = habitat_two(),
                   map_res = 10, sd_crs = "EPSG:32723", edge_dist_min = 100)
  expect_equal(nrow(pts), 6)
  expect_true(all(pts$Dist_edge >= 100))
})
