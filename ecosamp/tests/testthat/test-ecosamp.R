test_that("point_number_min = 0 spreads points across the map (issue 3.1)", {
  # Previously every point went to the bin with the least remaining area,
  # so all 20 points ended up in the small Pasture bin.
  set.seed(1)
  pts <- run_quiet(n_total = 20, map_habitat = habitat_two(),
                   map_res = 10, sd_crs = "EPSG:32723", point_number_min = 0)
  expect_equal(nrow(pts), 20)
  # Forest is 75% of the area; under random placement it gets ~15 points.
  expect_gte(sum(pts$Habitat == "Forest"), 8)
})

test_that("point_number_min is still honoured per bin", {
  set.seed(1)
  pts <- run_quiet(n_total = 6, map_habitat = habitat_two(),
                   map_res = 10, sd_crs = "EPSG:32723", point_number_min = 3)
  expect_equal(as.vector(table(pts$Habitat)[c("Forest", "Pasture")]), c(3, 3))
})

test_that("treatment map works without index_treatmt (issue 3.4)", {
  # Previously the default treatment index was only built when index_habitat
  # was NULL, which it never is for an sf habitat map, so this errored.
  treat <- make_map("Treatment", c("Control", "Plantation"),
                    list(make_rect(0, 0, 2000, 500), make_rect(0, 500, 2000, 1000)))
  index_habitat <- data.frame(Field = c("Forest", "Pasture"), Index = c(1, 2))
  set.seed(1)
  expect_no_error(
    pts <- run_quiet(n_total = 4, map_habitat = habitat_two(),
                     index_habitat = index_habitat, map_treatmt = treat,
                     map_res = 10, sd_crs = "EPSG:32723", point_number_min = 1)
  )
  expect_setequal(unique(pts$Treatment), c("Control", "Plantation"))
})

test_that("an exhausted bin only removes its own habitat x treatment (issue 3.5)", {
  # Pasture x Control is a 20 m square: it fits one point, then runs out.
  # Previously that also dropped Pasture x Plantation from the minimum list,
  # so it could end up with fewer than point_number_min points.
  treat <- make_map("Treatment", c("Plantation", "Control"),
                    list(sf::st_difference(make_rect(0, 0, 2000, 1000),
                                           make_rect(1900, 0, 1920, 20)),
                         make_rect(1900, 0, 1920, 20)))
  habitat <- make_map("Habitat", c("Forest", "Pasture"),
                      list(make_rect(0, 0, 1800, 1000), make_rect(1800, 0, 2000, 1000)))
  index_treatmt <- data.frame(Field = c("Plantation", "Control"), Index = c(1, 2))
  for (seed in 1:5) {
    set.seed(seed)
    pts <- run_quiet(n_total = 5, map_habitat = habitat, map_treatmt = treat,
                     index_treatmt = index_treatmt,
                     map_res = 10, sd_crs = "EPSG:32723",
                     point_number_min = 2, point_dist_min = 100)
    n_pasture_plantation <- sum(pts$Habitat == "Pasture" & pts$Treatment == "Plantation")
    expect_gte(n_pasture_plantation, 2)
  }
})

test_that("points respect the minimum spacing", {
  set.seed(1)
  pts <- run_quiet(n_total = 10, map_habitat = habitat_two(),
                   map_res = 10, sd_crs = "EPSG:32723", point_dist_min = 300)
  expect_gte(min(stats::dist(pts[, c("x", "y")])), 300)
})

test_that("max-distance sampling runs and respects the spacing", {
  set.seed(1)
  pts <- run_quiet(n_total = 6, map_habitat = habitat_two(),
                   map_res = 10, sd_crs = "EPSG:32723", point_dist_min = 300,
                   max_distance_sampling = TRUE)
  expect_equal(nrow(pts), 6)
  expect_gte(min(stats::dist(pts[, c("x", "y")])), 300)
})
