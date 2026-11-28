test_that("sample points match the reference outputs", {
  # reference-points.rds holds the outputs of the previous sampling code for
  # these settings (see fixtures/make-reference-points.R). The current code
  # must select exactly the same points for the same seed.
  # Longitude and latitude come from sp::spTransform(), whose last digits
  # depend on the PROJ library version, so they are compared with a small
  # tolerance; every other column (UTM coordinates, habitat, treatment,
  # distances) must be identical.
  ref <- readRDS(test_path("fixtures", "reference-points.rds"))
  settings <- reference_settings()
  expect_equal(ref$settings, settings)
  latlon <- c("Longitude", "Latitude", "latlon_coords")
  for (k in seq_len(nrow(settings))) {
    actual <- run_reference_setting(settings[k, ])
    expected <- ref$points[[k]]
    expect_identical(actual[setdiff(names(actual), latlon)],
                     expected[setdiff(names(expected), latlon)],
                     label = paste("setting", k))
    expect_equal(actual[c("Longitude", "Latitude")],
                 expected[c("Longitude", "Latitude")], tolerance = 1e-9,
                 label = paste("setting", k, "longitude/latitude"))
  }
})

test_that("the same seed gives the same points", {
  s <- reference_settings()[2, ]
  expect_identical(run_reference_setting(s), run_reference_setting(s))
})

test_that("sampling stops early when no eligible area is left", {
  set.seed(1)
  out <- NULL
  utils::capture.output(
    expect_message(
      out <- ecosamp(n_total = 50, map_habitat = habitat_two(), map_res = 10,
                     sd_crs = "EPSG:32723", point_dist_min = 600),
      "No further valid region"
    )
  )
  expect_lt(nrow(out), 50)
  expect_gte(min(stats::dist(out[, c("x", "y")])), 600)
})

test_that("max-distance sampling spreads points further than random sampling", {
  min_spacing <- function(max_distance) {
    set.seed(1)
    pts <- run_quiet(n_total = 8, map_habitat = habitat_two(), map_res = 10,
                     sd_crs = "EPSG:32723", max_distance_sampling = max_distance)
    min(stats::dist(pts[, c("x", "y")]))
  }
  expect_gt(min_spacing(TRUE), min_spacing(FALSE))
})
