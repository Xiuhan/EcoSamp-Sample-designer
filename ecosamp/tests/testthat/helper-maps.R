# Small synthetic maps (UTM 23S, metres) so tests run in seconds and need no
# external data. Each helper returns an sf object with one attribute column.

make_rect <- function(xmin, ymin, xmax, ymax) {
  sf::st_polygon(list(rbind(c(xmin, ymin), c(xmax, ymin), c(xmax, ymax),
                            c(xmin, ymax), c(xmin, ymin))))
}

make_map <- function(field, values, rects, crs = 32723) {
  df <- data.frame(values)
  names(df) <- field
  sf::st_sf(df, geometry = sf::st_sfc(rects, crs = crs))
}

# 2000 m x 1000 m site: Forest covers 75% of it, Pasture 25%.
habitat_two <- function() {
  make_map("Habitat", c("Forest", "Pasture"),
           list(make_rect(0, 0, 1500, 1000), make_rect(1500, 0, 2000, 1000)))
}

# Run ecosamp() quietly (it prints progress) and return its output.
run_quiet <- function(...) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(ecosamp::ecosamp(...)))
  out
}

# Four habitats (one excluded, like a road) and two treatments, used by the
# sampling tests and by fixtures/make-reference-points.R.
habitat_four <- function() {
  make_map("Habitat", c("Forest", "Pasture", "Scrub", "Road"),
           list(make_rect(0, 0, 1200, 1000), make_rect(1200, 0, 1800, 1000),
                make_rect(1800, 0, 2000, 300), make_rect(1800, 300, 2000, 1000)))
}
index_habitat_four <- function() {
  data.frame(Field = c("Forest", "Pasture", "Scrub", "Road"), Index = c(1, 2, 3, 0))
}
treatment_two <- function() {
  make_map("Treatment", c("Control", "Plantation"),
           list(make_rect(0, 0, 2000, 400), make_rect(0, 400, 2000, 1000)))
}
index_treatment_two <- function() {
  data.frame(Field = c("Control", "Plantation"), Index = c(1, 2))
}

# Settings covered by the reference-points test.
reference_settings <- function() {
  data.frame(
    treatment        = c(FALSE, TRUE,  FALSE, TRUE,  FALSE, TRUE,  TRUE,  TRUE),
    max_distance     = c(FALSE, FALSE, TRUE,  TRUE,  FALSE, FALSE, TRUE,  FALSE),
    area_weighting   = c(TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  FALSE, FALSE, TRUE),
    point_number_min = c(2,     2,     2,     2,     0,     2,     0,     3),
    point_dist_min   = c(250,   250,   250,   250,   250,   250,   0,     400),
    edge_dist_min    = c(50,    50,    50,    50,    0,     0,     50,    0),
    seed             = c(1,     2,     3,     4,     5,     6,     7,     8)
  )
}

# Run ecosamp() on habitat_four() with row k of reference_settings().
run_reference_setting <- function(s) {
  args <- list(n_total = 12, map_habitat = habitat_four(),
               index_habitat = index_habitat_four(),
               map_res = 10, sd_crs = "EPSG:32723",
               max_distance_sampling = s$max_distance,
               area_weighting = s$area_weighting,
               point_number_min = s$point_number_min,
               point_dist_min = s$point_dist_min,
               edge_dist_min = s$edge_dist_min)
  if (s$treatment) {
    args <- c(args, list(map_treatmt = treatment_two(),
                         index_treatmt = index_treatment_two()))
  }
  set.seed(s$seed)
  do.call(run_quiet, args)
}
