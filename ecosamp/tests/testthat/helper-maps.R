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
