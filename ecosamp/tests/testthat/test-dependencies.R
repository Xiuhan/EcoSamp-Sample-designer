test_that("ecosamp() runs with only ecosamp loaded (issue 3.3)", {
  skip_on_cran()
  skip_if_not_installed("callr")
  # Run in a fresh R session where raster, terra and dplyr are NOT attached.
  # Previously this failed with: could not find function "%>%".
  helper <- normalizePath(test_path("helper-maps.R"))
  # Under R CMD check the installed package is used; otherwise load the
  # source tree being tested.
  pkg_path <- if (testthat::is_checking()) NULL else normalizePath(test_path("../.."))
  res <- callr::r(function(helper, pkg_path) {
    if (!is.null(pkg_path)) pkgload::load_all(pkg_path, quiet = TRUE)
    source(helper, local = TRUE)
    set.seed(1)
    out <- run_quiet(n_total = 3, map_habitat = habitat_two(),
                     map_res = 10, sd_crs = "EPSG:32723")
    c(n = nrow(out), dplyr_attached = "package:dplyr" %in% search())
  }, args = list(helper = helper, pkg_path = pkg_path))
  expect_equal(unname(res["dplyr_attached"]), 0)
  expect_equal(unname(res["n"]), 3)
})
