# Generates reference-points.rds: the output of ecosamp() for each row of
# reference_settings() (see helper-maps.R).
#
# The file was generated with the sampling code from before it was rewritten
# to work on vectors of eligible pixels (commit 3ff15ae). The rewrite makes the
# same random draws in the same order, so it must reproduce these outputs
# exactly; test-sampling.R checks this.
#
# Only regenerate this file if a change to ecosamp() is meant to change which
# points are selected. Run from the package folder:
#   pkgload::load_all(); source("tests/testthat/helper-maps.R")
#   source("tests/testthat/fixtures/make-reference-points.R")
settings <- reference_settings()
reference <- lapply(seq_len(nrow(settings)), function(k) run_reference_setting(settings[k, ]))
saveRDS(list(settings = settings, points = reference),
        "tests/testthat/fixtures/reference-points.rds", version = 2)
