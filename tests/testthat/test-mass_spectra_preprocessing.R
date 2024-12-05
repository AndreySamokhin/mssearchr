test_that("PreprocessMassSpectra(), binning", {
  mz <- seq(51L, 55L)
  intst <- c(seq(200, 800, 200), 999)
  bin_bndry <- 0.649
  msp_objs <- list(
    list(name = "ms1",
         mz = rep(mz, each = 2) + bin_bndry - c(0.999, 0.001),
         intst = rep(intst / 2, each = 2) + c(-50, +50)),
    # 50.650 51.648 51.650 52.648 52.650 53.648 53.650 54.648 54.650 55.648
    # 50.0 150.0 150.0 250.0 250.0 350.0 350.0 450.0 449.5 549.5
    list(name = "ms2",
         mz = rep(mz, 2) + bin_bndry - rep(c(0.999, 0.001), each = length(mz)),
         intst = rep(intst / 2, 2) + rep(c(-50, +50), each = length(intst)))
    # 50.650 51.650 52.650 53.650 54.650 51.648 52.648 53.648 54.648 55.648
    # 50.0 150.0 250.0 350.0 449.5 150.0 250.0 350.0 450.0 549.5
  )
  pp_msp_objs <- PreprocessMassSpectra(msp_objs, bin_boundary = bin_bndry)
  for (x in pp_msp_objs) {
    expect_equal(x$mz, mz, tolerance = testthat_tolerance())
    expect_equal(x$intst, intst, tolerance = testthat_tolerance())
  }
})



test_that("PreprocessMassSpectra(), intensity normalization", {
  mz <- seq(51L, 55L)
  intst <- c(seq(200, 800, 200), 999)
  msp_objs <- list(
    list(name = "ms1", mz = mz, intst = intst / 999),
    list(name = "ms2", mz = mz, intst = 150 * intst)
  )
  pp_msp_objs <- PreprocessMassSpectra(msp_objs)
  for (x in pp_msp_objs) {
    expect_equal(x$intst, intst, tolerance = testthat_tolerance())
  }
})



test_that("PreprocessMassSpectra(), zero-intensity peaks", {
  mz <- seq(52L, 60L, 2L)
  intst <- c(seq(200, 800, 200), 999)
  msp_objs <- list(
    list(name = "ms1",
         mz = c(rbind(mz, seq(53L, 61L, 2L))),
         intst = c(rbind(intst, rep(0, length(intst))))),
    list(name = "ms2",
         mz = c(rbind(mz, seq(53L, 61L, 2L))),
         intst = c(rbind(intst, rep(0.49, length(intst)))))
  )
  pp_msp_objs <- PreprocessMassSpectra(msp_objs, remove_zeros = TRUE)
  for (x in pp_msp_objs) {
    expect_equal(x$mz, mz, tolerance = testthat_tolerance())
    expect_equal(x$intst, intst, tolerance = testthat_tolerance())
  }
})



test_that("PreprocessMassSpectra(), the 'preprocessed' attribute", {
  mz <- seq(51L, 55L)
  intst <- c(seq(200, 800, 200), 999)
  msp_objs <- list(list(name = "ms1", mz = mz, intst = intst))
  pp_msp_objs <- PreprocessMassSpectra(msp_objs)
  expect_true(attr(pp_msp_objs[[1]], "preprocessed"))
})


