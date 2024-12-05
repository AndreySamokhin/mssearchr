test_that("ReadMsp(), delimiters", {
  msp_objs <- ReadMsp(test_path("data/read_msp/delimiters.msp"))
  for (msp in msp_objs[-1L]) {
    expect_equal(c(msp_objs[[1L]]$mz, msp_objs[[1L]]$intst),
                 c(msp$mz, msp$intst))
  }
})



test_that("ReadMsp(), number format", {
  msp_objs <- ReadMsp(test_path("data/read_msp/number_format.msp"))
  for (msp in msp_objs[-1L]) {
    expect_equal(c(msp_objs[[1L]]$mz, msp_objs[[1L]]$intst),
                 c(msp$mz, msp$intst))
  }
})



test_that("ReadMsp(), duplicated (non-unique) fields are present", {
  expect_error(ReadMsp(test_path("data/read_msp/duplicated_fields.msp")),
               "Field names are not unique")
})



test_that("ReadMsp(), 'mz' field is present", {
  expect_error(ReadMsp(test_path("data/read_msp/mz_field_present.msp")),
               "The 'mz' field is present")
})



test_that("ReadMsp(), 'intst' field is present", {
  expect_error(ReadMsp(test_path("data/read_msp/intst_field_present.msp")),
               "The 'intst' field is present")
})



test_that("ReadMsp(), 'Name' field is missing", {
  expect_error(ReadMsp(test_path("data/read_msp/name_field_missing.msp")),
               "The msp-file is invalid")
})



test_that("ReadMsp(), 'Num Peaks' field is missing", {
  expect_warning(ReadMsp(test_path("data/read_msp/numpeaks_field_missing.msp")),
                 "The 'Num Peaks' field is missing")
})



test_that("ReadMsp(), 'Num Peaks' field is incorrect", {
  expect_warning(ReadMsp(test_path("data/read_msp/numpeaks_field_incorrect.msp")),
                 "The 'Num Peaks' field is incorrect")
})



test_that("ReadMsp(), read single value fields from an msp-file", {
  msp_objs <- ReadMsp(test_path("data/read_msp/read_single_value_fields.msp"))
  field_names <- names(msp_objs[[1]])
  mask <- (field_names %in% c("mz", "intst"))
  expect_identical(unname(unlist(msp_objs[[1]][!mask])),
                   paste0(field_names[!mask], "_value"))
})



test_that("ReadMsp(), read 'Synon' fileds", {
  msp_objs <- ReadMsp(test_path("data/read_msp/merge_synon_fields.msp"))
  expect_identical(unname(unlist(msp_objs[[1]]$synon)),
                   paste0("synon_value_", seq(1L, length(msp_objs[[1]]$synon))))
})



test_that("ReadMsp(), read 'CAS#' and 'NIST#' fileds", {
  msp_objs <- ReadMsp(test_path("data/read_msp/read_cas_nist_fields_1.msp"))
  for (x in msp_objs) {
    expect_identical(x$cas, "cas_no_value")
    expect_identical(x$nist, "nist_no_value")
  }
  msp_objs <- ReadMsp(test_path("data/read_msp/read_cas_nist_fields_2.msp"))
  for (x in msp_objs) {
    expect_identical(x$compound_rep, "compound_rep_value")
    expect_identical(x$nist, "nist_no_value")
  }
})


