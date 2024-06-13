test_that("ReadMsp(), delimiters (m/z values and intensities)", {
  msp_objs <- ReadMsp(test_path("data/ms_delimiters.msp"))
  mz_intst <- c(msp_objs[[1]]$mz, msp_objs[[1]]$intst)
  for (msp in msp_objs) {
    expect_equal(mz_intst, c(msp$mz, msp$intst))
  }
})



test_that("ReadMsp(), number format (m/z values and intensities)", {
  msp_objs <- ReadMsp(test_path("data/ms_number_format.msp"))
  mz_intst <- c(msp_objs[[1]]$mz, msp_objs[[1]]$intst)
  for (msp in msp_objs) {
    expect_equal(mz_intst, c(msp$mz, msp$intst))
  }
})



test_that("ReadMsp(), duplicated (non-unique) fields", {
  expect_error(ReadMsp(test_path("data/nonunique_field.msp")),
               "Field names are not unique")
})



test_that("ReadMsp(), the 'mz' field is present", {
  expect_error(ReadMsp(test_path("data/mz_presence.msp")),
               "The 'mz' field is present")
})



test_that("ReadMsp(), the 'intst' field is present", {
  expect_error(ReadMsp(test_path("data/intst_presence.msp")),
               "The 'intst' field is present")
})



test_that("ReadMsp(), the 'Name' field is absent", {
  expect_error(ReadMsp(test_path("data/name_absence.msp")),
               "The msp-file is invalid")
})



test_that("ReadMsp(), the 'Num Peaks' field is absent", {
  expect_warning(ReadMsp(test_path("data/numpeaks_absence.msp")),
                 "The 'Num Peaks' field is absent")
})



test_that("ReadMsp(), the 'Num Peaks' field is incorrect", {
  expect_warning(ReadMsp(test_path("data/numpeaks_incorrect.msp")),
                 "The 'Num Peaks' field is incorrect")
})



test_that("ReadMsp(), reading ordinary fields from an msp-file", {
  msp_objs <- ReadMsp(test_path("data/reading_ordinary_fields.msp"))
  field_names <- names(msp_objs[[1]])
  mask <- (field_names %in% c("mz", "intst"))
  expect_identical(unname(unlist(msp_objs[[1]][!mask])),
                   paste0(field_names[!mask], "_value"))
})



test_that("ReadMsp(), merging all 'Synon' fileds", {
  msp_objs <- ReadMsp(test_path("data/merging_synon_fields.msp"))
  expect_identical(unname(unlist(msp_objs[[1]]$synon)),
                   paste0("synon_value_", seq(1L, length(msp_objs[[1]]$synon))))
})



test_that("ReadMsp(), reading 'CAS#' and 'NIST#' fileds", {
  msp_objs <- ReadMsp(test_path("data/reading_cas_and_nist.msp"))
  for (x in msp_objs) {
    expect_identical(x$cas, "cas_no_value")
    expect_identical(x$nist, "nist_no_value")
  }
  msp_objs <- ReadMsp(test_path("data/reading_cas_and_nist2.msp"))
  for (x in msp_objs) {
    expect_identical(x$compound_rep, "compound_rep_value")
    expect_identical(x$nist, "nist_no_value")
  }
})



test_that("WriteMsp(), writing msp files", {
  msp_objs <- list(
    list(name = "name_value_1",
         synon = paste0("synon_value_", seq(1L, 3L)),
         db_no = "db_no_value_1",
         ion_mode = "ion_mode_value_1",
         mz = c(1L, 2L),
         intst = c(20L, 999L)),
    list(name = "name_value_2",
         db_no = "db_no_value_2",
         cas_no = "cas_no_value_1",
         mz = c(2L, 4L),
         intst = c(20L, 999L))
  )

  temp_file1 <- tempfile("mssearchr_", fileext = ".txt")
  on.exit(unlink(temp_file1), add = TRUE)
  WriteMsp(msp_objs, temp_file1)
  msp1 <- readLines(temp_file1)
  # cat(paste(paste0("\"", msp1[msp1 != ""], "\""), collapse = ",\n"))
  expected_msp1 <- c("name: name_value_1",
                     "synon: synon_value_1",
                     "synon: synon_value_2",
                     "synon: synon_value_3",
                     "db#: db_no_value_1",
                     "ion_mode: ion_mode_value_1",
                     "num peaks: 2",
                     "1 20",
                     "2 999",
                     "name: name_value_2",
                     "db#: db_no_value_2",
                     "cas#: cas_no_value_1",
                     "num peaks: 2",
                     "2 20",
                     "4 999")
  expect_identical(msp1[msp1 != ""], expected_msp1)

  temp_file2 <- tempfile("mssearchr_", fileext = ".txt")
  on.exit(unlink(temp_file2), add = TRUE)
  WriteMsp(msp_objs, temp_file2, fields = c("db_no", "cas_no"))
  msp2 <- readLines(temp_file2)
  # cat(paste(paste0("\"", msp2[msp2 != ""], "\""), collapse = ",\n"))
  expected_msp2 <- c("name: name_value_1",
                     "db#: db_no_value_1",
                     "num peaks: 2",
                     "1 20",
                     "2 999",
                     "name: name_value_2",
                     "db#: db_no_value_2",
                     "cas#: cas_no_value_1",
                     "num peaks: 2",
                     "2 20",
                     "4 999")
  expect_identical(msp2[msp2 != ""], expected_msp2)
})


