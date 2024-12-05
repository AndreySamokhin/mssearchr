test_that("WriteMsp(), write msp files", {
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


