test_that("Search(), compare with snapshot (synthetic mass spectra)", {
  load(test_path("data/library_search/library_search_synthetic.rda"))
  # The file contains two objects: 'input_args' and 'snapshot'
  res <- lapply(input_args, function(a1) { do.call("LibrarySearch", a1)})
  for (i in seq_along(res)) {
    for (j in seq_along(res[[i]])) {
      expect_equal(res[[i]][[j]], snapshot[[i]][[j]])
    }
  }
})


