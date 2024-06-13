test_that("Search(), compare with snapshot (synthetic test)", {
  load(test_path("data/library_search_synthetic.rda"))
  # The file contains: 'fun_name', 'input_args', 'snapshot'
  res <- do.call(fun_name, input_args)
  for (i in seq_along(res)) {
    expect_equal(res[[i]], snapshot[[i]])
  }
})



test_that("Search(), compare with snapshot (alkanes)", {
  load(test_path("data/library_search_alkanes.rda"))
  # The file contains: 'fun_name', 'input_ms', 'input_optns', 'snapshot'
  for(i in seq_along(input_optns)) {
    input_args <- c(input_ms, input_optns[[i]])
    res <- do.call(fun_name, input_args)
    for (j in seq_along(res)) {
      expect_equal(res[[j]], snapshot[[i]][[j]])
    }
  }
})


