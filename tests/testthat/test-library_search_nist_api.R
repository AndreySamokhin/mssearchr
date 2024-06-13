test_that(".ParseSrcreslt(), compare with snapshot", {
  hitlists <- .ParseSrcreslt(test_path("data/SRCRESLT.TXT"))
  load(test_path("data/srcreslt.rda"))
  # The file contains: 'snapshot'
  for (i in seq_along(hitlists)) {
    expect_equal(hitlists[[i]], snapshot[[i]])
  }
})


