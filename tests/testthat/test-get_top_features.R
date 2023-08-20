test_that("get_top_features works", {
  
    obj <- get_HPO()
    top_features <- get_top_features(obj = obj)
    testthat::expect_equal(nrow(top_features$data), 150)
    testthat::expect_true(methods::is(top_features$data,"data.table"))
    testthat::expect_true(methods::is(top_features$plot,"gg"))
})
