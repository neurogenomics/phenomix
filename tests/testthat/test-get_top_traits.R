test_that("get_top_traits works", {
  
    obj <- get_HPO()
    top_phenos <- get_top_traits(obj = obj)
    testthat::expect_equal(nrow(top_phenos$data),171)
    testthat::expect_true(methods::is(top_phenos$data,"data.table"))
    testthat::expect_true(methods::is(top_phenos$plot,"gg"))
})
