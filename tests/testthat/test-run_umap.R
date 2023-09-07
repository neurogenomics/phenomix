test_that("run_umap works", {
  
    
    obj <- get_HPO()[seq(50),seq(100)]
    um <- run_umap(obj)
    testthat::expect_true(methods::is(um$embedding,"matrix"))
    testthat::expect_equal(nrow(um$embedding), ncol(obj))
    testthat::expect_equal(dim(um$nn$euclidean$dist),
                           c(ncol(obj),formals(run_umap)$n_neighbors)
                           )
})
