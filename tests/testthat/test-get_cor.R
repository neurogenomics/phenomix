test_that("get_cor works", {
  
    obj <- get_HPO()[seq(100),]
    obj2 <- get_cor(obj = obj,
                    return_obj = TRUE,
                    keys = "pca")
    testthat::expect_equal(dim(obj2@graphs$pca_cor),rep(ncol(obj),2))
    
    cmat <- get_cor(obj = obj,
                    return_obj = FALSE,
                    keys = "pca")
    testthat::expect_true(methods::is(cmat,"Graph"))
    testthat::expect_equal(dim(cmat),rep(ncol(obj),2))
})
