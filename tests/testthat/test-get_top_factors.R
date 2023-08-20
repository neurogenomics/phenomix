test_that("get_top_factors works", {
  
    obj <- get_HPO()
    top_factors <- get_top_factors(obj = obj,
                                   term = c("parkinson","cardio"),
                                   search_col = "HPO_label")
    testthat::expect_equal(names(top_factors),"PC_3")
    
    top_factors <- get_top_factors(obj = obj,
                                   term = c("parkinson","cardio"),
                                   keys = "umap",
                                   search_col = "HPO_label")
    testthat::expect_equal(names(top_factors),"UMAP_1")
})
