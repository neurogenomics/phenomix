test_that("plot_trait_cor works", {
  
    obj <- get_HPO()[seq(100),]
    knn <- find_neighbors(obj = obj,
                          var1_search = "cardio",
                          label_col = "HPO_label")
    gg_cor <- plot_trait_cor(knn=knn,  top_n = 3)
    testthat::expect_true(methods::is(gg_cor,"plotly"))
})
