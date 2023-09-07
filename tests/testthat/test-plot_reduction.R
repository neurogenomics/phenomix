test_that("plot_reduction works", {
  
    obj <- get_HPO()
    gp <- plot_reduction(obj = obj)
    testthat::expect_true(methods::is(gp,"gg"))
    
    gp <- plot_reduction(obj = obj,
                         keys = "umap", 
                         color_var = "seurat_clusters")
    testthat::expect_true(methods::is(gp,"gg"))
})
