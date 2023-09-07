test_that("run_autoencoder works", {
  
    
    obj <- get_HPO()[seq(100),seq(50)] 
    #### 1 hidden layer with 2 nodes ####
    ae_res <- run_autoencoder(obj = obj,
                              hidden = 3,
                              color_var = "group_depth3")
    ##### 3 hidden layers with 164/32/164 nodes #####
    ae_res <- run_autoencoder(obj = obj, 
                              hidden = c(164,32,164),
                              color_var = "group_depth3")
    X <- scKirby::get_obsm(ae_res)[[1]]
    um <- run_umap(obj = X)
    
    
})
