run_autoencoder <- function(obj,
                            label_var = "Category"){ 
    # library(dplyr)
    # library(ggplot2)
    # gwas <- phenomix::get_DEGAS()
    # gwas <- Seurat::FindVariableFeatures(gwas, nfeatures = 5000)
    # select_features <- Seurat::VariableFeatures(gwas)
    # #### Subset data ####
    # obj <- gwas# [select_features,] 
    
    
    #### initialize H2O instance ####
    set.seed(2021)
    h2o::h2o.init(max_mem_size = "15g")  
    options("h2o.use.data.table"=TRUE)
    #### Prepare data ####
    features <- h2o::as.h2o(
        Matrix::t(Seurat::GetAssayData(obj,slot = "counts"))
        )
    labels <- h2o::as.h2o(obj[[label_var]]) 
    
    #### Train autoencoder ####
    ae1 <- h2o::h2o.deeplearning(
        x = seq_along(features),
        training_frame = features,
        autoencoder = TRUE,
        hidden = 2,
        activation = 'Tanh',
        sparse = TRUE,
        epochs = 10
    )
    #### Extract the deep features ####
    layer <- 3
    ae1_codings <- h2o::h2o.deepfeatures(object = ae1,
                                         data = features,
                                         layer = 2)
    plot_data <- data.table::data.table(as.matrix(ae1_codings),
                                        obj@meta.data
    )
    
    
    gg_dr <- ggplot(plot_data,
                    aes_string(x=names(ae1_codings)[1], y=names(ae1_codings)[2],
                               color=label_var,
                               label_phe = "label_phe",
                               label_phe_code = "label_phe_code",
                               Number_of_cases = "Number_of_cases",
                               nCount_genes="nCount_genes"
                          )
           ) +
        geom_point(size=1) +
        theme_bw()
    print(gg_dr)
    plotly::ggplotly(gg_dr)
    
}
