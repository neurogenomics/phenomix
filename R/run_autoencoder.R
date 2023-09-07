#' Run autoencoder
#'
#' Run a customisable deep autoencoder to reduce your data to N dimensions,
#' and extract feature importance scores.
#' Uses \link[h2o]{h2o.deeplearning}.
#' @param obs Phenotype metadata.
#' @param seed Seed passed to \link[base]{set.seed} 
#' for reproducibility between runs.
#' @param color_var Variable in \code{obs} to color points by.
#' @param label_var Variable \code{obs} to label points by.
#' @param normalise_method Normalisation method to apply to the data matrix
#'  before training the autoencoder.
#' @inheritParams scKirby::get_x
#' @inheritParams SeuratObject::CreateSeuratObject
#' @inheritParams h2o::h2o.deeplearning 
#' @inheritDotParams h2o::h2o.deeplearning
#' @returns A trained autoencoder.
#' Specifically, a named list containing:
#' \itemize{
#' \item{embedding : }{Latent space embedding from the smallest hidden layer.}
#' \item{model : }{Trained autoencoder model with parameters.}
#' \item{variable_importance : }{Feature importance data from
#'  \link[h2o]{h2o.varimp} in \pkg{data.table} format.}
#' \item{plot_latent : }{2D plot of the embedding. 
#' \emph{Note}: Only uses nodes 1 and 2, 
#' even when the embedding has more dimensions.}
#' \item{plot_importance : }{Feature importance plot from
#'  \link[h2o]{h2o.varimp_plot}.}
#' }
#'
#' @source \href{https://bradleyboehmke.github.io/HOML/autoencoders.html}{
#' autoencoder documentation}
#' 
#' @export
#' @importFrom Matrix t
#' @examples 
#' obj <- get_HPO()[seq(100),seq(50)] 
#' ae_res <- run_autoencoder(obj = obj, color_var = "group_depth3")
run_autoencoder <- function(obj,
                            transpose = TRUE, 
                            normalise_method = "log1p",
                            assay = NULL,
                            slot = NULL,
                            obs = NULL,
                            color_var = NULL,
                            label_var = NULL,
                            hidden = c(2),
                            activation = 'Tanh',
                            sparse = FALSE,
                            variable_importances = TRUE, 
                            epochs = 10,
                            seed = 2020,
                            verbose = TRUE,
                            ...){ 
    # devoptera::args2vars(run_autoencoder, reassign = TRUE)
    requireNamespace("h2o")
    #### Input needs to be in sample (trait) x feature (gene) format ####
    X <- scKirby::get_x(obj = obj,
                        assay = assay,
                        n = 1,
                        transpose = transpose, 
                        slot = slot) 
    if(is.null(obs)) obs <- scKirby::get_obs(obj = obj)
    if(!is.null(normalise_method)){
        X <- normalise(X = X,
                       method = normalise_method)
    }
    #### initialize H2O instance ####
    if (!is.null(seed)) set.seed(seed)
    h2o::h2o.init(max_mem_size = "15g")  
    options("h2o.use.data.table"=TRUE)
    #### Prepare data ####
    ## Need to convert to dense matrix, or else won't keep colnames ####
    features <- h2o::as.h2o(as.matrix(X))  
    # labels <- h2o::as.h2o(obj[[label_var]]) 
    
    #### Train autoencoder ####
    ae1 <- h2o::h2o.deeplearning(
        x = seq_along(features),
        training_frame = features,
        autoencoder = TRUE,
        hidden = hidden,
        activation = activation,
        # distribution = "quantile",
        sparse = sparse,
        variable_importances = variable_importances, 
        epochs = epochs
        # ...
    )
    
    #### Explore variable importance ####
    if(isTRUE(variable_importances)){
        var_importance <- data.table::data.table(h2o::h2o.varimp(ae1))  
        plot_importance <- h2o::h2o.varimp_plot(model = ae1)
        # ### Conduct gene set enrichment analysis with top N genes 
        # gres <- gprofiler2::gost(
        #     query = list(top10=var_importance$variable[1:10],
        #                  top50=var_importance$variable[1:50],
        #                  top100=var_importance$variable[1:100]))
        # gprofiler2::gostplot(gres)
    } else {
        var_importance <- plot_importance <- NULL 
    }
   
    #### Extract the deep features #### 
    ### Figure out which layer is the smallest ###
    compressed_layer <- which(ae1@parameters$hidden==
                                  min(ae1@parameters$hidden))[1]
    ae1_codings <- h2o::h2o.deepfeatures(object = ae1,
                                         data = features,
                                         ## Note: referring to HIDDEN layer
                                         layer = compressed_layer) 
    #### Extract predictions (decoded data) ####
    # ae1_codings <- h2o::h2o.predict(object = ae1, 
    #                                 newdata = features)
    obsm <- as.matrix(ae1_codings) |>
        `row.names<-`(rownames(X)) 
    #### Plot traits in latent space ####
    gg_latent <- plot_reduction(obj = obsm, 
                                x_dim = 1,
                                y_dim = 2,
                                fix_rownames = TRUE,
                                obs = obs,
                                color_var = color_var,
                                label_var = label_var,
                                labels = FALSE)
    ### Reduce further with UMAP ####
    # {
    #     umap_res <- run_umap(mat = latent_mat,
    #                          transpose = FALSE)
    #     gg_umap <- plot_reduction(obj = umap_res,
    #                               obs = obs,
    #                               fix_rownames = TRUE,
    #                               color_var = color_var,
    #                               label_var = label_var,
    #                               labels = FALSE)
    #     plotly::ggplotly(gg_umap)
    # }
    # 
    #### Compare to UMAP run on original matrix ####
    # gg_umap <- plot_reduction(obj,
    #                           reduction = "umap",
    #                           fix_rownames = TRUE,
    #                color_var = color_var,
    #                label_var = label_var,
    #                labels = FALSE)
    ##### Cluster ####
    # # obj <- Seurat::FindClusters(obj, resi)
    # # obj <- scNLP::run_tfidf(object = obj,
    # #                         reduction = "umap", 
    # #                         label_var = label_var,
    # #                         cluster_var = "seurat_clusters")
    # # scNLP::plot_tfidf(object = obj,
    # #                   point_palette = rep(pals::alphabet(),3),
    # #                   label = label_var)
    # plotly::ggplotly(gg_umap)
    # #### PCA ####
    # pca_res <- run_pca(mat = latent_mat, 
    #                    transpose = FALSE)  
    # plot_reduction(obj = pca_res,
    #               obs = obj@meta.data,
    #               color_var = color_var,
    #               label_var = label_var,
    #               labels = FALSE)
    # plot_reduction(obj, reduction = "pca",
    #                color_var = color_var,
    #                label_var = label_var,
    #                labels = FALSE)
    return(list(obsm = obsm,
                model = ae1,
                variable_importance = var_importance,
                plot_latent = gg_latent,
                plot_importance = plot_importance
               ))
}
