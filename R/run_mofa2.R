#' Run MOFA2
#'
#' Run Multi-Omics Factor Analysis v2 (MOFA2).
#'
#' Uses \link[MOFA2]{run_mofa} to reduce a single dataset,
#' or multiple views of the same dataset, into a
#' set number of factors.
#'
#' @param obj Named list of input data matrices, or a \pkg{Seurat} object.
#' @param obs Metadata associated with \code{obj}. 
#' Will be extracted automatically if \code{obj} is a \pkg{Seurat} object.
#' @param assay Assay to use in \code{obj}.
#' @param slot Slot to use in \code{obj}.
#' @param features Features to use in reduction. Can be:
#' \itemize{
#' \item{\<vector\> : }{A character vector of feature names.}
#' \item{"variable_features" : }{Will trigger the automatic extraction of
#'  precalculated variable features
#'   (via \link[scKirby]{get_variable_features}).}
#' \item{\code{NULL} : }{Uses all features.}
#' } 
#' @param obs_idcol The name of the column to add to \code{obs} 
#' after fixing IDs.
#' @param transpose First transpose each matrix in \code{mat_list}.
#' @param maxiter Maximum number of training iterations (\emph{DEFAULT} = 1000).
#' @param reductions In addition to returning MOFA factors,
#' further reduce the data by "umap" and/or "tsne" to the MOFA factors.
#' @param seed Seed passed to \link[base]{set.seed}
#' for reproducibility between runs.
#' @param verbose Print messages.
#' @inheritParams MOFA2::run_mofa
#' @inheritParams MOFA2::create_mofa
#' @inheritParams MOFA2::prepare_mofa
#' @inheritDotParams MOFA2::run_mofa
#'
#' @source \href{https://biofam.github.io/MOFA2/}{MOFA2 site}
#'
#' @export
#' @import scKirby
#' @importFrom Matrix t 
#' @examples
#' obj <- get_HPO()[seq(100),]
#' res <- run_mofa2(obj, maxiter=10) 
run_mofa2 <- function(obj,
                      groups = NULL,
                      obs = NULL,
                      assay = NULL,
                      slot = NULL,
                      features = "variable_features",
                      obs_idcol = "sample",
                      transpose = FALSE,
                      maxiter = 1000,
                      data_options = NULL,
                      model_options = NULL,
                      training_options = NULL,
                      stochastic_options = NULL,
                      mefisto_options = NULL,
                      outfile = tempfile(fileext = ".mofa_model.hdf5"),
                      reductions = c("umap"),
                      use_basilisk = TRUE,
                      seed = 2020,
                      verbose = TRUE,
                      ...) {
    requireNamespace("MOFA2")
    
    if (!is.null(seed)) set.seed(seed)
    #### Get assays ####
    if(is.null(assay) && 
       scKirby::is_class(obj,"seurat")){ 
        assay <- names(obj@assays)
    }
    #### Transpose matrices ####
    if (isTRUE(transpose) &&
        scKirby::is_class(obj,"matrix_list")) {
        for (i in seq(length(obj))) {
            obj[[i]] <- Matrix::t(obj[[i]])
        }
    }
    #### Fix colnames ####
    ## MOFA rules:
    # Can't be non-ASCI
    # Can't be >50 characters
    obj <- fix_colnames(
        obj = obj,
        width = 40,
        make_unique = TRUE
    )
    #### Prepare groups ####
    if (is.null(obs)){
        obs <- scKirby::get_obs(obj = obj, 
                                verbose = verbose)
    }
    groups <- if ((!is.null(groups)) & length(groups) == 1) obs[[groups]]
    #### Get features ####
    features <- run_mofa2_features(obj = obj,
                                   features = features,
                                   assay = assay)
    #### Find intersecting traits ####
    if(is.list(obj)){
        traits <- Reduce(base::intersect,lapply(obj,scKirby::get_obs_names)) 
        all_traits <- unique(mapply(obj,FUN=colnames))
        messager(formatC(length(traits),big.mark = ","),
                 "/",
                 formatC(length(all_traits),big.mark = ","),
                 "intersecting traits will be used.",v=verbose)
        obj <- lapply(obj,function(x){x[,traits]})
        obs <- obs[traits,]
    } 
    #### Create MOFA object ####
    
    MOFAobject <- MOFA2::create_mofa(
        data = obj,
        groups = groups,
        assay = assay,
        features = features,
        slot = slot
    )
    #### Model parameters ####
    if (is.null(data_options)) {
        data_options <- MOFA2::get_default_data_options(MOFAobject)
    }
    if (is.null(model_options)) {
        model_options <- MOFA2::get_default_model_options(MOFAobject)
    }
    if (is.null(training_options)) {
        training_options <- MOFA2::get_default_training_options(MOFAobject)
        training_options$seed <- seed
        training_options$maxiter <- maxiter
    }
    #### Prepare MOFA object ####
    MOFAobject <- MOFA2::prepare_mofa(
        object = MOFAobject,
        data_options = data_options,
        model_options = model_options,
        training_options = training_options,
        stochastic_options = stochastic_options,
        mefisto_options = mefisto_options
    )
    #### Run model ####
    model <- MOFA2::run_mofa(
        object = MOFAobject,
        outfile = outfile,
        use_basilisk = use_basilisk,
        ...
    )
    #### Add original names ####
    if(is.list(obj)){
        model@cache$obs_names <- lapply(obj,scKirby::get_obs_names) 
    } else {
        model@cache$obs_names <- scKirby::get_obs_names(obj)
    } 
    #### Add obs ####
    if (!is.null(obs)) {
        if (obs_idcol %in% colnames(obs)) {
            messager("Adding obs to model.", v = verbose)
            obs$sample <- obs[[obs_idcol]]
            MOFA2::samples_metadata(model) <- obs
        }
    }
    #### Reduce further ####
    model <- run_mofa2_reductions(
        model = model,
        reductions = reductions,
        verbose = verbose
    )
    #### Add MOFA/MOFA+UMAP reductions back into obj ####
    obj <- add_mofa2_dimred(
        obj = obj,
        model = model,
        assay = assay,
        keys = c("mofa",paste("mofa",reductions,sep="_")),
    )
    return(list(
        model = model,
        obj = obj
    ))
}
