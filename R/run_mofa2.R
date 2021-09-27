#' Run MOFA2
#'
#' Run Multi-Omics Factor Analysis v2 (MOFA2).
#'
#' Uses \link[MOFA2]{run_mofa} to reduce a single dataset,
#' or multiple views of the same dataset, into a
#' set number of factors.
#'
#' @param obj Named list of input data matrices, or a \pkg{Seurat} object.
#' @param transpose First transpose each matrix in \code{mat_list}.
#' @param maxiter Maximum number of training iterations (\emph{DEFAULT} = 1000).
#' @param reductions In addition to returning  MOFA factors,
#' further reduce the data with "umap" and/or "tsne".
#' @param seed Seed passed to \[base]{set.seed}
#' for reproducibility between runs.
#' @param verbose Print messages.
#' @param ... Additional parameters passed to \link[MOFA2]{run_mofa}.
#' @inheritParams MOFA2::run_mofa
#' @inheritParams MOFA2::create_mofa
#'
#' @source \href{https://biofam.github.io/MOFA2/}{MOFA2 site}
#'
#' @export
#' @importFrom MOFA2 create_mofa prepare_mofa samples_metadata
#' @importFrom MOFA2 get_default_data_options
#' @importFrom MOFA2 get_default_model_options
#' @importFrom MOFA2 get_default_training_options
#' @importFrom Matrix t
#' @importFrom methods is
run_mofa2 <- function(obj,
                      groups = NULL,
                      metadata = NULL,
                      assay = NULL,
                      slot = NULL,
                      features = NULL,
                      metadata_idcol = "sample",
                      transpose = FALSE,
                      maxiter = 1000,
                      data_options = NULL,
                      model_options = NULL,
                      training_options = NULL,
                      stochastic_options = NULL,
                      mefisto_options = NULL,
                      outfile = file.path(tempdir(), "model.hdf5"),
                      reductions = c("umap"),
                      seed = 2020,
                      verbose = TRUE,
                      ...) {
    if (!is.null(seed)) set.seed(seed)
    #### Transpose matrices ####
    if (transpose & methods::is(obj, "list")) {
        for (i in seq(1, length(mat_list))) {
            mat_list[[i]] <- Matrix::t(mat_list[[i]])
        }
    }
    #### Fix colnames ####
    # Can't be non-ASCI
    # Can't be >50 characters
    obj <- fix_colnames(
        obj = obj,
        width = 40,
        make_unique = TRUE
    )
    #### Prepare groups ####
    metadata <- if (is.null(metadata)) extract_metadata(obj)
    groups <- if ((!is.null(groups)) & length(groups) == 1) metadata[[groups]]
    # #### Use all features ####
    # if(features[1]=="all") {
    #     obj@assays[[assay]]@var.features <- rownames(
    #         obj@assays[[assay]]@meta.features)
    # }
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
        ...
    )
    #### Add metadata ####
    if (!is.null(metadata)) {
        if (metadata_idcol %in% colnames(metadata)) {
            messager("Adding metadata to model.", v = verbose)
            metadata$sample <- metadata[[metadata_idcol]]
            MOFA2::samples_metadata(model) <- metadata
        }
    }
    #### Reduce further ####
    model <- run_mofa2_reductions(
        model = model,
        reductions = reductions,
        verbose = verbose
    )
    #### Add MOFA/MOFA+UMAP reductions back into obj ####
    for (x in c("mofa2", "umap")) {
        obj <- add_mofa2_dimred(
            obj = obj,
            model = model,
            assay = assay,
            reduction = x
        )
    }
    return(list(
        model = model,
        obj = obj
    ))
}
