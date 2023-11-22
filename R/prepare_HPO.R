#' Prepare Human Phenotype Ontology
#'
#' Conver the Human Phenotype Ontology (HPO) to a Seurat object with
#' extensive metadata.
#' @export
#' @import data.table
#' @import HPOExplorer
#' @examples
#' obj <- prepare_hpo(id_types="hpo_id")
prepare_hpo <- function(dt_genes = HPOExplorer::load_phenotype_to_genes(1),
                        dt_annot = HPOExplorer::load_phenotype_to_genes(3),
                        ctd = NULL,
                        # ctd = MAGMA.Celltyping::get_ctd("ctd_DescartesHuman"),
                        ct_results =
                            MultiEWCE::load_example_results("Descartes_All_Results_extras.rds"),
                        id_types = c("disease_id","hpo_id"),
                        min_quantile = 35,
                        min_genes = NULL,
                        value.var = "evidence_score_mean",
                        vars.to.regress = "n_genes",
                        run_nlp = FALSE,
                        run_impute = FALSE,
                        workers = 1,
                        seed = 2023,
                        verbose = TRUE){
    
    # o=devoptera::args2vars(prepare_hpo)
    top_celltype <- gene_symbol <- disease_id <- disease_db <- id <- NULL
    
    set.seed(seed)
    #### Create annotations ####
    if(is.null(dt_annot)){
        dt_annot <- HPOExplorer::make_phenos_dataframe(
            phenotype_to_genes = dt_genes,
            add_disease_data = TRUE,
            add_hoverboxes = FALSE)
    }
    #### Construct and bind gene_symbol x disease/phenotype matrices ####
    dt_genes <- HPOExplorer::add_evidence(phenos = dt_genes,
                                          verbose = verbose)
    X_ref <- prepare_hpo_matrix(
        dt_genes = dt_genes,
        id_types = grep("CellType",id_types, value = TRUE, invert = TRUE),
        value.var = value.var,
        verbose = verbose)
    #### Construct CTD matrix ####
    if(!is.null(ctd) ){
        messager("Integrating CTD data.",v=verbose)
        X_ctd <- lapply(stats::setNames(ctd,
                                        paste0("level",seq_len(length(ctd)))),
                        function(x){
                            x$specificity_quantiles[x$specificity_quantiles<min_quantile] <- 0
                            x$specificity_quantiles
                        })
        ctd_meta <- lapply(ctd,
                           function(x){
                               if(!is.null(x$annot)){
                                   data.table::as.data.table(dt_annot)
                               } else{
                                   data.table::data.table(id=colnames(x$mean_exp))
                               }
                           }) |> data.table::rbindlist(fill = TRUE,
                                                       use.names = TRUE,
                                                       idcol = "annotLevel")
        ctd_meta[,id_type:="CellType"]
        X_ctd2 <- X_ctd[[1]]
        for(x in X_ctd[-1]){
            X_ctd2 <- Seurat::RowMergeSparseMatrices(X_ctd2,x)
        }
        remove(X_ctd,x)
        X_ref <- Seurat::RowMergeSparseMatrices(X_ref, X_ctd2)
        remove(X_ctd2)
    }
    #### Construct metadata ####
    {
        if("disease_id" %in% names(dt_annot)){
            dt_annot[,disease_db:=sapply(disease_id,function(x){
                strsplit(x,":")[[1]][1]})]
        }
        if("hpo_id" %in% names(dt_annot)){
            dt_annot$hpo_name <-
                HPOExplorer::harmonise_phenotypes(dt_annot$hpo_id,
                                                  keep_order = TRUE,
                                                  verbose = verbose) 
        }
        
    }
    dt_annot_melt <- data.table::melt.data.table(
        dt_annot,
        measure.vars = id_types,
        variable.name = "id_type",
        value.name = "id") |>
        data.table::setkeyv("id")
    ## Ensure one row per ID
    dt_annot_melt <- dt_annot_melt[, head(.SD, 1), by="id"][colnames(X_ref),]
    ## Add id_name column for labeling
    id_names <- gsub("_id$","_name",id_types)
    id_names <- id_names[id_names %in% names(dt_annot)]
    dt_annot_melt[,id_name:=data.table::fcoalesce(as.list(id_names,"id"))]
    #### Add per sample gene_symbol counts ####
    dt_genes_melt <- data.table::melt.data.table(
        dt_genes,
        measure.vars = id_types,
        variable.name = "id_type",
        value.name = "id"
    )[,list(n_genes=length(unique(gene_symbol))), by="id"]
    data.table::setkeyv(dt_annot_melt,"id")
    dt_annot_melt <- dt_annot_melt[dt_genes_melt]
    #### Add CTD metadata ####
    if(exists("ctd_meta")){
        dt_annot_melt <- data.table::rbindlist(list(dt_annot_melt,
                                                    ctd_meta), fill=TRUE)
    }
    #### Add celltype enrichment results ####
    ct_results$CellType_fixed <-  EWCE::fix_celltype_names(ct_results$CellType,
                                                           make_unique = FALSE)
    
    #### Get
    # ct_results <- dplyr::group_by(ct_results, hpo_id) |>
    #   dplyr::arrange(p,dplyr::desc(fold_change)) |>
    #     dplyr::slice_head(n = 1) |> data.table::data.table()
    
    ct_results_cast <- data.table::dcast.data.table(
        ct_results,
        formula = hpo_id ~ paste("q",CellType_fixed,sep="."),
        fun.aggregate = mean,
        value.var = "q",
        na.rm = TRUE)
    #### Add
    ct_cols <- grep("^q\\.",names(ct_results_cast), value = TRUE)
    ct_results_cast[,
                    top_celltype:=gsub(
                        "^q\\.","",
                        ct_cols[apply(ct_results_cast[,ct_cols, with=FALSE],1,which.min)])
    ]
    ### Add top tissue ####
    ct_results_cast$top_tissue <- MultiEWCE::map_tissues(
        ct_results_cast$top_celltype,
        collapse = ";")
    #### Merge with rest of annotations ####
    dt_annot_melt <- data.table::merge.data.table(dt_annot_melt,
                                                  ct_results_cast,
                                                  by.x = "id",
                                                  by.y = "hpo_id",
                                                  all.x = TRUE)
    #### Filter out any terms with <x genes ####
    if(!is.null(min_genes)){
        X_ref <- X_ref[,Matrix::colSums(X_ref)>=min_genes]
    }
    #### Construct Seurat obj ####
    obj <- list(data=list("freq"=X_ref),
                obs=data.frame(
                    dt_annot_melt,
                    row.names = dt_annot_melt$id))
    ref <- scKirby::process_seurat(obj = obj,
                                   nfeatures = nrow(obj),
                                   vars.to.regress = vars.to.regress,
                                   default_assay = "freq",
                                   workers = workers)
    #### Impute ####
    if(isTRUE(run_impute)){
        messager("Running imputation.",v=verbose)
        ref <- SeuratWrappers::RunALRA(ref)
        ref@assays$alra@counts <-ref@assays$alra@data
        ref <- scKirby::process_seurat(obj = list("alra",ref),
                                       vars.to.regress = vars.to.regress,
                                       nfeatures = nrow(ref),
                                       default_assay = "alra",
                                       workers = workers)
    }
    #### Run scNLP ####
    if(isTRUE(run_nlp)){
        scnlp_res <- scNLP::plot_tfidf(
            ref,
            label_var = "definition",
            terms_per_cluster = 1,
            size_var = "n_genes",
            point_size = .5,
            point_palette = pals::kovesi.cyclic_mrybm_35_75_c68_s25(
                length(unique(ref@meta.data[["seurat_clusters"]]))
            )
        )
        scnlp_res$obj <- ref
        return(scnlp_res)
    }
    return(ref)
}
