#' Prepare Human Phenotype Ontology
#'
#' Convert the Human Phenotype Ontology (HPO) to a \link{Seurat} object 
#' with extensive metadata.
#' @param dt_genes \link{data.table} of phenotype to gene associations.
#' @param dt_annot \link{data.table} of phenotype to disease annotations.
#' @param ctd \link{list} of \link{data.table} of cell type specificity data.
#' @param ct_results \link{data.table} of cell type enrichment results.
#' @param id_types \link{character} vector of identifiers to use for
#' disease and HPO terms.
#' @param min_quantile \link{numeric} minimum quantile for filtering
#' cell type enrichment results.
#' @param min_genes The minimum number of genes of phenotypes/diseases
#' to include from the data matrix.
#' @param min_value The minimum sum of columns to include from the data matrix.
#' @param celltype_col Name of the cell type column. 
#' @param run_nlp Label clusters using natural language processing via 
#' \link[scNLP]{plot_tfidf}.
#' @param run_impute Run imputation on the \link{Seurat} object using
#'  \link[SeuratWrappers]{RunALRA}.
#' @inheritParams HPOExplorer::hpo_to_matrix
#' @inheritParams HPOExplorer::make_phenos_dataframe
#' @inheritParams scKirby::process_seurat
#' @inheritDotParams scKirby::process_seurat
#' 
#' @export
#' @import data.table
#' @import HPOExplorer
#' @examples
#' obj <- prepare_hpo(id_types="hpo_id", min_genes=10)
prepare_hpo <- function(dt_genes = HPOExplorer::load_phenotype_to_genes(1),
                        dt_annot = HPOExplorer::load_phenotype_to_genes(3),
                        hpo = HPOExplorer::get_hpo(),
                        ctd = NULL, # MSTExplorer::example_ctd(),
                        ct_results =
                            MSTExplorer::load_example_results()|>
                            MSTExplorer::map_celltype(),
                        id_types = c("hpo_id"),#c("disease_id","hpo_id"),
                        min_quantile = 35,
                        min_genes = NULL,
                        min_value = NULL,
                        celltype_col=c("cl_name","cl_id","CellType"),
                        celltype_value="q",
                        value.var = "evidence_score_sum",
                        default_assay = "score",
                        nfeatures = NULL,
                        vars.to.regress = NULL, # paste0("nFeature_",default_assay),
                        run_nlp = FALSE,
                        run_impute = FALSE,
                        workers = 1,
                        seed = 2023,
                        verbose = TRUE,
                        save_path=NULL,
                        force_new=FALSE,
                        ...){ 
    top_celltype <- gene_symbol <- disease_id <- disease_db <- id <- 
        ancestor_name <- ancestor_name_abnormality <- NULL
    
    celltype_col <- celltype_col[1]
    #### Check for existing file ####
    if(!is.null(save_path) &&
       file.exists(save_path) && 
       isFALSE(force_new)){
        messager("Loading precomputed data:",save_path)
        return(readRDS(save_path))
    }
    set.seed(seed)
    #### Create annotations ####
    if(is.null(dt_annot)){
        dt_annot <- HPOExplorer::make_phenos_dataframe(
            phenotype_to_genes = dt_genes,
            hpo = hpo,
            add_disease_data = TRUE,
            add_hoverboxes = FALSE)
    }
    dt_annot <- HPOExplorer::add_hpo_definition(phenos = dt_annot,
                                                hpo = hpo)
    dt_annot <- HPOExplorer::add_ancestor(phenos = dt_annot,
                                          hpo = hpo)
    
    dt_annot[,ancestor_name_abnormality:=ifelse(
        grepl("abnormality",ancestor_name,
              ignore.case = TRUE),
        ancestor_name,
        NA
    )] 
    #### Construct and bind gene_symbol x disease/phenotype matrices ####
    dt_genes <- HPOExplorer::add_evidence(phenos = dt_genes)
    X_ref <- prepare_hpo_matrix(
        dt_genes = dt_genes,
        id_types = grep(celltype_col, id_types, value = TRUE, invert = TRUE),
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
                HPOExplorer::map_phenotypes(terms = dt_annot$hpo_id,
                                            hpo = hpo,
                                            keep_order = TRUE) 
        }
        
    }
    {
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
        dt_annot_melt[,name:=data.table::fcoalesce(.SD),.SDcols=c(id_names,"id")]
    }
    {
        #### Add per sample gene_symbol counts ####
        dt_genes_melt <- data.table::melt.data.table(
            dt_genes,
            measure.vars = id_types,
            variable.name = "id_type",
            value.name = "id"
        )[,list(n_genes=length(unique(gene_symbol))), by="id"]
        data.table::setkeyv(dt_annot_melt,"id")
        dt_annot_melt <- dt_annot_melt[dt_genes_melt]
    }
    #### Add CTD metadata ####
    if(exists("ctd_meta")){
        dt_annot_melt <- data.table::rbindlist(list(dt_annot_melt,
                                                    ctd_meta), fill=TRUE)
    }
    #### Add celltype enrichment results ####
    # ct_results$CellType <-  EWCE::fix_celltype_names(ct_results$CellType,
    #                                                  make_unique = FALSE)
    
    #### Get
    # ct_results <- dplyr::group_by(ct_results, hpo_id) |>
    #   dplyr::arrange(p,dplyr::desc(fold_change)) |>
    #     dplyr::slice_head(n = 1) |> data.table::data.table() 
    ct_results_cast <- data.table::dcast.data.table(
        ct_results,
        formula = hpo_id ~ paste(celltype_value,get(celltype_col),sep="."),
        fun.aggregate = mean,
        value.var = celltype_value,
        na.rm = TRUE)
    #### Add
    ct_cols <- grep(paste0("^",celltype_value,"\\."),names(ct_results_cast), value = TRUE)
    ct_results_cast[,
                    top_celltype:=gsub(
                        paste0("^",celltype_value,"\\."),"",
                        ct_cols[apply(ct_results_cast[,ct_cols, with=FALSE],1,
                                      which.min)])
    ]
    ### Add top tissue ####
    # tissues <- MSTExplorer::map_tissue(return_agg=TRUE)
    #### Merge with rest of annotations ####
    dt_annot_melt <- data.table::merge.data.table(dt_annot_melt,
                                                  ct_results_cast,
                                                  by.x = "id",
                                                  by.y = "hpo_id",
                                                  all.x = TRUE) 
    #### Filter out any terms with <x genes ####
    if(!is.null(min_genes)){
        X_ref <- X_ref[,Matrix::colSums(X_ref!=0)>min_genes]
    }
    if(!is.null(min_value)){
        X_ref <- X_ref[,Matrix::colSums(X_ref)>=min_value]
    }
    #### Construct Seurat obj #### 
    dt_annot_melt <- dt_annot_melt[!is.na(id)]
    shared_ids <- intersect(colnames(X_ref),
                            dt_annot_melt$id)
    obj <- list(data=list(X_ref[,shared_ids])|>`names<-`(default_assay),
                obs=data.frame(
                    dt_annot_melt,
                    row.names = dt_annot_melt$id)[shared_ids,]
                ) 
    ref <- scKirby::process_seurat(obj = obj,
                                   nfeatures = nfeatures,
                                   vars.to.regress = vars.to.regress,
                                   default_assay = default_assay,
                                   workers = workers,
                                   ...)
    ref$orig.ident <- "HPO"
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
    #### Save ####
    KGExplorer::cache_save(ref,save_path)
    #### Run scNLP ####
    if(isTRUE(run_nlp)){ 
        ref$name_definition <- paste(ref$name,
                                     ref$definition)
        scnlp_res <- scNLP::plot_tfidf(
            ref,
            label_var = "name_definition",
            terms_per_cluster = 1,
            size_var = paste0("nFeature_",default_assay),
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
