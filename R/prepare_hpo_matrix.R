prepare_hpo_matrix <- function(
        dt_genes = HPOExplorer::load_phenotype_to_genes(1),
        id_types,
        value.var,
        fill = 0,
        verbose = TRUE,
        ...){
    
    X_ref <- lapply(id_types,
                    function(id_type){
                        HPOExplorer::hpo_to_matrix(
                            phenotype_to_genes = dt_genes,
                            formula = paste("gene_symbol ~",id_type),
                            value.var = value.var,
                            fill = fill,
                            verbose = verbose,
                            ...)
                    })
    if(length(X_ref)==1){
        X_ref <- X_ref[[1]]
    } else {
        X_ref <- X_ref |> do.call(what = SeuratObject::RowMergeSparseMatrices)
    }
    return(X_ref)
}