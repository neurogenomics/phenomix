#' Reconstruct original matrix from DeGAs components
#'
#' @keywords internal
#' @importFrom stats quantile
#' @importFrom scales rescale
#' @importFrom dplyr %>%
reconstruct_matrix <- function(u_name = "contribution_phe",
                               v_name = "contribution_var",
                               d_name = NULL,
                               u_labels = "label_phe_code",
                               v_labels = if(gene_level){
                                   "label_gene"
                                   } else {"label_var"},
                               npz_obj = NULL,
                               annot = NULL,
                               gene_level = FALSE,
                               remove_nongenes = TRUE,
                               filter_quantiles = FALSE,
                               norm = TRUE,
                               norm_scale = c(1, 100),
                               replace_underscores = c(u=TRUE, v=TRUE),
                               translate_snp_positions = FALSE,
                               invert_UV_labels = FALSE,
                               verbose = TRUE) {
    requireNamespace("reticulate")
    if (is.null(npz_obj)){
        py <- reticulate::py
        npz_obj <- py$DEGAS
    }
    label_phe_code <- npz_obj[u_labels]
    if (gene_level) {
        ## NOTE:  Genes symbols are provided directly in label_gene
        messager("Annotating at the gene-level.",v=verbose)
        label_var_gene <- npz_obj[v_labels]
    } else {
        ## NOTE: Variants must be translated from positions to RSIDs
        ## Confirmed that none of label_var start with "rs"
        if(translate_snp_positions){
            ## Only necessary for DEGAS (not dPRS)
            label_var_gene <- annotate_vector(
                vec = npz_obj[v_labels],
                annot = annot
            ) %>% unname()
        } else {
            label_var_gene <- npz_obj[v_labels]
        }
    } 
    if(invert_UV_labels){
        messager("Inverting U and V labels.",v=verbose)
        label_phe_code_tmp <- label_phe_code
        label_var_gene_tmp <- label_var_gene
        label_phe_code <- label_var_gene_tmp
        label_var_gene <- label_phe_code_tmp
    }
    if (any(replace_underscores)) {
        message("Replacing _ with -.")
        ## Seurat requires this, and tries 
        ## to enforce it but doesn't do it very well
        ## https://github.com/satijalab/seurat/issues/1296
        if("u" %in% names(replace_underscores) &&
           isTRUE(replace_underscores["u"])){
            label_var_gene <- gsub("_", "-", label_var_gene)
        }
        if("v" %in% names(replace_underscores) &&
           isTRUE(replace_underscores["v"])){
            label_phe_code <- gsub("_", "-", label_phe_code)
        } 
    }
    messager("Preparing input matrices.", v=verbose) 
    u <- npz_obj[u_name]
    label_dim <- paste0("dim", seq(1,ncol(u)))
    u <- u %>%
        `rownames<-`(label_phe_code) %>%
        `colnames<-`(label_dim)
    v <- npz_obj[v_name] %>%
        `rownames<-`(label_var_gene) %>%
        `colnames<-`(label_dim)
    if (gene_level & remove_nongenes) {
        ## NOTE: Most of the "genes" do not have any gene symbols. Remove them.
        messager("Removing non-genes.",v=verbose)
        yes_genes <- grep(paste(paste0("^", seq(1,23), "-"), collapse = "|"), 
                          rownames(v), value = TRUE, invert = TRUE)
        v <- v[yes_genes, ]
    }
    if (is.null(d_name)) {
        messager("Reconstructing matrix from `u %*% t(v)`.",v=verbose)
        M2 <- u %*% t(v)
    } else {
        messager("Reconstructing matrix from `u %*% diag(d) %*% t(v)`.",
                 v=verbose)
        d <- npz_obj[d_name]
        M2 <- u %*% diag(d) %*% t(v)
    }
    if (filter_quantiles) {
        messager("Removing values in the lowest percentile.",v=verbose)
        quants <- stats::quantile(abs(M2), 
                                  probs = seq(0, 1, length.out = 100))
        M2_remove <- abs(M2) < unname(quants[2])
        M2[abs(M2) < quants[2]] <- NA
    }

    if (norm) {
        messager("Normalizing data.",v=verbose)
        M2_norm <- scales::rescale(M2, to = norm_scale)
    } else {
        M2_norm <- NULL
    }
    return(list(
        M2 = M2,
        M2_norm = M2_norm,
        embeddings = u,
        loadings = v,
        stdev = if (is.null(d_name)) NULL else d
    ))
}
