#' Reconstruct original matrix from DeGAs components
#'
#' @keywords internal
#' @importFrom stats quantile
#' @importFrom scales rescale
#' @importFrom dplyr %>%
#' @importFrom reticulate py
reconstruct_matrix <- function(u_name = "contribution_phe",
                               v_name = "contribution_var",
                               d_name = NULL,
                               DEGAS = NULL,
                               annot = NULL,
                               gene_level = FALSE,
                               remove_nongenes = TRUE,
                               filter_quantiles = FALSE,
                               norm = TRUE,
                               norm_scale = c(1, 100),
                               replace_underscores = TRUE) {
    if (is.null(DEGAS)){
        py <- reticulate::py
        DEGAS <- py$DEGAS
    }
    label_phe_code <- DEGAS["label_phe_code"]
    if (gene_level) {
        ## NOTE:  Genes symbols are provided directly in label_gene
        message("+ Annotating at the gene-level")
        label_var_gene <- DEGAS["label_gene"]
    } else {
        ## NOTE: Variants must be translated from positions to RSIDs
        ## Confirmed that none of label_var start with "rs"
        label_var_gene <- annotate_vector(
            vec = DEGAS["label_var"],
            annot = annot
        ) %>% unname()
    } 
    if (replace_underscores) {
        message("+ Replacing _ with -")
        ## Seurat requires this, and tries 
        ## to enforce it but doesn't do it very well
        ## https://github.com/satijalab/seurat/issues/1296
        label_var_gene <- gsub("_", "-", label_var_gene)
    }
    message("+ Preparing input matrices")
    u <- DEGAS[u_name]
    label_dim <- paste0("dim", seq(1,ncol(u)))
    u <- u %>%
        `rownames<-`(label_phe_code) %>%
        `colnames<-`(label_dim)
    v <- DEGAS[v_name] %>%
        `rownames<-`(label_var_gene) %>%
        `colnames<-`(label_dim)
    if (gene_level & remove_nongenes) {
        ## NOTE: Most of the "genes" do not have any gene symbols. Remove them.
        message("+ Removing non-genes")
        yes_genes <- grep(paste(paste0("^", seq(1,23), "-"), collapse = "|"), 
                          rownames(v), value = TRUE, invert = TRUE)
        v <- v[yes_genes, ]
    }
    if (is.null(d_name)) {
        message("+ Reconstructing matrix from `u %*% t(v)`")
        M2 <- u %*% t(v)
    } else {
        print("+ Reconstructing matrix from `u %*% diag(d) %*% t(v)`")
        d <- DEGAS[d_name]
        M2 <- u %*% diag(d) %*% t(v)
    }
    if (filter_quantiles) {
        print("+ Removing values in the lowest percentile")
        quants <- stats::quantile(abs(M2), 
                                  probs = seq(0, 1, length.out = 100))
        M2_remove <- abs(M2) < unname(quants[2])
        M2[abs(M2) < quants[2]] <- NA
    }

    if (norm) {
        print("+ Normalizing data")
        M2 <- scales::rescale(M2, to = norm_scale)
    }
    return(list(
        M2 = M2,
        embeddings = u,
        loadings = v,
        stdev = if (is.null(d_name)) NULL else d
    ))
}
