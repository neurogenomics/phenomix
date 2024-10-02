#' Run sparse SVD
#'
#' Run sparse Singular Value Decomposition (sparseSVD).
#'
#' Uses \link[sparsesvd]{sparsesvd}.
#'
#' @param mat Matrix to run sparseSVD on.
#' @param transpose Whether to transpose the matrix first.
#' @param add_names Add colnames and rownames to embeddings and loadings.
#' @inheritParams sparsesvd::sparsesvd
#'
#' @importFrom Matrix t
#' @importFrom sparsesvd sparsesvd
#'
#' @export
#' @examples
#' obj <- get_HPO()[,1:50]
#' obj2 <-run_sparsesvd(obj, rank=50L)
run_sparsesvd <- function(obj,
                          transpose = TRUE,
                          add_names = TRUE,
                          rank = 50L,
                          tol = 1e-15,
                          kappa = 1e-6,
                          assay=NULL,
                          layer=NULL) {
    ## Extract matrix 
    X <- scKirby::get_x(obj,
                        n=1,
                        simplify = TRUE,
                        as_sparse = TRUE,
                        layer = layer,
                        assay = assay,
                        transpose = transpose)
    ## Run SSVD
    messager("Running SparseSVD.")
    ssvd <- sparsesvd::sparsesvd(
        M = X,
        rank = rank,
        tol = tol,
        kappa = kappa
    )
    ## Add names to submatrices
    if (add_names) {
        ssvd_names <- paste0("SSVD_", seq(1, ncol(ssvd$u)))
        ssvd$u <- ssvd$u |>
            `colnames<-`(ssvd_names) |>
            `row.names<-`(rownames(X))
        ssvd$v <- ssvd$v |>
            `colnames<-`(ssvd_names) |>
            `row.names<-`(colnames(X))
        ssvd$d <- setNames(ssvd$d, ssvd_names)
    }
    ## Add DimReducObject
    if(scKirby::is_class(obj,"seurat")){
        obj[["ssvd"]] = SeuratObject::CreateDimReducObject(
            embeddings = ssvd$u,
            loadings = ssvd$v,
            key = "ssvd_",
            assay=Seurat::DefaultAssay(obj)
        )
        return(obj)
    } else {
        return(ssvd)
    }
    # obs <- merge(ssvd$u, obj@meta.data[,c("nFeature_score"),drop=FALSE], by="row.names")[,-1]
    # Xcor=WGCNA::cor(obs)[,'nFeature_score']
    # sort(Xcor)
    return(ssvd)
}
