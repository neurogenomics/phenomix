prepare_HPO <- function(save_path = tempfile(
                            pattern = "HPO_seurat",
                            fileext = ".rds"
                        ),
                        as_seurat = TRUE) {
    HPO <- data.table::fread(
        file.path(
            "https://ci.monarchinitiative.org/view/hpo/job",
            "hpo.annotations/lastSuccessfulBuild/artifact",
            "rare-diseases/util/annotation/phenotype_to_genes.txt"
        ),
        nThread = 10
    )
    #### Fix names ####
    colnames(HPO) <- gsub(
        "-|[ ]", "_",
        stringr::str_split(
            gsub("#Format: ", "", colnames(HPO)[1]),
            "<tab>"
        )[[1]]
    )
    HPO <- HPO %>% dplyr::mutate(
        HPOid = HPO_id,
        HPO_id = gsub("[:]", ".", HPO_id)
    )
    #### Prepare metadata ####
    meta <- HPO %>%
        dplyr::group_by(HPO_id, HPOid, HPO_label) %>%
        dplyr::summarise(n_genes = dplyr::n_distinct(entrez_gene_symbol)) %>%
        data.frame()
    #### Get ontology depths ####
    depths <- get_ontology_depths(ontology = hpo)
    meta$depth <- unname(depths[meta$HPOid])
    #### Assign groups ####
    meta_groups <- get_ontology_groups(
        ontology = hpo,
        depth_levels = seq(2, 4)
    )
    meta <- merge(meta, meta_groups, by.x = "HPOid", by.y = "id")
    # Ancestors
    n_ancestors <- get_n_ancestors(ontology = hpo)
    meta$n_ancestors <- unname(n_ancestors[meta$HPOid])
    # Descendents
    n_descendants <- get_n_descendants(ontology = hpo)
    meta$n_descendants <- unname(n_descendants[meta$HPOid])

    ## Created by Bobby: https://github.com/ovrhuman/RD_EWCE_Website_and_apps/blob/master/Cell_select_interactive/data/disease_descriptions.rds
    info <- readRDS("/Desktop/disease_descriptions.rds")
    meta <- merge(meta, info, by.x = "HPOid", by.y = "HPO_id")
    rownames(meta) <- meta$HPO_id
    #### Prepare matrix ####
    HPO$val <- 1
    mat <- data.table::dcast.data.table(HPO,
        formula = entrez_gene_symbol ~ HPO_id,
        value.var = "val",
        fill = 0,
        fun.aggregate = max
    ) %>%
        data.frame() %>%
        tibble::column_to_rownames(var = "entrez_gene_symbol") %>%
        data.frame()
    mat <- as(as.matrix(mat), "sparseMatrix")

    if (as_seurat) {
        obj <- scNLP::seurat_pipeline(
            counts = mat,
            meta.data = meta,
            add_specificity = TRUE
        )
    } else {
        obj <- list(
            mat = mat,
            metadata = ancest_df
        )
    }
    if (!is.null(save_path)) {
        message("Saving HPO Seurat object to ==>", save_path)
        dir.create(dirname(save_path),
            showWarnings = FALSE,
            recursive = TRUE
        )
        saveRDS(obj, file = save_path)
    }
    return(obj)
}
