compare_snps2genes_methods <- function(res1,
                                       key1 = "GENE",
                                       res2,
                                       key2 = "SYMBOL",
                                       show_plot = TRUE,
                                       verbose = TRUE) {
    #### MAGMA ####
    # magma <- data.table::fread(magma_files$`finn-a-G6_ALZHEIMER`)
    # magma <- translate_magma_geneids(magma_dt = magma,
    #                                  gene_col = "GENE")
    # magma <- adjust_zstat_magma(magma)
    # data.table::setkey(magma,"GENE")
    #
    # res1 <- magma
    # res2 <- gene_hits

    #### Set keys ####
    data.table::setnames(res1, key1, "GENE")
    data.table::setkeyv(res1, "GENE")
    data.table::setnames(res2, key2, "GENE")
    data.table::setkeyv(res2, "GENE")

    #### Subset to shared genes ####
    shared_genes <- unique(intersect(res1$GENE, res2$GENE))
    messager(formatC(length(shared_genes), big.mark = ","),
        "shared genes identified.",
        v = verbose
    )
    res1 <- res1[shared_genes, ]
    res2 <- res2[shared_genes, ]
    #### Merge data ####
    res <- cbind(
        res1 = res1[, c("GENE", "ZSTAT", "ADJ_ZSTAT")],
        res2 = res2[, c("GENE", "ZSTAT", "ADJ_ZSTAT")]
    )
    res[, mean.ADJ_ZSTAT := (res1.ADJ_ZSTAT + res2.ADJ_ZSTAT) / 2]

    #### plot ####
    gp <- ggplot(res, aes(x = res1.ADJ_ZSTAT, y = res2.ADJ_ZSTAT, color = mean.ADJ_ZSTAT)) +
        geom_point(alpha = .5) +
        # ggecho::stat_echo(alpha=.1,
        #           geom = "point", size_increment = 2, alpha_factor = 0.5,
        #           x_offset = 0, y_offset = 0, n = 10) +
        geom_smooth(method = "lm") +
        scale_color_viridis_c() +
        ggpubr::stat_cor() +
        theme_bw()
    if (show_plot) print(gp)
    return(gp)
}
