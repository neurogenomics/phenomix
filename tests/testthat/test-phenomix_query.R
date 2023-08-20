test_that("phenomix_query works", {
  
    query_granges <- GenomicRanges::GRanges(c("2:15000-17000","3:60000-70000"))
    res <- phenomix_query(ids="bbj-a-1",
                          query_granges=query_granges)
    testthat::expect_equal(nrow(res$`bbj-a-1`),22)
    
    # ids <- tail(phenomix::OpenGWAS[!grepl("^finn",id),]$id,3)
    # res_ukb <- phenomix_query(ids=ids,
    #                           query_granges="cs2g_ukb", 
    #                           query_method="conda", 
    #                           overlapping_only = TRUE) 
    
    meta <- opengwas_meta()
    meta <- meta[grepl("eqtl-a-",id),][seq(5)]
    res2 <- phenomix_query(meta = meta, 
                           data_type = "magma",
                           as_matrix = TRUE)
    testthat::expect_equal(dim(res2),c(5043,5))
})
