test_that("map_snps2genes works", {
  
    dat <- MungeSumstats::formatted_example()
    #### txdb: gene-level ####
    dat2 <- map_snps2genes(dat,
                           method = "txdb", 
                           agg_var = "SYMBOL")
    testthat::expect_equal(nrow(dat2), 103)
    testthat::expect_true("ADJ_ZSTAT" %in% names(dat2))
    #### txdb: SNP-level ####
    dat3 <- map_snps2genes(dat,
                           method = "txdb", 
                           agg_var = c("SYMBOL","SNP"))
    testthat::expect_equal(nrow(dat3),104)
    testthat::expect_true("ADJ_ZSTAT" %in% names(dat3)) 
})
