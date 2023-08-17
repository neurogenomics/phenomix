test_that("adjust_zstat works", {
  
    dat <- MungeSumstats::formatted_example()
    dat2 <- map_snps2genes(dat, adjust_z=FALSE)
    dat3 <- adjust_zstat(dat2)
    testthat::expect_false("ADJ_ZSTAT" %in% names(dat2))
    testthat::expect_true("ADJ_ZSTAT" %in% names(dat3))
})
