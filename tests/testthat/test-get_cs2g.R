test_that("get_cs2g works", {
  
    #### UKB ####
    dataset <- "finemapping_cS2G_UKBB"
    out <- get_cs2g(dataset=dataset)
    testthat::expect_true(methods::is(out$data,'data.table'))
    testthat::expect_true(methods::is(out$obs,'data.table'))
    out_mtx <- get_cs2g(dataset=dataset, 
                        as_matrix=TRUE)
    testthat::expect_true(methods::is(out_mtx$data,'sparseMatrix'))
    out_gr <- get_cs2g(dataset=dataset, 
                       as_granges=TRUE)
    testthat::expect_true(methods::is(out_gr$data,'GRanges'))
    
    #### GWAS Catalog ####
    dataset <- "gwas_catalog_cS2G"
    out <- get_cs2g(dataset=dataset)
    testthat::expect_true(methods::is(out$data,'data.table'))
    testthat::expect_true(methods::is(out$obs,'data.table'))
    out_mtx <- get_cs2g(dataset=dataset, 
                        as_matrix=TRUE)
    testthat::expect_true(methods::is(out_mtx$data,'sparseMatrix'))
    out_gr <- get_cs2g(dataset=dataset, 
                       as_granges=TRUE)
    testthat::expect_true(methods::is(out_gr$data,'GRanges'))
})
