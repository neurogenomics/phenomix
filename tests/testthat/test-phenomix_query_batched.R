test_that("phenomix_query_batched works", {
  
    query_granges <- GenomicRanges::GRanges(c("2:15000-17000","3:60000-70000"))
    ids <- tail(phenomix::OpenGWAS[!grepl("^finn",id),]$id,4)
    X <- phenomix_query_batched(ids=ids,
                                query_granges=query_granges,
                                batch_size=2)
    testthat::expect_equal(dim(X),c(16,3))
    testthat::expect_true(methods::is(X,"sparseMatrix"))
})
