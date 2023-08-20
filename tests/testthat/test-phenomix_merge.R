test_that("phenomix_merge works", {

    query_granges <- GenomicRanges::GRanges(c("2:15000-17000","3:60000-70000"))
    ids <- tail(phenomix::OpenGWAS[!grepl("^finn",id),]$id,3)
    res <- phenomix_query(ids=ids, 
                          query_granges=query_granges)
    X <- phenomix_merge(res)
    testthat::expect_equal(dim(X),c(95,2))
    Xdt <- phenomix_merge(res, as_matrix = FALSE)
    testthat::expect_equal(dim(Xdt),c(95,3))
})
