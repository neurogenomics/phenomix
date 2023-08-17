test_that("find_neighbors works", {
  
    obj <- get_HPO()[seq(100),]
    top_neighbors <- find_neighbors(
        obj = obj,
        var1_search = "parkinson",
        label_col = "HPO_label"
    )
    testthat::expect_equal(nrow(top_neighbors),70)
    testthat::expect_equal(top_neighbors$trait1[1],"Parkinsonism")
    testthat::expect_equal(top_neighbors$trait2[1],"Hallucinations")
    #### Error when col not present ####
    testthat::expect_error(
        find_neighbors(
            obj = obj, 
            label_col = "typoooo"
        )
    )
})
