melt_matrix <- function(X,
                        row.name = "feature",
                        variable.name = "variable",
                        value.name = "value"){
    X |>
        as.matrix() |>
        data.table::as.data.table(keep.rownames = row.name) |> 
        data.table::melt.data.table(id.vars = row.name, 
                                    variable.name = variable.name,
                                    value.name = value.name) |>
        data.table::setkeyv(row.name)
}