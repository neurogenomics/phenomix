fix_colnames <- function(obj,
                         width = NULL,
                         make_unique = TRUE) {
    if (methods::is(obj, "list")) {
        for (i in length(obj)) {
            cnames <- scKirby::get_obs_names(obj = obj[[i]])
            cnames <- fix_names(
                names = cnames,
                width = width,
                make_unique = make_unique
            )
            obj[[i]] <- assign_colnames(
                obj = obj[[i]],
                cnames = cnames
            )
        }
    } else {
        cnames <- scKirby::get_obs_names(obj = obj)
        cnames <- fix_names(
            names = cnames,
            width = width,
            make_unique = make_unique
        )
        obj <- assign_colnames(
            obj = obj,
            cnames = cnames
        )
    }
    return(obj)
}
