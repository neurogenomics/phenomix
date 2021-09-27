fix_names <- function(names,
                      width = NULL,
                      make_unique = TRUE) {
    new_names <- iconv(gsub(" |[+]|[%]", "-", names),
        "latin1", "ASCII",
        sub = ""
    )
    # MOFA2 doesn't let you have IDs > 50 characters
    if (!is.null(width)) {
        new_names <- stringr::str_trunc(new_names,
            width = width
        )
    }
    if (make_unique) new_names <- make.unique(new_names)
    return(new_names)
}
