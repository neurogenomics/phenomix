report_time <- function(start,
                        verbose = TRUE) {
    messager(utils::capture.output(difftime(Sys.time(), start)), v = verbose)
}
