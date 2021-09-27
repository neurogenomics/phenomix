
gprofiler_plots <- function(gp_list) {
    lapply(names(gp_list), function(x) {
        print(x)
        if ("result" %in% names(gp_list[[x]])) {
            return(gprofiler2::gostplot(gp_list[[x]]))
        } else {
            return(NULL)
        }
    })
}
