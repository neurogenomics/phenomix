#' Get data
#' 
#' Get data stored on GitHub Releases via \pkg{piggyback}.
#' @param piggyback_cache_duration piggyback cache duration.
#' @returns R object.
#' 
#' @inheritParams piggyback::pb_download
#' @keywords internal
#' @importFrom piggyback pb_download
#' @importFrom tools R_user_dir
get_data <- function(file,
                     repo = "neurogenomics/phenomix",
                     tag = "latest",
                     save_dir = tools::R_user_dir(package = "phenomix",
                                                  which = "cache"),
                     .token = gh::gh_token(),
                     piggyback_cache_duration = 10,
                     overwrite = FALSE) {
    
    requireNamespace("gh")
    
    tmp <- file.path(save_dir, file) 
    if (!file.exists(tmp)) {
        dir.create(save_dir,showWarnings = FALSE,recursive = TRUE)
        Sys.setenv("piggyback_cache_duration" = piggyback_cache_duration)
        piggyback::pb_download(
            file = file,
            dest = save_dir,
            repo = repo,
            tag = tag,
            overwrite = overwrite
        )
    }
    return(tmp)
}
