#' Generate metadata
#'
#' Generate metadata for \href{https://gwas.mrcieu.ac.uk/}{OpenGWAS} datasets.
#'
#' @param munged_files List of paths to summary statostics files
#' that have been munged with \link[MungeSumstats]{import_sumstats}.
#' @param infer_builds Infer genome builds of each file.
#' @param drop_na_N Drop rows without "N" (total sample size).
#' @param trait_groups [Optional] A vector of strings.
#' Creates a new column "trait_groups"
#' based on a series of substring searches.
#' @param force_new Create a new metadata file even when one
#' already exists in \code{save_dir}.
#' @param save_dir Directory to save metadata to.
#' @param verbose Print messages.
#'
#' @return data.frame of OpenGWAS metadata.
#'
#' @export
#' @importFrom MungeSumstats find_sumstats get_genome_builds
#' @importFrom utils read.csv write.csv
#' @importFrom  dplyr mutate %>%
generate_metadata <- function(munged_files,
                              infer_builds = TRUE,
                              drop_na_N = TRUE,
                              trait_groups = NULL,
                              force_new = FALSE,
                              save_dir = tempdir(),
                              verbose = TRUE) {
    sample_size <- ncase <- ncontrol <- N <- NULL;
    save_path <- file.path(save_dir, "OpenGWAS_metadata.csv")
    #### Import existing file ####
    if (file.exists(save_path) & (!force_new)) {
        messager("Importing existing metadata: ", save_path)
        metagwas_all <- utils::read.csv(save_path, row.names = 1)
        return(metagwas_all)
    } 
    #### Collected munged G#WAS ####
    munged_ids <- gsub(".tsv.gz", "", basename(munged_files))
    names(munged_files) <- munged_ids
    #### Search OpenGWAS for thei IDs ####
    metagwas_all <- MungeSumstats::find_sumstats(ids = munged_ids)
    #### Add munged GWAS path into metadata ####
    metagwas_all$path <- munged_files[metagwas_all$id]
    #### Infer N ####
    metagwas_all <- metagwas_all %>%
        dplyr::mutate(N = ifelse(is.na(sample_size),
            sum(ncase, ncontrol, na.rm = TRUE), sample_size
        ))
    if (drop_na_N) {
        metagwas_all <- subset(metagwas_all, !is.na(N))
    }
    #### Creare trait groups ####
    if (!is.null(trait_groups)) {
        metagwas_all$trait_group <- stringr::str_extract(
            tolower(metagwas_all$trait),
            paste0("(", tolower(trait_groups), ")", collapse = "|")
        )
        metagwas_all[is.na(metagwas_all$trait_group), "trait_group"] <- "other"
    }
    #### Set rownames ####
    metagwas_all <- data.frame(metagwas_all) %>%
        `rownames<-`(gsub("-|[:]|[.]", "_", metagwas_all$id))
    #### Infer genome ####
    if (infer_builds & (!"build_inferred" %in% colnames(metagwas_all))) {
        builds <- MungeSumstats::get_genome_builds(
            sumstats_list = munged_files,
            header_only = TRUE
        )
        metagwas_all$build_inferred <- unlist(builds)
    }
    #### Save ####
    if (!is.null(save_dir)) {
        utils::write.csv(metagwas_all, save_path)
    }
    return(metagwas_all)
}
