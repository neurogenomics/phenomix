# 
# ewce_pipeline <- function(marker_list,
#                           ctd = ewceData::ctd(),
#                           level = length(ctd),
#                           bg = ewceData::all_hgnc(),
#                           sctSpecies = "human",
#                           genelistSpecies = "human",
#                           gsub_celltype = "mouse.Zeisel2018.",
#                           idcol = "trait") {
#     boot_res <- lapply(names(marker_list), function(x) {
#         print(x)
#         markers <- marker_list[[x]]
#         hits <- markers[(markers %in% bg) &
#             (markers %in% row.names(ctd[[level]]$mean_exp))]
#         res <- EWCE::bootstrap_enrichment_test(
#             sct_data = ctd,
#             annotLevel = level,
#             hits = hits,
#             bg = bg,
#             sctSpecies = sctSpecies,
#             genelistSpecies = genelistSpecies
#         )
#         res_df <- res$results %>%
#             dplyr::mutate(
#                 FDR = p.adjust(p = p, method = "fdr"),
#                 CellType = gsub(gsub_celltype, "", CellType)
#             )
#         return(res_df)
#     }) %>%
#         `names<-`(names(marker_list)) %>%
#         data.table::rbindlist(idcol = idcol)
#     return(boot_res)
# }
