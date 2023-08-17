# 
# find_corr_celltypes <- function(integrated,
#                                 target_pheno,
#                                 show_plot = T) {
#     celltypes <- colnames(integrated)[integrated$type == "celltype"]
#     phenotypes <- grep(target_pheno, colnames(integrated), ignore.case = T, value = T)
# 
#     cor_integrated <- run_corr(GetAssayData(integrated)[, c(celltypes, phenotypes)],
#         pval_thresh = .05 / (length(celltypes) * length(phenotypes))
#     )
#     diag(cor_integrated) <- NA
# 
#     cor_integrated_melt <- reshape2::melt(cor_integrated) |> `colnames<-`(c("trait1", "trait2", "r"))
#     top_celltypes <- subset(cor_integrated_melt, (trait1 %in% phenotypes) & (trait2 %in% celltypes)) |>
#         dplyr::group_by(trait1) |>
#         dplyr::slice_max(order_by = r, n = 5) |>
#         data.frame()
#     ### Get mean r for ordering
#     mean_r <- top_celltypes |>
#         dplyr::group_by(trait2) |>
#         dplyr::summarise(r = median(r, na.rm = T)) |>
#         dplyr::arrange(desc(r))
# 
#     top_celltypes$trait2 <- factor(top_celltypes$trait2, levels = mean_r$trait2, ordered = T)
# 
#     gp <- ggplot(top_celltypes, aes(x = trait2, y = r)) +
#         geom_boxplot() +
#         # facet_wrap(facets = .~trait1, ncol = 1) +
#         scale_fill_gradient(low = "blue", high = "red") +
#         # scale_fill_brewer(palette="BuPu") +
#         labs(title = "Specificity correlation", subtitle = target_pheno, x = "Cell-type") +
#         theme_bw() +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1))
#     if (show_plot) print(gp)
#     return(list(
#         data = top_celltypes,
#         plot = gp
#     ))
# }
