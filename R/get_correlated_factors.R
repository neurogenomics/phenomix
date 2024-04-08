get_correlated_factors <- function(obj,
                                   keys,
                                   metadata_var,
                                   p_threshold=.05,
                                   r_threshold=.1){
    requireNamespace("Hmisc")
    obsm <- scKirby::get_obsm(obj,
                              keys = keys,
                              n=1)
    obsm <- cbind(obsm, 
                  nFeature_score=obj$nFeature_score[rownames(obsm)])
    obsm_rcor <- Hmisc::rcorr(obsm)
    obsm_cor_sig <- sort(
        abs(obsm_rcor$r[metadata_var,][obsm_rcor$P[metadata_var,]<p_threshold]), 
        decreasing = TRUE
    ) 
    omit_dims <- as.numeric(gsub(".*_","",
                                 names(obsm_cor_sig[obsm_cor_sig>r_threshold]))
                            )
    return(omit_dims)
}