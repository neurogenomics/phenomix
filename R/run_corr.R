run_corr <- function(mat, 
                     type="spearman",
                     target_str=NULL,
                     target_str_in_rows=T,
                     pval_thresh=1){ 
    cor_res <- Hmisc::rcorr(as.matrix(mat))  
    if(!is.null(target_str)){
        target_traits <- grepl(target_str,row.names(cor_res$P), ignore.case = T)
        pmat <- cor_res$P[,target_traits]
    } else {pmat <- cor_res$P}
    pmat_sig <- pmat[DelayedArray::rowMins(pmat, na.rm = T)<pval_thresh, ]
    rmat_sig <- cor_res$r[rownames(cor_res$r) %in% rownames(pmat_sig) , colnames(pmat_sig)] 
    if(target_str_in_rows==F){
        rmat_sig <- rmat_sig[(!rownames(rmat_sig) %in% colnames(rmat_sig)),]
    }
    return(rmat_sig)
}