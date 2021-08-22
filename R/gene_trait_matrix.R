#' Construct gene x trait matrix from MAGMA results 
#' 
#' Once you've run \link[MAGMA.Celltyping]{map.snps.to.genes} on multiple GWAS, 
#' you can aggregate the gene-wise results into a single gene x trait matrix.
#' 
#' @param magma_files Paths to \emph{genes.out} files
#'  produced by \link[MAGMA.Celltyping]{map.snps.to.genes}.
#' @param metric Which metric to fill the matrix with 
#' \itemize{
#'  \item{\code{"ADJ_ZSTAT"} : }{MAGMA Z-statistic after adjusting for
#'  gene length and other confounds (\emph{DEFAULT}).}
#'  \item{\code{"ZSTAT"} : }{Unadjusted MAGMA Z-statistic.}
#'  \item{\code{"P"} : }{MAGMA P-value.}
#' }
#' @param save_path Path to save results as RDS object.
#' @inheritParams adjust_zstat
#' 
#' @return Sparse matrix 
#' @export
gene_trait_matrix <- function(magma_files,
                              metric=c("ADJ_ZSTAT","ZSTAT","P"),
                              drop_MHC=TRUE,
                              save_path=NULL,
                              mc.cores=1, 
                              fillna=TRUE,
                              agg_FUN="mean"){
    
    metric <- toupper(metric[1])
    fill_value <- if(metric=="P") 1 else 0
    
    #### Merge all results ####
    magma_dt <- parallel::mclapply(names(magma_files), 
                                   function(x,
                                            .metric=metric,
                                            .drop_MHC=drop_MHC){
        message_parallel(x)
        dat <- data.table::fread(magma_files[[x]], nThread = 1)
        #### Adjust ZSTAT ####
        if(.metric=="ADJ_ZSTAT"){
            dat <- adjust_zstat(dat = dat,
                                drop_MHC = .drop_MHC,
                                log_vars=c("NSNPS","NPARAM","GENELEN"),
                                formula=ZSTAT ~ NSNPS + logNSNPS + NPARAM + 
                                      logNPARAM + GENELEN + logGENELEN,
                                verbose = verbose)
        }
        data.table::setkey(dat, "GENE") 
        data.table::setnames(dat, old = .metric, new = x)
        return(dat[,c("GENE",..x)])
    }, mc.cores = mc.cores) %>% 
        base::Reduce(f = function(x,y){merge(x,y,all.x=TRUE,all.y=TRUE)})
    #### Makes colnames compatible with data.table and Matrix format ####
    message("Replacing '-'/':'/'.' with '_' in colnames ",
            "to make compatible with sparse matrix format")
    old_colnames <- colnames(magma_dt)  
    new_colnames <- gsub("-|[:]|[.]","_",colnames(magma_dt))
    data.table::setnames(magma_dt, old_colnames, new_colnames)
    #### Fill NAs ####
    if(fillna){
        message("Replacing NAs with ",fill_value,".")
        #### Very efficient NA replacement  ####  
        na.replace = function(v,value=fill_value) { v[is.na(v)] = value; v }
        for (i in names(magma_dt)){
            eval(parse(text=paste("magma_dt[,",i,":=na.replace(",i,")]")))
        } 
    }
    #### Translate gene IDs ####
    magma_dt <- translate_geneids_magma(magma_dt = magma_dt, 
                                        gene_col = "GENE")
    #### Convert to sparse matrix ####
    magma_matrix <- Matrix::Matrix(as.matrix(magma_dt[,-1], 
                                             rownames = magma_dt$GENE),
                                   sparse = TRUE)
    #### -log transform ####
    if(metric=="P"){
        message("Applying -log10 transformation to P-value matrix.")
        magma_matrix <- Matrix::Matrix(-log10(magma_matrix), sparse=TRUE)
    } 
    
   if(!is.null(agg_FUN) & (sum(duplicated(rownames(magma_matrix)))>0) ){
       message("Aggregating duplicated gene rows by ",agg_FUN)
       magma_matrix <-  orthogene:::aggregate_rows(X = magma_matrix, 
                                          groupings = rownames(magma_matrix), 
                                          FUN = agg_FUN, 
                                          as_DelayedArray = FALSE)
   } 
    #### Save ####
    if(!is.null(save_path)){
        message("Saving results ==> ",save_path)
        saveRDS(magma_matrix, save_path)
    } 
    return(magma_matrix)
}
