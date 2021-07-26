get_trait_markers <- function(seurat,
                              trait_set=NULL, 
                              assay1="MAGMA",
                              assay1_slot="scale.data",
                              assay2="specificity",
                              assay2_slot="counts",
                              n_genes=100){
    trait_set <- grep(paste(trait_set,collapse = '|'), colnames(seurat), 
                      ignore.case = T, value = T)
    print(paste(length(trait_set),"traits identified."))
    if(!is.null(assay1)){ 
        expressed_genes <- head( sort(rowMeans(Seurat::GetAssayData(seurat, assay=assay1, slot=assay1_slot)[,trait_set]),
                                      decreasing = T), n_genes)
    }else{expressed_genes <- NULL}
    if(!is.null(assay2)){
        specific_genes <- head( sort(rowMeans(Seurat::GetAssayData(seurat, assay=assay2, slot=assay2_slot)),
                                     decreasing = T), n_genes)
    }else {specific_genes <- NULL} 
    
    
    #### Return marker genes ####
    if(sum(!is.null(assay1), !is.null(assay2))==2){
        marker_genes <- intersect(names(expressed_genes), names(specific_genes))
    } else{
        marker_genes <- if(!is.null(assay1)) expressed_genes else specific_genes
    }
    print(paste(length(marker_genes),"marker genes identified."))
    return(marker_genes)
}