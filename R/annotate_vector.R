
annotate_vector <- function(vec, 
                            annot=NULL,
                            gene_level=F){
    if(is.null(annot)){
        message("+ Missing `annot` argument. Returning input `vec`")
        return(vec)
    }
    if(gene_level){
        print("+ Annotating at the gene-level")
        pos_dict <- setNames(annot$Gene_symbol, annot$CHROM_POS) 
    }else {
        print("+ Annotating at the SNP-level")
        pos_dict <- setNames(annot$ID, annot$CHROM_POS)  
    }
    # sum(vec %in% names(pos_dict))
    vec_trans <- pos_dict[vec]
    return(vec_trans)
}