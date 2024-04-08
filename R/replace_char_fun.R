#' Map ID separators 
#' 
#' Replace the separator in an ID in order to standardise it.
#' @export
#' @examples
#' map_id_sep(c("HP:0001","HP.0001","HP_0001"))
map_id_sep <- function(nms,
                             replace_char=list("."=":",
                                               "_"=":")){
    if(length(replace_char)>0){
        for(r in names(replace_char)){
            nms <- gsub(r,replace_char[[r]],nms, fixed = TRUE)
        }
    }
    return(nms)
}
