replace_char_fun <- function(nms,
                             replace_char=list("."=":")){
    if(length(replace_char)>0){
        for(r in names(replace_char)){
            nms <- gsub(r,replace_char[[r]],nms, fixed = TRUE)
        }
    }
    return(nms)
}