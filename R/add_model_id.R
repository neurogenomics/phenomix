add_model_id <- function(res,i){ 
    mid <- gsub("file",paste0("model",i,"_"),basename(tempfile()))
    res[,model_id:=mid,]  
}