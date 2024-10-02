#' Run integration
#' 
#' Run integration pipeline on a \link{Seurat} object.
#' @inheritParams Seurat::IntegrateLayers
#' @inheritParams Seurat::SplitObject
#' @inheritParams Seurat::ScaleData
#' @inheritParams Seurat::CCAIntegration
#' @inheritParams scKirby::process_seurat
#' @inheritDotParams Seurat::IntegrateLayers
#' @export
#' @examples
#' obj <- scKirby::example_obj("seurat")
#' obj2 <- run_integration(obj, split.by="groups")
run_integration <- function(obj, 
                            split.by,
                            assay = Seurat::DefaultAssay(obj),
                            nfeatures = nrow(obj),
                            vars.to.regress = NULL,
                            dims = seq(100),
                            save_path = tempfile(fileext = ".rds"),
                            pipeline = c("seuratv5","seuratv4"),
                            method = Seurat::CCAIntegration, 
                            k.weight = 100,
                            orig.reduction = "pca", 
                            new.reduction = "cca",
                            cluster_reduction = "umap",
                            max_mem = 8000*1024^2,
                            workers=1,
                            seed=2020,
                            verbose = TRUE,
                            ...){
    # SeuratWrappers::scVIIntegration()
    # SeuratWrappers::FastMNNIntegration()
    
    
    pipeline <- match.arg(pipeline)
    if(file.exists(save_path)){
        message("Loading cached file: ",save_path)
        obj <- readRDS(save_path)
    } else {
        requireNamespace("future")
        set.seed(seed) 
        future::plan(strategy = "multicore", workers = workers)
        options(future.globals.maxSize = max_mem)
        force(obj)
        force(split.by)
        ## remove traits with no gene associations 
        # obj <- obj[,obj$nFeature_score!=0]
        
        # obj_hpo <- Seurat::FindVariableFeatures(obj_hpo, nfeatures = nrow(obj_hpo)) 
        # var_features <- lapply(Seurat::SplitObject(obj_ot, split.by = "orig.ident"),
        #                        function(x){
        #   cat(dim(x),"\n")
        #   if(ncol(x)<2) return(NULL)
        #   x <- Seurat::FindVariableFeatures(x, nfeatures = nrow(x))
        #   Seurat::VariableFeatures(x)
        # })
        # ### Get union of top 2000 genes from each gene lists ###
        # var_features <- Reduce(union, lapply(c(var_features,
        #                                         HPO=Seurat::VariableFeatures(obj_hpo)
        #                                        ),
        #                                      head, 1000)) 
       
        max_dim <- min(
            data.table::data.table(obj@meta.data)[,.N,by=c(split.by)]$N
            ) - 1 
        if(max(dims)>max_dim){
            messager("Reducing number of dimensions to match split data:",
                     max_dim)
            dims <- dims[seq(max_dim)]
        }
        #### Seurat v5 strategy ####
        if(pipeline=="seuratv5"){ 
            obj[[assay]] <- split(obj[[assay]],
                                  f = obj@meta.data[[split.by]])
            obj <- Seurat::NormalizeData(obj)
            obj <- Seurat::FindVariableFeatures(obj, 
                                                nfeatures = nfeatures)
            obj <- Seurat::ScaleData(obj, 
                                     vars.to.regress = vars.to.regress)
            if(!orig.reduction %in% names(obj@reductions)){
                if(orig.reduction=="pca"){
                    messager(
                        "PCA not found in merged object.",
                        "Running PCA within each object split.by group now.")
                    obj <- Seurat::RunPCA(obj, 
                                          npcs=max(dims))
                } else{
                    stop("orig.reduction not found in object.")
                }
            }
            obj <- Seurat::IntegrateLayers(object = obj, 
                                           method = method, 
                                           orig.reduction = orig.reduction, 
                                           new.reduction = new.reduction,
                                           features = rownames(obj),
                                           dims = dims,
                                           k.weight = min(k.weight,max_dim),
                                           verbose = verbose)
            # re-join layers after integration
            obj[[assay]] <- SeuratObject::JoinLayers(obj[[assay]]) 
            obj <- Seurat::FindNeighbors(obj, 
                                         reduction = new.reduction,
                                         dims = dims)
            obj <- Seurat::FindClusters(obj)
            max_dims <- ncol(obj@reductions[[new.reduction]])
            if(max(dims)>max_dims){
                messager("Reducing number of dimensions to match reduced data:",
                         max_dims)
                dims <- dims[seq(max_dims)]
            }
            obj <- Seurat::RunUMAP(obj, 
                                   reduction = new.reduction, 
                                   dims = dims, 
                                   return.model = TRUE)
            if(cluster_reduction!=new.reduction){
                max_dims <- ncol(obj@reductions[[cluster_reduction]])
                obj <- Seurat::FindNeighbors(obj, 
                                             reduction = cluster_reduction,
                                             dims = seq(max_dims))
                obj <- Seurat::FindClusters(obj)
            }
        } else if(pipeline=="seuratv4"){ 
            #### Seurat <v5 strategy ####
            ## Merge then split to ensure the union of all genes is used during integration 
            object.list <- Seurat::SplitObject(obj,
                                               split.by = split.by)
            object.list <- lapply(object.list,
                                  Seurat::FindVariableFeatures,
                                  nfeatures = nfeatures)
            anchor.features <- Reduce(f = union,
                                      lapply(object.list, 
                                             Seurat::VariableFeatures) 
                                      )
            #### Integrate data with Seurat ####
            anchors <- Seurat::FindIntegrationAnchors( 
                object.list = object.list,
                dims = dims,
                anchor.features = anchor.features
            )
            obj <- Seurat::IntegrateData(
                anchorset = anchors,
                dims = dims,
                features = anchors@anchor.features,
                features.to.integrate = anchors@anchor.features
            )
            obj <- Seurat::FindVariableFeatures(obj,  nfeatures = nrow(obj)) 
            #### Process Seurat obj ####
            obj <- scKirby::process_seurat(obj,
                                           reduction = new.reduction,
                                           cluster_reduction = cluster_reduction,
                                           ## Regressing out source/orig.ident 
                                           ## seems to give a MUCH worse integration in UMAP space. 
                                           vars.to.regress = vars.to.regress,
                                           nfeatures = NULL,
                                           dims = dims,#setdiff(dims,omit_dims),
                                           normalize_data = FALSE) 
        }  
        #### Save obj ####
        KGExplorer::cache_save(obj = obj, 
                               save_path = save_path) 
    }
    return(obj)
}