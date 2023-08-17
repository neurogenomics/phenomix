# phenomix 0.99.5

## New features

* Add `rworkflows`
* Offload element-extraction functions to `scKirby`
    - `extract_metadata` -> `get_obs`
    - `extract_embeddings` -> `get_obsm`
    - `extract_loadings` -> `get_varm`
    - `extract_graph` -> `get_graphs`
    - `extract_matrix` -> `get_x`
    - `extract_colnames` -> `get_obs_names`
    - `is_seurat` -> `is_class("seurat")`
    - `is_matrix` -> `is_class("matrix")` 
* `extract_cor` -> `get_cor`
* Change args globally:
    - `reduction` -> `keys`
    - `drop_MHC` -> `drop_mhc`
* Remove functions:
    - `create_DT`

# phenomix 0.99.4

## New features

* Added a `NEWS.md` file to track changes to the package.
* Updated GHA, added *Dockerfile*, removed *docs* folder. 


## Bug fixes 

* `gwas_matrix_magma`: general improvements. 

