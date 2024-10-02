# phenomix 0.99.6

## New features

* Transfer `phenoRx::prepare_hpo` to `phenomix`
* Remove all ontology-related functions (now within `HPOExplorer`)
* Explore `bigstatsr`/`bigreadr` a bit.
* New functions, introduced while writing up thesis:
    - `clip_outliers`
    - `run_cor`
    - `run_knn_overlap`
    - `un_ontological_similarity`
    - `run_phate`
    - `run_preseveration`
    - `run_pseudotime_subtypes`
    - `seurat_to_ontological_similarity`

## Bug fixes

* Fix GH token in rworkflows yml

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
* `gwas_matrix_magma` -> `magma_matrix`
    - Replace `parallel` with `BiocParallel`
* Change args globally:
    - `reduction` -> `keys`
    - `drop_MHC` -> `drop_mhc`
* Remove functions:
    - `create_DT`
* New functions:
    - `phenomix_query`
    - `phenomix_matrix`
    - `map_phenotypes`
* New data: `OpenGWAS`
    - Supplied/prepared by `opengwas_meta`
* Condense `get_cs2g_ukb`/`get_cs2g_gwascatalog` into `get_cs2g`.
* Convert to using `data.table`
    - Remove all `reshape2` functions.
    - Remove some `dplyr` functions.
* New function: `get_ctd`
    
## Bug fixes

* Get rid of namespace conflicts: `Warning messages:`
    - `1: replacing previous import ‘data.table::melt’ by ‘reshape2::melt’ when loading ‘phenomix’`
    - `2: replacing previous import ‘Matrix::head’ by ‘utils::head’ when loading ‘phenomix’`

# phenomix 0.99.4

## New features

* Added a `NEWS.md` file to track changes to the package.
* Updated GHA, added *Dockerfile*, removed *docs* folder. 

## Bug fixes 

* `gwas_matrix_magma`: general improvements. 
