NULL [![License: GPL
(\>=3)](https://img.shields.io/badge/license-GPL%20(%3E=3)-blue.svg)](https://cran.r-project.org/web/licenses/GPL%20(%3E=3))
[![](https://img.shields.io/badge/devel%20version-0.99.5-black.svg)](https://github.com/neurogenomics/phenomix)
[![](https://img.shields.io/github/languages/code-size/neurogenomics/phenomix.svg)](https://github.com/neurogenomics/phenomix)
[![](https://img.shields.io/github/last-commit/neurogenomics/phenomix.svg)](https://github.com/neurogenomics/phenomix/commits/master)
<br> [![R build
status](https://github.com/neurogenomics/phenomix/workflows/rworkflows/badge.svg)](https://github.com/neurogenomics/phenomix/actions)
[![](https://codecov.io/gh/neurogenomics/phenomix/branch/master/graph/badge.svg)](https://app.codecov.io/gh/neurogenomics/phenomix)
<br>
<a href='https://app.codecov.io/gh/neurogenomics/phenomix/tree/master' target='_blank'><img src='https://codecov.io/gh/neurogenomics/phenomix/branch/master/graphs/icicle.svg' title='Codecov icicle graph' width='200' height='50' style='vertical-align: top;'></a>  
<h4>  
Authors: <i>Brian Schilder</i>  
</h4>
<h4>  
Most recent update: <i>Aug-17-2023</i>  
</h4>

## R package for the exploration and analysis of many genotype-phenotype datasets at once.

## Installation

``` r
if(require("remotes")} install.packages("remotes")
remotes::install_github("neurogenomics/phenomix")
```

## Documentation

### [Documentation site](https://neurogenomics.github.io/phenomix)

### [Vignette: Celltype enrichment](https://neurogenomics.github.io/phenomix/articles/celltype_enrichment.html)

# Session Info

<details>

``` r
utils::sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] here_1.0.1          rprojroot_2.0.3     digest_0.6.31      
    ##  [4] utf8_1.2.3          BiocFileCache_2.6.1 R6_2.5.1           
    ##  [7] stats4_4.2.1        RSQLite_2.3.1       evaluate_0.21      
    ## [10] httr_1.4.6          ggplot2_3.4.2       pillar_1.9.0       
    ## [13] yulab.utils_0.0.6   rworkflows_0.99.10  biocViews_1.66.3   
    ## [16] rlang_1.1.1         curl_5.0.0          data.table_1.14.8  
    ## [19] rstudioapi_0.14     whisker_0.4.1       blob_1.2.4         
    ## [22] DT_0.28             RUnit_0.4.32        rmarkdown_2.22     
    ## [25] desc_1.4.2          readr_2.1.4         stringr_1.5.0      
    ## [28] htmlwidgets_1.6.2   dlstats_0.1.7       BiocPkgTools_1.16.1
    ## [31] igraph_1.5.0.1      RCurl_1.98-1.12     bit_4.0.5          
    ## [34] munsell_0.5.0       compiler_4.2.1      xfun_0.39          
    ## [37] pkgconfig_2.0.3     BiocGenerics_0.44.0 rorcid_0.7.0       
    ## [40] htmltools_0.5.5     tidyselect_1.2.0    tibble_3.2.1       
    ## [43] httpcode_0.3.0      XML_3.99-0.14       fansi_1.0.4        
    ## [46] dplyr_1.1.2         tzdb_0.4.0          dbplyr_2.3.2       
    ## [49] bitops_1.0-7        rappdirs_0.3.3      crul_1.4.0         
    ## [52] grid_4.2.1          RBGL_1.74.0         jsonlite_1.8.4     
    ## [55] gtable_0.3.3        lifecycle_1.0.3     DBI_1.1.3          
    ## [58] magrittr_2.0.3      scales_1.2.1        graph_1.76.0       
    ## [61] cli_3.6.1           stringi_1.7.12      cachem_1.0.8       
    ## [64] renv_0.17.3         fauxpas_0.5.2       xml2_1.3.4         
    ## [67] rvcheck_0.2.1       filelock_1.0.2      generics_0.1.3     
    ## [70] vctrs_0.6.2         gh_1.4.0            RColorBrewer_1.1-3 
    ## [73] tools_4.2.1         bit64_4.0.5         Biobase_2.58.0     
    ## [76] glue_1.6.2          hms_1.1.3           fastmap_1.1.1      
    ## [79] yaml_2.3.7          colorspace_2.1-0    BiocManager_1.30.20
    ## [82] rvest_1.0.3         memoise_2.0.1       badger_0.2.3       
    ## [85] knitr_1.43

</details>
