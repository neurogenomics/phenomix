run_phate <- function(){
    echo
    echoconda::activate_env("phate")
    phateR::install.phate()
    reticulate::py_module_available("phate")
}