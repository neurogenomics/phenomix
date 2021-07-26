

setup_spark <- function(){
    # https://community.rstudio.com/t/how-to-read-large-json-file-in-r/13486/32
    Sys.setenv(SPARK_HOME="/usr/lib/spark")
    # Configure cluster (c3.4xlarge 30G 16core 320disk)
    conf <- sparklyr::spark_config()
    conf$'sparklyr.shell.executor-memory' <- "7g"
    conf$'sparklyr.shell.driver-memory' <- "7g"
    conf$spark.executor.cores <- 20
    conf$spark.executor.memory <- "7G"
    conf$spark.yarn.am.cores  <- 20
    conf$spark.yarn.am.memory <- "7G"
    conf$spark.executor.instances <- 20
    conf$spark.dynamicAllocation.enabled <- "false"
    conf$maximizeResourceAllocation <- "true"
    conf$spark.default.parallelism <- 32
    sc <- sparklyr::spark_connect(master = "local", config = conf, version = '2.2.0')
    return(sc)
}