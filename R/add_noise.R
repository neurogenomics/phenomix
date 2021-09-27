#' Add small amounts of random noise to a matrix
#' 
#' @keywords internal
#' @importFrom stats runif
add_noise <- function(data) {
    if (is.vector(data)) {
        noise <- stats::runif(length(data), -0.00001, 0.00001)
        noisified <- data + noise
    } else {
        length <- dim(data)[1] * dim(data)[2]
        noise <- matrix(stats::runif(length, -0.0001, 0.00001), dim(data)[1])
        noisified <- data + noise
    }
    return(noisified)
}
