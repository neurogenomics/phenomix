

add_noise <- function(data) {
    if (is.vector(data)) {
        noise <- runif(length(data), -0.00001, 0.00001)
        noisified <- data + noise
    } else {
        length <- dim(data)[1] * dim(data)[2]
        noise <- matrix(runif(length, -0.0001, 0.00001), dim(data)[1])
        noisified <- data + noise
    }
    return(noisified)
}
