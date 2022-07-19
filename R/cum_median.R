cum_median <- function
### Efficient log-linear cumulative median.
(data.vec,
### Numeric vector of data.
  weight.vec=rep(1, length(data.vec))
### Numeric vector of weights.
){
  cum_median_interface(data.vec, weight.vec)
}
