test.list <- atime::atime_grid(
  list(DIST=c("mean_norm", "poisson", "meanvar_norm", "l1")),
  "binseg(1:N,maxSegs=N/2)"=list(
    N=2^seq(2, 20),
    setup={
      max.segs <- as.integer(N/2)
      data.vec <- 1:N
    },
    expr=binsegRcpp::binseg(DIST, data.vec, max.segs)))

