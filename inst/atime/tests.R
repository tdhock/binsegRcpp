test.list <- c(
  atime::atime_grid(
    list(DIST=c("l1", "poisson", "meanvar_norm")),
    "binseg(1:N,maxSegs=N/2)"=list(
      N=2^seq(2, 20),
      setup={
        max.segs <- as.integer(N/2)
        data.vec <- 1:N
      },
      expr=binsegRcpp::binseg(DIST, data.vec, max.segs))),
  atime::atime_grid(
    "binseg_normal(1:N,maxSegs=N/2)"=list(
      seconds.limit=0.1,
      N=2^seq(2, 20),
      setup={
        max.segs <- as.integer(N/2)
        data.vec <- 1:N
      },
      expr=binsegRcpp::binseg_normal(data.vec, max.segs))))

