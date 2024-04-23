## three tests which are the same except for DIST arg.
partial.list <- lapply(atime::atime_grid(
  list(DIST=c("l1", "poisson", "meanvar_norm")),
  "binseg(1:N,maxSegs=N/2)"=atime::atime_test(
    expr=binsegRcpp::binseg(DIST, data.vec, max.segs))),
  eval)
test.list <- atime::atime_test_list(
  N=2^seq(2, 20),
  setup={
    max.segs <- as.integer(N/2)
    data.vec <- 1:N
  },
  tests=partial.list,
  ## A test using a different expr, limit, historical commit.
  "binseg_normal(1:N,maxSegs=N/2)"=atime::atime_test(
    expr=binsegRcpp::binseg_normal(data.vec, max.segs),
    seconds.limit=0.1,
    Before="083c6827aeae00568e31614bbec9bbb31dca1081"),
  NULL)
