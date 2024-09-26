N.tests.preview <- 2
N.setup <- atime::atime_grid(
  N=2^seq(2, 20),
  setup={
    max.segs <- as.integer(N/2)
    data.vec <- 1:N
  })
expr.list <- atime::atime_grid(
  list(DIST=c("l1", "poisson", "meanvar_norm")),
  "binseg(1:N,maxSegs=N/2)"=binsegRcpp::binseg(DIST, data.vec, max.segs))
partial.list <- lapply(expr.list, function(expr)list(expr=expr))
partial.list[["binseg_normal(1:N,maxSegs=N/2)"]] <- list(
  expr=quote(binsegRcpp::binseg_normal(data.vec, max.segs)),
  seconds.limit=0.1,
  Before="15becb3f0d730da7a1e0e71ad4c3539ee7858d73")
test.list <- lapply(partial.list, function(plist)c(plist, N.setup))
