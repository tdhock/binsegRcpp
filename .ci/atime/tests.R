test.list <- atime::atime_test_list(
  binseg_normal_best=atime::atime_test(
    setup={
      max.segs <- as.integer(N/2)
      data_vec <- 1:N
    },
    expr=binsegRcpp::binseg_normal(data_vec, max.segs)
  )
)
