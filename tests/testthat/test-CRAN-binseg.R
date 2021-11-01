library(binsegRcpp)
library(testthat)

test_that("one data point has zero loss", {
  fit <- binsegRcpp::binseg_normal(5)
  expect_identical(fit[["loss"]], 0)
  expect_identical(fit[["before.mean"]], 5)
})

sloss <- function(m, x){
  sum(m*(m-2*x)+x^2)
}
test_that("equal split cost is ok", {
  x <- c(0, 0.1, 1, 1.1, 0, 0.1)
  L <- binsegRcpp::binseg_normal(x, length(x))
  expect_equal(sort(L$end[1:3]), c(2, 4, 6))
  m1 <- mean(x)
  expect_equal(L$before.mean[1], m1)
  expect_equal(L$loss[1], sloss(m1, x))
  expect_equal(L$loss[6], 0)
  m3 <- c(0.05, 0.05, 1.05, 1.05, 0.05, 0.05)
  expect_equal(L$loss[3], sloss(m3, x))
})

test_that("error for 0 data", {
  x <- double()
  expect_error({
    binsegRcpp::binseg_normal(x, 5L)
  }, "need at least one data point")
})

test_that("error for 0 segments", {
  x <- c(4.1, 4, 1.1, 1)
  expect_error({
    binsegRcpp::binseg_normal(x, 0L)
  }, "kmax must be positive")
})

test_that("error for too many segments", {
  x <- c(4.1, 4, 1.1, 1)
  expect_error({
    binsegRcpp::binseg_normal(x, 10L)
  }, "too many segments")
})

x <- c(0.2, 0, 1, 1.4, 3.6, 3)
pos.dt <- binseg_normal(x)
neg.dt <- binseg_normal(-x)
test_that("binseg_normal means ok for negative data", {
  expect_equal(pos.dt[["loss"]], neg.dt[["loss"]])
  expect_equal(pos.dt[["end"]], neg.dt[["end"]])
  expect_equal(pos.dt[["before.mean"]], -neg.dt[["before.mean"]])
  expect_equal(pos.dt[["after.mean"]][-1], -neg.dt[["after.mean"]][-1])
  expect_equal(pos.dt[["before.size"]], neg.dt[["before.size"]])
  expect_equal(pos.dt[["after.size"]], neg.dt[["after.size"]])
  expect_equal(pos.dt[["invalidates.index"]], neg.dt[["invalidates.index"]])
  expect_equal(pos.dt[["invalidates.after"]], neg.dt[["invalidates.after"]])
})

test_that("error for invalid coef segments", {
  expect_error({
    coef(pos.dt, 12.5)
  }, "segments must be a vector of unique integers between 1 and 6")
  expect_error({
    coef(pos.dt, -1L)
  }, "segments must be a vector of unique integers between 1 and 6")
  expect_error({
    coef(pos.dt, 5:10)
  }, "segments must be a vector of unique integers between 1 and 6")
  segs.vec <- 5:6
  segs.dt <- coef(pos.dt, segs.vec)
  expect_equal(nrow(segs.dt), sum(segs.vec))
})

pos.list <- binsegRcpp::binseg_normal(x, length(x))
neg.list <- binsegRcpp::binseg_normal(-x, length(x))
test_that("rcpp_binseg_normal means ok for negative data", {
  expect_equal(pos.list[["loss"]], neg.list[["loss"]])
  expect_equal(pos.list[["end"]], neg.list[["end"]])
  expect_equal(pos.list[["before.mean"]], -neg.list[["before.mean"]])
  expect_equal(pos.list[["after.mean"]][-1], -neg.list[["after.mean"]][-1])
  expect_equal(pos.list[["before.size"]], neg.list[["before.size"]])
  expect_equal(pos.list[["after.size"]], neg.list[["after.size"]])
  expect_equal(pos.list[["invalidates.index"]], neg.list[["invalidates.index"]])
  expect_equal(pos.list[["invalidates.after"]], neg.list[["invalidates.after"]])
})
