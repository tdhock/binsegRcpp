library(binsegRcpp)
library(testthat)
context("rcpp_binseg_normal")

test_that("equal split cost is ok", {
  x <- c(0, 0.1, 1, 1.1, 0, 0.1)
  L <- rcpp_binseg_normal(x, length(x))
  expect_equal(sort(L$end[1:3]), c(1, 3, 5))
  m <- mean(x)
  expect_equal(L$before.mean[1], m)
  const.term <- -sum(x^2)
  ## cost does not include the constant x^2 term.
  expect_equal(L$loss[1], sum(m*(m-2*x)))
  expect_equal(L$loss[6], const.term)
  m <- c(0.05, 0.05, 1.05, 1.05, 0.05, 0.05)
  expect_equal(L$loss[3], sum(m*(m-2*x)))
})

test_that("error for 0 data", {
  x <- double()
  expect_error({
    rcpp_binseg_normal(x, 5L)
  }, "no data")
})

test_that("error for 0 segments", {
  x <- c(4.1, 4, 1.1, 1)
  expect_error({
    rcpp_binseg_normal(x, 0L)
  }, "kmax must be positive")
})

test_that("error for too many segments", {
  x <- c(4.1, 4, 1.1, 1)
  expect_error({
    rcpp_binseg_normal(x, 10L)
  }, "too many segments")
})

test_that("means ok for negative data", {
  x <- c(0.2, 0, 1, 1.4, 3.6, 3)
  pos.dt <- binseg_normal(x)
  neg.dt <- binseg_normal(-x)
  expect_equal(pos.dt[["loss"]], neg.dt[["loss"]])
  expect_equal(pos.dt[["end"]], neg.dt[["end"]])
  expect_equal(pos.dt[["before.mean"]], -neg.dt[["before.mean"]])
  expect_equal(pos.dt[["after.mean"]], -neg.dt[["after.mean"]])
  expect_equal(pos.dt[["before.size"]], neg.dt[["before.size"]])
  expect_equal(pos.dt[["after.size"]], neg.dt[["after.size"]])
  expect_equal(pos.dt[["invalidates.index"]], neg.dt[["invalidates.index"]])
  expect_equal(pos.dt[["invalidates.after"]], neg.dt[["invalidates.after"]])
})
