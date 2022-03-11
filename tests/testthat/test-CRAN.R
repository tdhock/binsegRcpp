library(binsegRcpp)
library(testthat)

test_that("one data point has zero loss", {
  L <- binsegRcpp::binseg_normal(5)
  fit <- L$splits
  expect_identical(fit[["loss"]], 0)
  expect_identical(fit[["before.mean"]], 5)
})

sloss <- function(m, x){
  sum(m*(m-2*x)+x^2)
}
test_that("equal split cost is ok", {
  x <- c(0, 0.1, 1, 1.1, 0, 0.1)
  L <- binsegRcpp::binseg_normal(x, length(x))
  fit <- L$splits
  expect_equal(sort(fit$end[1:3]), c(2, 4, 6))
  m1 <- mean(x)
  expect_equal(fit$before.mean[1], m1)
  expect_equal(fit$loss[1], sloss(m1, x))
  expect_equal(fit$loss[6], 0)
  m3 <- c(0.05, 0.05, 1.05, 1.05, 0.05, 0.05)
  expect_equal(fit$loss[3], sloss(m3, x))
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
pos.L <- binseg_normal(x)
pos.dt <- pos.L$splits
neg.L <- binseg_normal(-x)
neg.dt <- neg.L$splits
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
    coef(pos.L, 12.5)
  }, "segments must be a vector of unique integers between 1 and 6")
  expect_error({
    coef(pos.L, -1L)
  }, "segments must be a vector of unique integers between 1 and 6")
  expect_error({
    coef(pos.L, 5:10)
  }, "segments must be a vector of unique integers between 1 and 6")
  segs.vec <- 5:6
  segs.dt <- coef(pos.L, segs.vec)
  expect_equal(nrow(segs.dt), sum(segs.vec))
})

test_that("rcpp_binseg_normal means ok for negative data", {
  expect_equal(pos.dt[["loss"]], neg.dt[["loss"]])
  expect_equal(pos.dt[["end"]], neg.dt[["end"]])
  expect_equal(pos.dt[["before.mean"]], -neg.dt[["before.mean"]])
  expect_equal(pos.dt[["after.mean"]][-1], -neg.dt[["after.mean"]][-1])
  expect_equal(pos.dt[["before.size"]], neg.dt[["before.size"]])
  expect_equal(pos.dt[["after.size"]], neg.dt[["after.size"]])
  expect_equal(pos.dt[["invalidates.index"]], neg.dt[["invalidates.index"]])
  expect_equal(pos.dt[["invalidates.after"]], neg.dt[["invalidates.after"]])
})

test_that("validation loss ok for simple example", {
  is.validation <-
    c(1,0,1,1,    1,  0,     0,  1,  1)
  position <-
    c(1,2,3,4,  101,102,   201,202,203)
  data.vec <-
    c(1,1,1,1,    2,  2,    30, 30, 30)
  kmax <- sum(!is.validation)
  L <- binsegRcpp::binseg_normal(data.vec, kmax, is.validation, position)
  fit <- L$splits
  subtrain.vec <- data.vec[is.validation==0]
  validation.vec <- data.vec[is.validation==1]
  m1 <- mean(subtrain.vec)
  expect_equal(fit$validation.loss[1], sum((validation.vec-m1)^2))
  m2 <- rep(c(1.5, 30), c(4, 2))
  expect_equal(fit$validation.loss[2], sum((validation.vec-m2)^2))
  expect_equal(fit$validation.loss[3], 0)
})

test_that("error for no subtrain data", {
  expect_error({
    binsegRcpp::binseg_normal(5, is.validation.vec=1)
  }, "need at least one subtrain data")
})

test_that("error for two segments with one subtrain", {
  expect_error({
    binsegRcpp::binseg_normal(1:2, 2, is.validation.vec=0:1)
  }, "too many segments")
})

test_that("two data with one subtrain and one segment is ok", {
  L <- binsegRcpp::binseg_normal(1:2, 1, is.validation.vec=0:1)
  fit <- L$splits
  expect_equal(fit$loss, 0)
  expect_equal(fit$validation.loss, 1)
  expect_equal(fit$before.mean, 1)
  L <- binsegRcpp::binseg_normal(2:1, 1, is.validation.vec=0:1)
  fit <- L$splits
  expect_equal(fit$loss, 0)
  expect_equal(fit$validation.loss, 1)
  expect_equal(fit$before.mean, 2)
})

test_that("error for positions not increasing", {
  expect_error({
    binsegRcpp::binseg_normal(1:2, position.vec=2:1)
  }, "positions must increase")
})

test_that("error for NA data", {
  expect_error({
    binsegRcpp::binseg("mean_norm", c(1, 4.3, NA, 5))
  }, "data must be finite")
})

test_that("error for unrecognized distribution", {
  expect_error({
    binsegRcpp::binseg("foo", c(1, 4.3, 5))
  }, "unrecognized distribution")
})
