library(binsegRcpp)
library(testthat)

test_that("one data point has zero loss", {
  L <- binsegRcpp:::binseg_interface(5, weight_vec=1, max_segments=1, min_segment_length=1, distribution_str="mean_norm", container_str = "multiset", is_validation_vec = FALSE, position_vec = 1)
  expect_identical(L[["loss"]], 0)
  expect_identical(L[["before.param.mat"]], cbind(mean=5))
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
  }, "max_segments must be positive")
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
  expect_warning({
    L <- binsegRcpp::binseg_normal(data.vec, kmax, as.logical(is.validation), position)
  }, "some consecutive data values are identical in set=validation")
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
    binsegRcpp::binseg_normal(5, is.validation.vec=TRUE)
  }, "need at least one subtrain data")
})

test_that("error for two segments with one subtrain", {
  expect_error({
    binsegRcpp::binseg_normal(1:2, 2, is.validation.vec=c(FALSE,TRUE))
  }, "too many segments")
})

test_that("two data with one subtrain and one segment is ok", {
  L <- binsegRcpp::binseg_normal(1:2, 1, is.validation.vec=c(FALSE,TRUE))
  fit <- L$splits
  expect_equal(fit$loss, 0)
  expect_equal(fit$validation.loss, 1)
  expect_equal(fit$before.mean, 1)
  L <- binsegRcpp::binseg_normal(2:1, 1, is.validation.vec=c(FALSE,TRUE))
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


ploss <- function (count, seg.mean, weight = 1){
    stopifnot(is.numeric(count))
    stopifnot(is.numeric(seg.mean))
    stopifnot(is.numeric(weight))
    n.data <- length(count)
    if (length(seg.mean) == 1) {
        seg.mean <- rep(seg.mean, n.data)
    }
    if (length(weight) == 1) {
        weight <- rep(weight, n.data)
    }
    stopifnot(n.data == length(seg.mean))
    stopifnot(n.data == length(weight))
    if (any(weight < 0)) {
        stop("PoissonLoss undefined for negative weight")
    }
    if (any(seg.mean < 0)) {
        stop("PoissonLoss undefined for negative segment mean")
    }
    not.integer <- round(count) != count
    not.positive <- count < 0
    loss <- ifelse(not.integer | not.positive, Inf, ifelse(seg.mean == 
        0, ifelse(count == 0, 0, Inf), seg.mean - count * log(seg.mean)))
    sum(loss * weight)
}
test_that("poisson ok with identity weights", {
  data.vec <- c(3, 4, 15, 20)
  fit <- binsegRcpp::binseg("poisson", data.vec)
  expect_equal(fit$splits$end, c(4, 2, 3, 1))
  expected.segs <- list(
    list(data.vec),
    list(c(3,4), c(15, 20)),
    list(c(3,4), 15, 20),
    list(3, 4, 15, 20))
  expected.loss <- sapply(expected.segs, function(L){
    mean.vec <- rep(sapply(L, mean), sapply(L, length))
    ploss(data.vec, mean.vec)
  })
  expect_equal(fit$splits$loss, expected.loss)
})

test_that("poisson ok with non-identity weights", {
  data.vec <- c(3, 4, 15, 20)
  w.each <- 2
  weight.vec <- rep(w.each, length(data.vec))
  fit <- binsegRcpp::binseg("poisson", data.vec, weight.vec=weight.vec)
  expect_equal(fit$splits$end, c(4, 2, 3, 1))
  expected.segs <- list(
    list(data.vec),
    list(c(3,4), c(15, 20)),
    list(c(3,4), 15, 20),
    list(3, 4, 15, 20))
  expected.loss <- sapply(expected.segs, function(L){
    mean.vec <- rep(sapply(L, mean), sapply(L, length))
    ploss(data.vec, mean.vec)*w.each
  })
  expect_equal(fit$splits$loss, expected.loss)
})

test_that("at least one distribution", {
  dist.df <- binsegRcpp::get_distribution_info()
  expect_gt(nrow(dist.df), 0)
})

test_that("min seg length enforced", {
  data.vec <- c(3,4,10,20)
  (fit1 <- binsegRcpp::binseg("poisson", data.vec, weight.vec=c(1,1,1,10)))
  expect_equal(nrow(fit1$splits), 4)
  (fit1.segs2 <- coef(fit1, 2L))
  expect_equal(fit1.segs2$end, c(3,4))
  (fit2 <- binsegRcpp::binseg("poisson", data.vec, weight.vec=c(1,1,1,10), min.segment.length=2L))
  expect_equal(nrow(fit2$splits), 2)
  (fit2.segs2 <- coef(fit2, 2L))
  expect_equal(fit2.segs2$end, c(2,4))
})

test_that("no crash for max_segs < n_data", {
  n.segs <- 5L
  seg.mean.vec <- 1:n.segs
  data.mean.vec <- rep(seg.mean.vec, each=20)
  set.seed(1)
  n.data <- length(data.mean.vec)
  data.vec <- rnorm(n.data, data.mean.vec, 0.2)
  fit <- binsegRcpp::binseg("mean_norm", data.vec, n.segs)
  segs.dt <- coef(fit, n.segs)
  expect_equal(nrow(segs.dt), n.segs)
})

test_that("error for invalid min seg length", {
  expect_error({
    binsegRcpp::binseg("poisson", c(3,4,10,20), min.segment.length=0L)
  }, "min.segment.length must be a positive integer", fixed=TRUE)
})
    
test_that("either two or three segments for max.segments=3", {
  fit2 <- binsegRcpp::binseg("mean_norm", 1:6, min.segment.length=2L, max.segments=3L)
  expect_equal(fit2$splits$end, c(6, 3))
  fit3 <- binsegRcpp::binseg("mean_norm", c(0,0.1, 1,1.1, 0.5,0.6), min.segment.length=2L, max.segments=3L)
  expect_equal(fit3$splits$end, c(6, 2, 4))
})

test_that("error for incompatible max.segments/min.segment.length", {
  expect_error({
    binsegRcpp::binseg("mean_norm", 1:5, min.segment.length=2L, max.segments=3L)
  }, "too many segments")
})

test_that("warning for consecutive data", {
  expect_warning({
    binsegRcpp::binseg("mean_norm", c(1,1))
  }, "some consecutive data values are identical in set=subtrain")
})

test_that("error for unrecognized container", {
  expect_error({
    binsegRcpp::binseg("mean_norm", 1:4, container.str="foo")
  }, "unrecognized container")
})

test_that("variance estimates and loss correct", {
  x <- c(0,0.1, 1,1.2)
  fit <- binsegRcpp::binseg("meanvar_norm", x, max.segments=2L)
  myvar <- function(y)mean((y-mean(y))^2)
  expect_equal(fit$splits$before.mean, c(mean(x), mean(x[1:2])))
  expect_equal(fit$splits$after.mean, c(NA, mean(x[3:4])))
  expect_equal(fit$splits$before.var, c(myvar(x), myvar(x[1:2])))
  expect_equal(fit$splits$after.var, c(NA, myvar(x[3:4])))
  nll <- function(y)-sum(dnorm(y, mean(y), sqrt(myvar(y)), log=TRUE))
  expect_equal(fit$splits$loss, c(nll(x), nll(x[1:2])+nll(x[3:4])))
  seg.dt <- coef(fit, 2L)
  expect_equal(seg.dt$end, c(2,4))
  expect_equal(seg.dt$mean, c(mean(x[1:2]),mean(x[3:4])))
  expect_equal(seg.dt$var, c(myvar(x[1:2]),myvar(x[3:4])))
})

test_that("meanvar_norm does not have segs with size 1", {
  data.per.seg <- 10
  sim <- function(mu,sigma)rnorm(data.per.seg,mu,sigma)
  set.seed(1)
  data.vec <- c(sim(10,1), sim(0, 5))
  fit <- binsegRcpp::binseg("meanvar_norm", data.vec)
  expect_lte(nrow(fit$splits), data.per.seg)
})

test_that("l1loss param is median", {
  data.vec <- c(1.3, 1.0, 1.1, 2.0, 2.1, 3.1)
  l1fit <- binsegRcpp::binseg("l1", data.vec, max.segments=2L)
  seg.dt <- coef(l1fit)
  expected.median <- c(
    median(data.vec),
    median(data.vec[1:3]),
    median(data.vec[4:6]))
  expect_equal(seg.dt$median, expected.median)
  expected.loss <- sum(abs(median(data.vec)-data.vec))
  expect_equal(l1fit$splits$loss[1], expected.loss)
})

test_that("laplace params median,scale", {
  data.vec <- c(1.3, 1.0, 1.1, 2.0, 2.1, 3.1)
  l1fit <- binsegRcpp::binseg("laplace", data.vec, max.segments=2L)
  seg.dt <- coef(l1fit)
  expected.median <- c(
    median(data.vec),
    median(data.vec[1:3]),
    median(data.vec[4:6]))
  expect_equal(seg.dt$median, expected.median)
  sum.abs.dev <- sum(abs(median(data.vec)-data.vec))
  N.data <- length(data.vec)
  est.scale <- sum.abs.dev/N.data
  expect_equal(l1fit$splits$before.scale[1], est.scale)
  expected.loss <- N.data*log(2*est.scale)+sum.abs.dev/est.scale
  expect_equal(l1fit$splits$loss[1], expected.loss)
})
