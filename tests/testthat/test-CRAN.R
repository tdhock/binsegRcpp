library(binsegRcpp)
library(testthat)
myvar <- function(y)mean((y-mean(y))^2)

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
  }, "number of data must be at least min segment length")
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
  subtrain.pos <- position[is.validation==0]
  validation.vec <- data.vec[is.validation==1]
  m1 <- mean(subtrain.vec)
  expect_equal(fit$validation.loss[1], sum((validation.vec-m1)^2))
  m2 <- rep(c(1.5, 30), c(4, 2))
  expect_equal(fit$validation.loss[2], sum((validation.vec-m2)^2))
  expect_equal(fit$validation.loss[3], 0)
  expected.borders <- c(0.5, 52.5, 151.5, 203.5)
  expect_equal(L$subtrain.borders, expected.borders)
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

test_that("one data point, subtrain borders", {
  L <- binsegRcpp::binseg_normal(1)
  expect_equal(L$subtrain.borders, c(0.5, 1.5))
})

test_that("two data, subtrain borders", {
  L <- binsegRcpp::binseg_normal(1:2)
  expect_equal(L$subtrain.borders, c(0.5, 1.5, 2.5))
})

test_that("two data, one valid, subtrain borders", {
  L <- binsegRcpp::binseg_normal(1:2, is.validation.vec = c(TRUE,FALSE))
  expect_equal(L$subtrain.borders, c(0.5, 2.5))
})

test_that("three data, middle valid, subtrain borders", {
  L <- binsegRcpp::binseg_normal(1:3, is.validation.vec = c(FALSE,TRUE,FALSE))
  expect_equal(L$subtrain.borders, c(0.5, 2.5, 3.5))
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

container.values <- c("priority_queue", "multiset", "list")
test_that("variance estimates and loss correct", {
  x <- c(0,0.1, 1,1.2)
  for(container in container.values){
    fit <- binsegRcpp::binseg("meanvar_norm", x, max.segments=2L, container.str = container)
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
  }
})

test_that("meanvar_norm does not have segs with size 1", {
  data.per.seg <- 10
  sim <- function(mu,sigma)rnorm(data.per.seg,mu,sigma)
  set.seed(1)
  data.vec <- c(sim(10,1), sim(0, 5))
  for(container in container.values){
    fit <- binsegRcpp::binseg("meanvar_norm", data.vec, container.str=container)
    expect_lte(nrow(fit$splits), data.per.seg)
  }
})

test_that("l1loss param is median", {
  data.vec <- c(1.3, 1.0, 1.1, 2.0, 2.1, 3.1)
  for(container in container.values){
    l1fit <- binsegRcpp::binseg("l1", data.vec, max.segments=2L, container.str=container)
    seg.dt <- coef(l1fit)
    expected.median <- c(
      median(data.vec),
      median(data.vec[1:3]),
      median(data.vec[4:6]))
    expect_equal(seg.dt$median, expected.median)
    expected.loss <- sum(abs(median(data.vec)-data.vec))
    expect_equal(l1fit$splits$loss[1], expected.loss)
  }
})

test_that("laplace params median,scale", {
  data.vec <- c(1.3, 1.0, 1.1, 2.0, 2.1, 3.1)
  for(container in container.values){
    l1fit <- binsegRcpp::binseg("laplace", data.vec, max.segments=2L, container.str=container)
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
  }
})

test_that("laplace validation loss ok", {
  data.vec <- c(1.3, 1.0, 1.1, 2.0, 2.1, 3.1)
  is.validation.vec <- c(FALSE,FALSE,TRUE,TRUE,FALSE,FALSE)
  for(container in container.values){
    l1fit <- binsegRcpp::binseg(
      "laplace", data.vec, max.segments=2L,
      container.str=container,
      is.validation.vec = is.validation.vec)
    subtrain.vec <- data.vec[!is.validation.vec]
    est.median <- median(subtrain.vec)
    est.scale <- sum(abs(est.median-subtrain.vec))/length(subtrain.vec)
    sum.abs.dev <- sum(abs(est.median-data.vec[is.validation.vec]))
    vloss1 <- sum(is.validation.vec)*log(2*est.scale)+sum.abs.dev/est.scale
    expect_equal(l1fit$splits$validation.loss[1], vloss1)
  }
})

test_that("laplace correct split for int data", {
  laplace <- function(N){
    y=1:N
    m=median(y)
    d=abs(y-m)
    b=mean(d)
    sum(d/b+log(2*b))
  }
  for(container in container.values){
    fit <- binseg("laplace", 1:8, container.str=container)
    splits <- fit[["splits"]]
    computed.loss <- splits[["loss"]]
    expected.loss <- c(
      laplace(8),
      laplace(3)+laplace(5),
      laplace(3)+laplace(3)+laplace(2))
    expect_equal(computed.loss, expected.loss)
  }
})

test_that("meanvar_norm validation loss ok", {
  data.vec <- c(1.3, 1.0, 1.1, 2.0, 2.1, 3.1)
  is.validation.vec <- c(FALSE,FALSE,TRUE,TRUE,FALSE,FALSE)
  for(container in container.values){
    l1fit <- binsegRcpp::binseg(
      "meanvar_norm", data.vec, max.segments=2L,
      container.str=container,
      is.validation.vec = is.validation.vec)
    subtrain.vec <- data.vec[!is.validation.vec]
    vloss1 <- -sum(dnorm(
      data.vec[is.validation.vec],
      mean(subtrain.vec),
      sqrt(myvar(subtrain.vec)),
      log=TRUE))
    expect_equal(l1fit$splits$validation.loss[1], vloss1)
  }
})

test_that("l1 validation loss ok", {
  data.vec <- c(1.3, 1.0, 1.1, 2.0, 2.1, 3.1)
  is.validation.vec <- c(FALSE,FALSE,TRUE,TRUE,FALSE,FALSE)
  for(container in container.values){
    l1fit <- binsegRcpp::binseg(
      "l1", data.vec, max.segments=2L,
      container.str=container,
      is.validation.vec = is.validation.vec)
    vloss1 <- sum(abs(
      data.vec[is.validation.vec]-median(data.vec[!is.validation.vec])))
    expect_equal(l1fit$splits$validation.loss[1], vloss1)
  }
})

test_that("poisson loss ok for simple ex with zero", {
  data.vec <- 0:1
  N.data <- length(data.vec)
  mu <- mean(data.vec)
  expected.loss <- c(
    N.data*mu - log(mu)*sum(data.vec),
    1)
  for(container in container.values){
    fit <- binsegRcpp::binseg("poisson", data.vec, container.str=container)
    expect_equal(fit$splits$end, 2:1)
    expect_equal(fit$splits$loss, expected.loss)
    segs <- coef(fit)
    expect_equal(segs$mean, c(mu, data.vec))
  }
})

test_that("error for poisson loss with bad data", {
  expect_error({
    binsegRcpp::binseg("poisson", 0.1)
  }, "data must be integer for poisson loss")
  expect_error({
    binsegRcpp::binseg("poisson", -2L)
  }, "data must be non-negative for poisson loss")
})

test_that("get_complexity respects min.segment.length", {
  n.data <- 8
  zero.one <- rep(0:1, l=n.data)
  for(container in container.values){
    fit <- binsegRcpp::binseg("mean_norm", zero.one, container.str=container)
    clist <- binsegRcpp::get_complexity(fit)
    worst <- clist$iterations[case=="worst"]
    expect_equal(worst$splits, seq(n.data-1, 0))
    zero.ten <- rep(c(0,1,10,11), l=n.data)
    mvfit <- binsegRcpp::binseg("meanvar_norm", zero.ten)
    mvlist <- binsegRcpp::get_complexity(mvfit)
    mvworst <- mvlist$iterations[case=="worst"]
    expect_equal(mvworst$splits, c(5,3,1,0))
  }
})

test_that("empirical splits not negative", {
  for(container in container.values){
    fit <- binsegRcpp::binseg("meanvar_norm", 1:8, container.str=container)
    clist <- binsegRcpp::get_complexity(fit)
    esplits <- clist$iterations[case=="empirical", splits]
    expect_equal(esplits, c(5,2,0,0))
  }
})

test_that("first three splits result in four segments, two data each", {
  for(container in container.values){
    fit <- binsegRcpp::binseg("mean_norm", 1:8, max.segments=4, container.str=container)
    expect_equal(sort(fit$splits$end), c(2,4,6,8))
  }
})

test_that("l1 loss chooses even split if equal loss", {
  N.max <- 8
  data.vec <- 1:N.max
  seg <- function(first,last){
    sdata <- data.vec[first:last]
    sum(abs(median(sdata)-sdata))
  }
  two.segs <- function(end){
    sum(c(seg(1,end),seg(end+1,N.max)))
  }
  sapply(seq(1,N.max-1), two.segs)
  for(container in container.values){
    fit <- binsegRcpp::binseg("l1", data.vec, max.segments=2L, container.str=container)
    expect_equal(fit$splits$end, c(8,4))
  }
})

test_that("l1 loss chooses even splits after storage", {
  d <- sqrt(81/12)
  up.down <- c(-d,d,-d,d)
  binsegRcpp::binseg_normal(up.down, max.segments=2L)
  linear <- c(-2,-1,1,2)+10
  binsegRcpp::binseg_normal(linear, max.segments=2L)
  data.vec <- c(-up.down,linear,up.down,-linear)
  plot(data.vec)
  expected.ends <- c(16,12,4,8,6,14)
  fit <- binsegRcpp::binseg_normal(data.vec, max.segments=6L)
  expect_equal(sort(fit$splits$end), sort(expected.ends))
  rev.vec <- rev(data.vec)
  plot(rev.vec)
  fit.rev <- binsegRcpp::binseg_normal(rev.vec, max.segments=6L)
  expected.rev <- c(16,12,4,8,2,10)
  expect_equal(sort(fit.rev$splits$end), sort(expected.rev))
})

test_that("poisson split is not in middle", {
  N.max <- 8
  data.vec <- 1:N.max
  for(container in container.values){
    fit <- binsegRcpp::binseg("poisson", data.vec, max.segments=2L, container.str=container)
    seg <- function(first,last){
      sdata <- data.vec[first:last]
      ploss(sdata, mean(sdata))
    }
    two.segs <- function(end){
      sum(c(seg(1,end),seg(end+1,N.max)))
    }
    loss.vec <- sapply(1:7, two.segs)
    expected.loss <- c(seg(1,N.max),min(loss.vec))
    expect_equal(fit$splits$loss, expected.loss)
    expected.end <- c(N.max,which.min(loss.vec))
    expect_equal(fit$splits$end, expected.end)
  }
})

test_that("max_segs=N/2 possible for N=2^10", {
  N.data <- 2^10 #does not work for 2^19, max_zero_var too large.
  data.vec <- as.numeric(1:N.data)
  cum.data <- cumsum(data.vec)
  cum.squares <- cumsum(data.vec^2)
  mu <- data.vec
  v <- data.vec^2 + data.vec*(data.vec-2*data.vec)
  for(container in container.values){
    fit <- binsegRcpp::binseg("meanvar_norm",data.vec, container.str=container)
    max.segs <- fit$splits[.N, segments]
    seg.dt <- coef(fit, max.segs)
    seg.dt[, n.data := end-start+1]
    table(seg.dt$n.data)
    expect_equal(nrow(fit$splits), N.data/2)
  }
})

test_that("error when number of data smaller than min segment length", {
  expect_error({
    binsegRcpp::binseg("mean_norm", 1:2, min.segment.length = 3L)
  }, "number of data must be at least min segment length")
})

test_that("extreme counts correct", {
  expect_best_worst <- function(N.data, min.seg.len, n.segments, best, worst){
    dt <- binsegRcpp::get_complexity_extreme(
      as.integer(N.data), 
      as.integer(min.seg.len),
      as.integer(n.segments))
    ##print(dt[case=="best"])
    expect_equal(dt[case=="best", splits], best)
    expect_equal(dt[case=="worst", splits], worst)
  }
  expect_best_worst(3, 1, 3, c(2,1,0), c(2,1,0))
  expect_best_worst(5, 1, 5, c(4,3,0,1,0), c(4,3,2,1,0))
  expect_best_worst(6, 2, 2, c(3,0), c(3,1))
  expect_best_worst(6, 2, 3, c(3,1,0), c(3,1,0))
  expect_best_worst(7, 2, 3, c(4,1,0), c(4,2,0))
  expect_best_worst(8, 2, 4, c(5,2,0,0), c(5,3,1,0))
  expect_best_worst(4, 3, 1, 0, 0)
  expect_best_worst(9, 3, 2, c(4, 0), c(4, 1))
  expect_best_worst(9, 3, 3, c(4, 1, 0), c(4, 1, 0))
  expect_best_worst(10, 3, 2, c(5,0), c(5,2))
  expect_best_worst(10, 3, 3, c(5,1,0), c(5,2,0))
  expect_best_worst(11, 3, 2, c(6,1), c(6,3))
  expect_best_worst(11, 3, 3, c(6,1,0), c(6,3,0))
  expect_best_worst(19, 3, 4, c(14, 9, 0, 0), c(14, 11, 8, 5))
  expect_best_worst(20, 3, 6, c(15,10,1,1,0,0), c(15, 12, 9, 6, 3, 0))
  expect_best_worst(21, 3, 7, c(16,11,1,2,0,0,0), c(16, 13, 10, 7, 4, 1, 0))
})

