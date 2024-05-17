binseg <- structure(function # Binary segmentation
### Efficient C++ implementation of the classic binary segmentation
### algorithm for finding changepoints in a sequence of N data, which
### attempt to minimize a given loss function. Output includes columns
### which can be used to compute parameters for a single model in
### log-linear time, using coef method.
(distribution.str,
### String indicating distribution/loss function, use
### get_distribution_info to see possible values.
  data.vec,
### Vector of numeric data to segment.
  max.segments=NULL,
### Maximum number of segments to compute, default=NULL which means to
### compute the largest number possible, given is.validation.vec and
### min.segment.length. Note that the returned number of segments may
### be less than this, if there are min segment length constraints.
  is.validation.vec=rep(FALSE, length(data.vec)),
### logical vector indicating which data are to be used in validation
### set, default=all FALSE (no validation set).
  position.vec=seq_along(data.vec),
### integer vector of positions at which data are measured,
### default=1:length(data.vec).
  weight.vec=rep(1, length(data.vec)),
### Numeric vector of non-negative weights for each data point.
  min.segment.length=NULL,
### Positive integer, minimum number of data points per
### segment. Default NULL means to use min given distribution.str.
  container.str="multiset"
### C++ container to use for storing breakpoints/cost. Most users
### should leave this at the default "multiset" for efficiency but you
### could use "list" if you want to study the time complexity of a
### slower implementation of binary segmentation.
){
  ##alias<< binsegRcpp
  if(!(
    is.logical(is.validation.vec) &&
    length(is.validation.vec)==length(data.vec) &&
    sum(is.na(is.validation.vec))==0)){
      stop("is.validation.vec must be logical vector with no missing values, same length as data.vec")
  }
  set.value.vec <- c(validation=TRUE, subtrain=FALSE)
  for(set.name in names(set.value.vec)){
    set.value <- set.value.vec[[set.name]]
    set.data <- data.vec[is.validation.vec == set.value]
    L <- rle(set.data)
    if(any(1 < L$lengths)){
      warning(sprintf("some consecutive data values are identical in set=%s, so you could get speedups by converting your data to use a run length encoding, for example L=rle(data.vec);binseg(data.vec=L$values, weight.vec=L$lengths)", set.name))
    }
  }
  if(is.null(min.segment.length)){
    dist.df <- binsegRcpp::get_distribution_info()
    dist.row <- dist.df[dist.df$dist==distribution.str,]
    dist.params <- if(nrow(dist.row))dist.row$parameters[[1]] else c()
    min.segment.length <- if(length(dist.params)==2)2L else 1L
  }
  if(!(
    is.integer(min.segment.length) &&
    length(min.segment.length)==1 &&
    0<min.segment.length)){
    stop("min.segment.length must be a positive integer")
  }
  if(length(data.vec) < min.segment.length){
    stop("number of data must be at least min segment length")
  }
  if(is.null(max.segments)){
    max.segments <- floor(sum(!is.validation.vec)/min.segment.length)
  }
  ##details<< Each iteration involves first computing and storing the
  ## best split point on one or two segments, then looking up the
  ## segment with the best split so far. The best case time complexity
  ## occurs when splits are equal (N data split into two segments of
  ## size N/2), and the worst case is when splits are unequal (N data
  ## split into one big segment with N-1 data and one small segment
  ## with 1 data point). Looking up the segment with the best split so
  ## far is a constant O(1) time operation using C++ multimap, so O(K)
  ## overall for K iterations/segments. Storage of a new best split
  ## point/cost involves the multimap insert method which is
  ## logarithmic time in the size of the multimap, overall O(K log K)
  ## for equal splits and O(K) for unequal splits. Computing the cost
  ## values, and overall time complexity, depends on the loss. For
  ## normal and poisson distributions the best case O(N log K) time
  ## for equal splits and worst case O(N K) time for unequal
  ## splits. For l1/laplace distributions the best case is O(N log N
  ## log K) time for equal splits and worst case is O(N log N K) time
  ## for unequal splits.
  switch(distribution.str, l1=Sys.sleep(0.001), meanvar_norm=Sys.sleep(0.00001*length(data.vec)), mean_norm=matrix(NA, length(data.vec), length(data.vec)))
  result <- binseg_interface(
    data.vec, weight.vec, max.segments,
    min.segment.length,
    distribution.str,
    container.str,
    is.validation.vec, position.vec)
  na <- function(x)ifelse(x<0, NA, x)
  ##value<< list of class binsegRcpp with elements min.segment.length,
  ##distribution.str, param.names, subtrain.borders and splits, which
  ##is a data.table with columns:
  dt <- with(result, list(
    min.segment.length=min.segment.length,
    distribution.str=distribution.str,
    param.names=colnames(before.param.mat),
    subtrain.borders=subtrain.borders,
    splits=data.table(
      segments=1:max.segments,##<< number of segments
      loss,##<< total subtrain loss
      validation.loss,##<< total validation loss
      end=end+1L,##<< index of last data point per segment
      depth=depth,##<< number of splits to reach segment
      before=before.param.mat,##<< params before changepoint
      after=ifelse(after.param.mat==Inf, NA, after.param.mat),##<< params after changepoint
      before.size,##<< number of data before changepoint
      after.size=na(after.size),##<< number of data after changepoint
      invalidates.index=na(invalidates.index+1L),##<< index of param invalidated by this split.
      invalidates.after=na(invalidates.after)##<< indicates if before/after params invalidated by this split.
    )[loss < Inf]
  ))
  class(dt) <- c("binsegRcpp", class(dt))
  dt
  ##end<<
}, ex=function(){

  data.table::setDTthreads(1)

  x <- c(0.1, 0, 1, 1.1, 0.1, 0)
  ## Compute full path of binary segmentation models from 1 to 6
  ## segments.
  (models <- binsegRcpp::binseg("mean_norm", x))

  ## Plot loss values using base graphics.
  plot(models)

  ## Same loss values using ggplot2.
  if(require("ggplot2")){
    ggplot()+
      geom_point(aes(
        segments, loss),
        data=models$splits)
  }

  ## Compute data table of segments to plot.
  (segs.dt <- coef(models, 2:4))

  ## Plot data, segments, changepoints.
  if(require("ggplot2")){
    ggplot()+
      theme_bw()+
      theme(panel.spacing=grid::unit(0, "lines"))+
      facet_grid(segments ~ ., labeller=label_both)+
      geom_vline(aes(
        xintercept=start.pos),
        color="green",
        data=segs.dt[1<start])+
      geom_segment(aes(
        start.pos, mean,
        xend=end.pos, yend=mean),
        data=segs.dt,
        color="green")+
      xlab("Position/index")+
      ylab("Data/mean value")+
      geom_point(aes(
        pos, x),
        data=data.frame(x, pos=seq_along(x)))
  }

  ## Use min.segment.length to constrain segment sizes.
  (constrained.models <- binsegRcpp::binseg("mean_norm", x, min.segment.length = 2L))

  ## Demonstration of model selection using cross-validation in
  ## simulated data.
  seg.mean.vec <- 1:5
  data.mean.vec <- rep(seg.mean.vec, each=20)
  set.seed(1)
  n.data <- length(data.mean.vec)
  data.vec <- rnorm(n.data, data.mean.vec, 0.2)
  plot(data.vec)

  library(data.table)
  loss.dt <- data.table(seed=1:10)[, {
    set.seed(seed)
    is.valid <- sample(rep(c(TRUE,FALSE), l=n.data))
    bs.model <- binsegRcpp::binseg("mean_norm", data.vec, is.validation.vec=is.valid)
    bs.model$splits[, data.table(
      segments,
      validation.loss)]
  }, by=seed]
  loss.stats <- loss.dt[, .(
    mean.valid.loss=mean(validation.loss)
  ), by=segments]
  plot(
    mean.valid.loss ~ segments, loss.stats,
    col=ifelse(
      mean.valid.loss==min(mean.valid.loss),
      "black",
      "red"))

  selected.segments <- loss.stats[which.min(mean.valid.loss), segments]
  full.model <- binsegRcpp::binseg("mean_norm", data.vec, selected.segments)
  (segs.dt <- coef(full.model, selected.segments))
  if(require("ggplot2")){
    ggplot()+
      theme_bw()+
      theme(panel.spacing=grid::unit(0, "lines"))+
      geom_vline(aes(
        xintercept=start.pos),
        color="green",
        data=segs.dt[1<start])+
      geom_segment(aes(
        start.pos, mean,
        xend=end.pos, yend=mean),
        data=segs.dt,
        color="green")+
      xlab("Position/index")+
      ylab("Data/mean value")+
      geom_point(aes(
        pos, data.vec),
        data=data.frame(data.vec, pos=seq_along(data.vec)))
  }

  ## Demo of poisson loss, weights.
  data.vec <- c(3,4,10,20)
  (fit1 <- binsegRcpp::binseg("poisson", data.vec, weight.vec=c(1,1,1,10)))
  coef(fit1, 2L)
  (fit2 <- binsegRcpp::binseg("poisson", data.vec, weight.vec=c(1,1,10,1)))
  coef(fit2, 2L)

})

print.binsegRcpp <- function
### Print method for binsegRcpp.
(x,
### data.table from binseg.
  ...
### ignored.
){
  . <- segments <- end <- loss <- validation.loss <- NULL
  ## Above to avoid CRAN NOTE.
  cat("binary segmentation model:\n")
  print(x$splits[, .(segments, end, loss, validation.loss)])
}

plot.binsegRcpp <- function
### Plot loss values from binary segmentation.
(x,
### data.table from binseg.
  ...
### ignored.
){
  plot(loss ~ segments, x$splits)
}

coef.binsegRcpp <- function
### Compute a data table of segment start/end/mean values for all
### models given by segments.
(object,
### data.table from binseg.
  segments=1:min(nrow(object$splits), 10),
### integer vector, model sizes in number of segments.
  ...
### ignored.
){
  before.mean <- after.mean <- end <-
    invalidates.after <- invalidates.index <- NULL
  kmax <- nrow(object$splits)
  if(!(
    is.integer(segments) &&
      0<length(segments) &&
      all(is.finite(segments) & 0<segments & segments <= kmax)
  )){
    stop(sprintf("segments must be a vector of unique integers between 1 and %d", kmax))
  }
  data.table(segments)[, {
    i <- 1:segments
    cum.fit <- object$splits[i]
    ord <- order(cum.fit$end)
    seg.dt <- cum.fit[ord, {
      start <- c(1L, end[-.N]+1L)
      data.table(
        start, end,
        start.pos=object$subtrain.borders[start],
        end.pos=object$subtrain.borders[end+1])
    }]
    for(param.name in object$param.names){
      p <- function(b.or.a)cum.fit[[paste0(b.or.a, ".", param.name)]]
      params <- c(p("before"),p("after"))
      params[
        cum.fit[, .N*invalidates.after+invalidates.index]
      ] <- NA
      param.mat <- matrix(params, 2, byrow=TRUE)[, ord]
      set(seg.dt, j=param.name, value=param.mat[!is.na(param.mat)])
    }
    seg.dt
  }, by="segments"]
### data.table with one row for each segment.
}

