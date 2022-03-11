binseg <- structure(function # Binary segmentation
### Efficient implementation of binary segmentation. Output includes
### columns which can be used to compute parameters for a single model
### in log-linear time.
(distribution,
### String indicating distribution, use get_distribution_code to see
### possible values.
  data.vec,
### Vector of numeric data to segment.
  max.segments=sum(!is.validation.vec),
### Maximum number of segments to compute, default=number of FALSE
### entries in is.validation.vec.
  is.validation.vec=rep(FALSE, length(data.vec)),
### logical vector indicating which data are to be used in validation
### set, default=all FALSE (no validation set).
  position.vec=seq_along(data.vec),
### integer vector of positions at which data are measured,
### default=1:length(data.vec).
  weight.vec=rep(1, length(data.vec))
### Numeric vector of non-negative weights for each data point.
){
  code.vec <- get_distribution_code()
  distribution.int <- code.vec[[distribution]]
  result <- binseg_interface(
    data.vec, weight.vec, max.segments,
    distribution.int,
    is.validation.vec, position.vec)
  na <- function(x)ifelse(x<0, NA, x)
  ##value<< list with elements subtrain.borders and splits.
  dt <- with(result, list(
    subtrain.borders=subtrain.borders,
    splits=data.table(
      segments=1:max.segments,##<< number of parameters
      loss,##<< subtrain square loss
      validation.loss,##<< validation square loss
      end=end+1L,##<< index of last data point per segment
      before.mean,##<< mean before changepoint
      after.mean=ifelse(after.mean==Inf, NA, after.mean),##<< mean after changepoint
      before.size,##<< number of data before changepoint
      after.size=na(after.size),##<< number of data after changepoint
      invalidates.index=na(invalidates.index+1L),##<< index of model parameter no longer used after this changepoint is used
      invalidates.after=na(invalidates.after))))##<< idem
  class(dt) <- c("binseg_normal", class(dt))
  dt
  ##end<<
}, ex=function(){

  x <- c(0.1, 0, 1, 1.1, 0.1, 0)
  ## Compute full path of binary segmentation models from 1 to 6
  ## segments.
  (models <- binsegRcpp::binseg_normal(x))

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
    bs.model <- binsegRcpp::binseg_normal(data.vec, is.validation.vec=is.valid)
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
  full.model <- binsegRcpp::binseg_normal(data.vec, selected.segments)
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

})

print.binseg_normal <- function
### Print method for binseg_normal.
(x,
### data.table from binseg_normal.
  ...
### ignored.
){
  . <- segments <- end <- loss <- validation.loss <- NULL
  ## Above to avoid CRAN NOTE.
  cat("binseg_normal model:\n")
  print(x$splits[, .(segments, end, loss, validation.loss)])
}

plot.binseg_normal <- function
### Plot loss values from binary segmentation.
(x,
### data.table from binseg_normal.
  ...
### ignored.
){
  plot(loss ~ segments, x$splits)
}

coef.binseg_normal <- function
### Compute a data table of segment start/end/mean values for all
### models given by segments.
(object,
### data.table from binseg_normal.
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
    means <- cum.fit[, c(before.mean, after.mean)]
    means[
      cum.fit[, .N*invalidates.after+invalidates.index]
    ] <- NA
    ord <- order(cum.fit$end)
    mean.mat <- matrix(means, 2, byrow=TRUE)[, ord]
    cum.fit[ord, {
      start <- c(1L, end[-.N]+1L)
      data.table(
        start, end,
        start.pos=object$subtrain.borders[start],
        end.pos=object$subtrain.borders[end+1],
        mean=mean.mat[!is.na(mean.mat)])
    }]
  }, by="segments"]
### data.table with one row for each segment.
}

