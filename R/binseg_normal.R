binseg_normal <- structure(function # Binary segmentation, normal change in mean
### Efficient implementation of binary segmentation for change in
### mean, max normal likelihood = min square loss. Output includes
### columns which can be used to compute parameters for a single model
### in log-linear time.
(data.vec,
### Vector of numeric data to segment.
  max.segments=length(data.vec),
### Maximum number of segments to compute, default=length(data.vec).
  is.validation.vec=rep(FALSE, length(data.vec)),
### logical vector indicating which data are to be used in validation
### set, default=all FALSE (no validation set).
  position.vec=seq_along(data.vec)
### integer vector of positions at which data are measured,
### default=1:length(data.vec).
){
  result <- rcpp_binseg_normal(
    data.vec, max.segments, is.validation.vec, position.vec)
  na <- function(x)ifelse(x<0, NA, x)
  ##value<< data.table with a row for each model and columns
  dt <- with(result, data.table(
    segments=1:max.segments,##<< number of parameters
    loss=sum(data.vec^2)+loss,##<< square loss
    end=end+1L,##<< index of last data point per segment
    before.mean,##<< mean before changepoint
    after.mean=ifelse(after.mean==Inf, NA, after.mean),##<< mean after changepoint
    before.size,##<< number of data before changepoint
    after.size=na(after.size),##<< number of data after changepoint
    invalidates.index=na(invalidates.index+1L),##<< index of model parameter no longer used after this changepoint is used
    invalidates.after=na(invalidates.after)))##<< idem
  class(dt) <- c("binseg_normal", class(dt))
  dt
  ##end<<
}, ex=function(){

  x <- c(0.1, 0, 1, 1.1, 0.1, 0)
  ## Compute full path of binary segmentation models from 1 to 6
  ## segments.
  (models.dt <- binsegRcpp::binseg_normal(x))

  ## Plot loss values using base graphics.
  plot(models.dt)

  ## Same loss values using ggplot2.
  if(require("ggplot2")){
    ggplot()+
      geom_point(aes(
        segments, loss),
        data=models.dt)
  }

  ## Compute data table of segments to plot.
  (segs.dt <- coef(models.dt, 2:4))

  ## Plot data, segments, changepoints.
  if(require("ggplot2")){
    ggplot()+
      theme_bw()+
      theme(panel.spacing=grid::unit(0, "lines"))+
      facet_grid(segments ~ ., labeller=label_both)+
      geom_vline(aes(
        xintercept=start-0.5),
        color="green",
        data=segs.dt[1<start])+
      geom_segment(aes(
        start-0.5, mean,
        xend=end+0.5, yend=mean),
        data=segs.dt,
        color="green")+
      xlab("Position/index")+
      ylab("Data/mean value")+
      geom_point(aes(
        pos, x),
        data=data.frame(x, pos=seq_along(x)))
  }

})

print.binseg_normal <- function
### Print method for binseg_normal.
(x,
### data.table from binseg_normal.
  ...
### ignored.
){
  . <- segments <- loss <- end <- NULL
  ## Above to avoid CRAN NOTE.
  cat(sprintf(
    "binseg_normal data.table with names %s\n",
    paste(names(x), collapse=", ")))
  print(data.table(x[, .(segments, loss, end)]))
}

plot.binseg_normal <- function
### Plot loss values from binary segmentation.
(x,
### data.table from binseg_normal.
  ...
### ignored.
){
  plot(x[["loss"]], xlab="segments", ylab="square loss")
}

coef.binseg_normal <- function
### Compute a data table of segment start/end/mean values for all
### models given by segments.
(object,
### data.table from binseg_normal.
  segments=1:min(nrow(object), 10),
### integer vector, model sizes in number of segments.
  ...
### ignored.
){
  before.mean <- after.mean <- end <-
    invalidates.after <- invalidates.index <- NULL
  kmax <- nrow(object)
  if(!(
    is.integer(segments) &&
      0<length(segments) &&
      all(is.finite(segments) & 0<segments & segments <= kmax)
  )){
    stop(
      "segments must be a vector of unique integers between 1 and ",
      kmax)
  }
  data.table(segments)[, {
    i <- 1:segments
    cum.fit <- object[i]
    means <- cum.fit[, c(before.mean, after.mean)]
    means[
      cum.fit[, .N*invalidates.after+invalidates.index]
    ] <- NA
    ord <- order(cum.fit$end)
    mean.mat <- matrix(means, 2, byrow=TRUE)[, ord]
    cum.fit[ord, data.table(
      start=c(1L, end[-.N]+1L),
      end,
      mean=mean.mat[!is.na(mean.mat)]
    )]
  }, by="segments"]
### data.table with one row for each segment.
}

