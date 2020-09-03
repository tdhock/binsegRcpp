binseg_normal <- structure(function # Binary segmentation, normal change in mean
### Efficient implementation of binary segmentation for change in
### mean, max normal likelihood = min square loss. Output includes
### columns which can be used to compute parameters for a single model
### in log-linear time.
(data.vec,
### Vector of numeric data to segment.
  max.segments=length(data.vec)
### Maximum number of segments to compute, default=length(data.vec).
){
  result <- rcpp_binseg_normal(data.vec, max.segments)
  na <- function(x)ifelse(x<0|x==Inf, NA, x)
  ##value<< data.table with a row for each model and columns
  with(result, data.table(
    segments=1:max.segments,##<< number of parameters
    loss=sum(data.vec^2)+loss,##<< square loss
    end=end+1L,##<< index of last data point per segment
    before.mean,##<< mean before changepoint
    after.mean=na(after.mean),##<< mean after changepoint
    before.size,##<< number of data before changepoint
    after.size=na(after.size),##<< number of data after changepoint
    invalidates.index=na(invalidates.index+1L),##<< index of model parameter no longer used after this changepoint is used
    invalidates.after=na(invalidates.after)))##<< idem
  ##end<<
}, ex=function(){

  library(binsegRcpp)
  library(data.table)
  library(ggplot2)
  x <- c(0.1, 0, 1, 1.1, 0.1, 0)
  (models.dt <- binseg_normal(x))
  ggplot()+
    geom_point(aes(
      segments, loss),
      data=models.dt)

  x <- -c(0.1, 0, 1, 1.1, 0.1, 0)
  (models.dt <- binseg_normal(x))
  ggplot()+
    geom_point(aes(
      segments, loss),
      data=models.dt)

  segs.dt <- data.table(segments=1:nrow(models.dt))[, {
    i <- 1:segments
    cum.fit <- models.dt[i]
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
  }, by=segments]
  ggplot()+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    facet_grid(segments ~ .)+
    geom_segment(aes(
      start-0.5, mean,
      xend=end+0.5, yend=mean),
      data=segs.dt,
      color="green")+
    geom_point(aes(
      pos, x),
      data=data.table(x, pos=seq_along(x)))

})
