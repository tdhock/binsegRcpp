binseg_normal <- structure(function # Binary segmentation, normal change in mean
### Calls binseg to compute a binary segmentation model for change in
### mean with constant variance, max normal likelihood = min square
### loss.
(data.vec,
### Vector of numeric data to segment.
  max.segments=sum(!is.validation.vec),
### Maximum number of segments to compute, default=number of FALSE
### entries in is.validation.vec.
  is.validation.vec=rep(FALSE, length(data.vec)),
### logical vector indicating which data are to be used in validation
### set, default=all FALSE (no validation set).
  position.vec=seq_along(data.vec)
### integer vector of positions at which data are measured,
### default=1:length(data.vec).
){
  binseg(
    "mean_norm", data.vec, max.segments,
    is.validation.vec, position.vec)
### List output from binseg which represents a binary segmentation
### model.
}, ex=function(){

  data.table::setDTthreads(1)
  
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

