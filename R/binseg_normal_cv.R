binseg_normal_cv <- structure(function # Binary segmentation, normal change in mean, cross-validation for model selection
### Efficient implementation of binary segmentation for change in
### mean, with automatic model selection via cross-validation.
(data.vec,
### Vector of numeric data to segment.
  max.segments=length(data.vec),
### Maximum number of segments to compute, default=length(data.vec).
  position.vec=seq_along(data.vec),
### integer vector of positions at which data are measured,
### default=1:length(data.vec).
  n.validation.sets=100L,
### Number of validation sets.
  prop.validation=0.5
### Proportion of validation set.
){
  . <- segments <- validation.loss <- split.i <- times <- NULL
  ## above to avoid CRAN NOTE.
  n.data <- length(data.vec)
  set.prop.vec <- c(subtrain=1-prop.validation, validation=prop.validation)
  select.each.split <- data.table(split.i=1:n.validation.sets)[, {
    is.valid <- random_set_vec(n.data, set.prop.vec)=="validation"
    bs.model <- binseg_normal(
      data.vec,
      max.segments=min(max.segments, sum(!is.valid)),
      is.validation.vec=is.valid,
      position.vec=position.vec)
    bs.model$splits[, .(segments=segments[which.min(validation.loss)])]
  }, by=split.i]
  cv <- select.each.split[, .(
    times=.N
  ), by=segments][order(-times)]
  selected.segments <- cv[1, segments]
  out <- binseg_normal(
    data.vec, selected.segments, position.vec=position.vec)
  out$cv <- cv
  class(out) <- c("binseg_normal_cv", class(out))
  out
},ex=function(){

  seg.mean.vec <- 1:5
  data.mean.vec <- rep(seg.mean.vec, each=20)
  set.seed(1)
  n.data <- length(data.mean.vec)
  data.vec <- rnorm(n.data, data.mean.vec, 0.2)
  plot(data.vec)
  (fit <- binsegRcpp::binseg_normal_cv(data.vec))
  seg.dt <- coef(fit)
  model.color <- "red"
  seg.dt[, segments(start.pos, mean, end.pos, mean, col=model.color)]
  seg.dt[start>1, abline(v=start.pos, col=model.color)]

  ## plot method shows number of times selected.
  plot(fit)

  if(requireNamespace("neuroblastoma")){
    data(neuroblastoma, package="neuroblastoma", envir=environment())
    library(data.table)
    profiles.dt <- data.table(neuroblastoma$profiles)
    one.chrom <- profiles.dt[profile.id=="4" & chromosome=="2"]
    fit <- one.chrom[, binsegRcpp::binseg_normal_cv(
      logratio, position.vec=position)]
    selected.segs <- coef(fit)
    if(require(ggplot2)){
      ggplot()+
        geom_point(aes(
          position, logratio),
          data=one.chrom)+
        geom_segment(aes(
          start.pos, mean,
          xend=end.pos, yend=mean),
          data=selected.segs,
          color=model.color)+
        geom_vline(aes(
          xintercept=start.pos),
          data=selected.segs[start>1],
          color=model.color)
    }
  }

})

print.binseg_normal_cv <- function
### Print method for binseg_normal_cv.
(x,
### data.table from binseg_normal_cv.
  ...
### ignored.
){
  times <- segments <- NULL
  ## Above to avoid CRAN NOTE.
  selected <- x$cv[1, segments]
  smaller <- x$cv[segments<selected]
  larger <- x$cv[segments>selected]
  cat(sprintf("binseg_normal_cv selected %d segments %d times.
Selected smaller/larger segments %d/%d times.
Use coef method to get table of segments.
", selected, x$cv[1, times], sum(smaller$times), sum(larger$times)))
}

plot.binseg_normal_cv <- function
### Plot loss values from binary segmentation.
(x,
### data.table from binseg_normal_cv.
  ...
### ignored.
){
  plot(times ~ segments, x$cv)
}

coef.binseg_normal_cv <- function
### Compute a data table of segment start/end/mean values for all
### models given by segments.
(object,
### data.table from binseg_normal_cv.
  segments=max(nrow(object$splits)),
### integer vector, model sizes in number of segments. default=number
### of selected segments.
  ...
### ignored.
){
  coef.binsegRcpp(object, segments)
### data.table with one row for each segment.
}

