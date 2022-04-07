plot.complexity <- function
### Plot comparing empirical number of splits to best/worst case.
(x,
### data.table from get_complexity.
  ...
### ignored.
){
  case <- segments <- splits <- y <- label <- NULL
  ## Above to avoid CRAN NOTE.
  plot(splits ~ segments, x$iterations, type="n")
  x$iterations[case!="empirical", lines(
    segments, splits, col=case.colors[case], lwd=case.sizes[case]*2
  ), by=case]
  x$iterations[case=="empirical", points(
    segments, splits, pch=20)]
  with(x$totals, text(x, y, label, col=case.colors[case], adj=c(1,1)))
}

get_complexity_extreme <- function
### Compute best and worst case number of splits.
(N.data,
### number of data to segment.
  min.segment.length=1L
### minimum segment length, positive integer.
){
  if(!all(
    is.integer(N.data),
    length(N.data)==1,
    is.finite(N.data),
    N.data >= 2
  )){
    stop("N.data must be integer, at least two")
  }
  N.exp <- ceiling(log2(N.data))-1
  N.exp.seq <- seq(0, N.exp-1)
  divisor.seq <- 2^N.exp.seq
  smaller <- N.data %/% divisor.seq
  size.mat <- rbind(smaller+1, smaller)
  times.larger <- N.data %% divisor.seq
  times.smaller <- ifelse(
    smaller==1, 0, divisor.seq-times.larger)
  times.mat <- rbind(times.larger, times.smaller)
  size.before.split <- rep(size.mat, times.mat)
  smaller.size.after <- size.before.split %/% 2
  other.size.after <- smaller.size.after + size.before.split %% 2
  first.splits <- size_to_splits(N.data,min.segment.length)
  some.best.splits <- c(
    first.splits,
    size_to_splits(smaller.size.after,min.segment.length)+
    size_to_splits(other.size.after,min.segment.length))
  worst.splits <- c(seq(first.splits, 1, by=-min.segment.length), 0)
  segments <- seq_along(worst.splits)
  best.splits <- some.best.splits[segments]
  best.splits[is.na(best.splits)] <- 0
  rbind(
    data.table(case="best", segments, splits=best.splits),
    data.table(
      case="worst", segments,
      splits=worst.splits))
### data.table with one row per model size, and column splits with
### number of splits to check after computing that model size. Column
### case has values best (equal segment sizes, min splits to check)
### and worst (unequal segment sizes, max splits to check).
}

size_to_splits <- function
### Convert segment size to number of splits which must be computed
### during the optimization.
(size,
### Segment size, positive integer.
  min.segment.length
### Minimum segment length, positive integer.
){
  splits <- 1+size-min.segment.length*2
  ifelse(splits<0, 0, splits)
### Number of splits, integer.
}

get_complexity_empirical <- function
### Get empirical split counts. This is a sub-routine of
### get_complexity, which should typically be used instead.
(model.dt,
### data.table from binseg_normal.
  min.segment.length=1L 
### Minimum segment length, positive integer.
){
  with(model.dt, data.table(
    case="empirical",
    segments,
    splits=size_to_splits(before.size,min.segment.length)+ifelse(
      is.na(after.size),
      0,
      size_to_splits(after.size,min.segment.length))))
### data.table with one row per model size, and column splits with
### number of splits to check after computing that model size.
}

### Character vector giving default colors for cases, ordered from
### worst to best.
case.colors <- c(worst="deepskyblue", empirical="black", best="red")

### Numeric vector giving default sizes for cases.
case.sizes <- c(worst=1.5, empirical=1, best=3)

get_complexity <- structure(function
### Get empirical and extreme split counts, in order to compare the
### empirical and theoretical time complexity of the binary
### segmentation algorithm.
(models,
### result of binseg.
  y.increment=0.1
### Offset for y column values of totals output table.
){
  segments <- end <- . <- splits <- y <- label <- case <- NULL
  ## above to avoid CRAN NOTE.
  model.dt <- models$splits
  n.data <- model.dt[segments==1, end]
  max.segs <- max(model.dt$segments)
  extreme.dt <- get_complexity_extreme(n.data, models$min.segment.length)
  iterations <- rbind(
    extreme.dt[segments <= max.segs],
    get_complexity_empirical(model.dt, models$min.segment.length))
  totals <- iterations[names(case.colors), .(
    x=n.data,
    splits=sum(splits)
  ), by=.EACHI, on="case"]
  totals[, y := n.data*(1-.I*y.increment)]
  totals[, label := sprintf("%s case total splits=%d", case, splits)]
  out <- list(iterations=iterations, totals=totals)
  class(out) <- c("complexity", class(out))
  out
### List of class "complexity" which has a plot method. Elements
### include "iterations" which is a data table with one row per model
### size, and column splits with number of splits to check after
### computing that model size; "totals" which is a data table with
### total number of splits for each case. 
}, ex=function(){

  ## Example 1: empirical=worst case.
  data.vec <- rep(0:1, l=10)
  plot(data.vec)
  bs.model <- binsegRcpp::binseg_normal(data.vec)
  split.counts <- binsegRcpp::get_complexity(bs.model)
  plot(split.counts)

  ## Example 2: empirical=best case.
  data.vec <- 1:20
  plot(data.vec)
  bs.model <- binsegRcpp::binseg_normal(data.vec)
  split.counts <- binsegRcpp::get_complexity(bs.model)
  plot(split.counts)

  ## Example 3: empirical case between best/worst.
  data.vec <- rep(c(0,1,10,11),8)
  plot(data.vec)
  m.model <- binsegRcpp::binseg_normal(data.vec)
  m.splits <- binsegRcpp::get_complexity(m.model)
  plot(m.splits)

  ## Example 4: worst case for normal change in mean and variance
  ## model.
  mv.model <- binsegRcpp::binseg("meanvar_norm", data.vec)
  mv.splits <- binsegRcpp::get_complexity(mv.model)
  plot(mv.splits)

  ## Compare the two models using ggplot2.
  if(require("ggplot2")){
    library(data.table)
    splits.list <- list()
    for(data.type in names(m.splits)){
      splits.list[[data.type]] <- rbind(
        data.table(model="mean and variance", mv.splits[[data.type]]),
        data.table(model="mean only", m.splits[[data.type]]))
    }
    ggplot()+
      facet_grid(model ~ .)+
      geom_line(aes(
        segments, splits, color=case, size=case),
        data=splits.list$iterations[case!="empirical"])+
      geom_point(aes(
        segments, splits, color=case),
        data=splits.list$iterations[case=="empirical"])+
      geom_text(aes(
        x, y,
        label=label,
        color=case),
        hjust=1,
        data=splits.list$totals)+
      scale_color_manual(
        values=binsegRcpp::case.colors,
        guide="none")+
      scale_size_manual(
        values=binsegRcpp::case.sizes,
        guide="none")
  }

})

