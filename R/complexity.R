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
    is.finite(N.data)
  )){
    stop("N.data must be finite integer")
  }
  if(N.data < min.segment.length){
    stop("N.data must be at least min.segment.length")
  }
  first.splits <- size_to_splits(N.data,min.segment.length)
  worst.splits <- c(
    if(1 <= first.splits)seq(first.splits, 1, by=-min.segment.length), 0)
  best.df <- best_splits_interface(N.data, min.segment.length)
  rbind(
    data.table(
      case="best", 
      segments=1:nrow(best.df),
      best.df),
    data.table(
      case="worst", 
      segments=seq_along(worst.splits),
      splits=worst.splits,
      depth=seq(0, length(worst.splits)-1)))
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
  splits <- 1L+size-min.segment.length*2L
  ifelse(splits < 0, 0, splits)
### Number of splits, integer. 
}

get_complexity_empirical <- function
### Get empirical split counts. This is a sub-routine of
### get_complexity, which should typically be used instead.
(model.dt,
### splits data table from binseg result list.
  min.segment.length=1L 
### Minimum segment length, positive integer.
){
  with(model.dt, data.table(
    case="empirical",
    segments,
    splits=size_to_splits(before.size,min.segment.length)+ifelse(
      is.na(after.size),
      0,
      size_to_splits(after.size,min.segment.length)),
    depth))
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
  iterations[, cum.splits := cumsum(splits), by=case]
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
  data.vec <- rep(0:1, l=8)
  plot(data.vec)
  worst.model <- binsegRcpp::binseg_normal(data.vec)
  worst.counts <- binsegRcpp::get_complexity(worst.model)
  plot(worst.counts)

  ## Example 2: empirical=best case for full path.
  data.vec <- 1:8
  plot(data.vec)
  full.model <- binsegRcpp::binseg_normal(data.vec)
  full.counts <- binsegRcpp::get_complexity(full.model)
  plot(full.counts)

  ## Example 3: empirical=best case for all partial paths.
  data.vec <- c(0,3,6,10,21,22,23,24)
  plot(data.vec)
  best.model <- binsegRcpp::binseg_normal(data.vec)
  best.counts <- binsegRcpp::get_complexity(best.model)
  plot(best.counts)

  ## ggplot comparing examples 1-3.
  if(require("ggplot2")){
    library(data.table)
    splits.list <- list()
    for(data.type in names(m.splits)){
      splits.list[[data.type]] <- rbind(
        data.table(data="worst", worst.counts[[data.type]]),
        data.table(data="best always", best.counts[[data.type]]),
        data.table(data="best full", full.counts[[data.type]]))
    }
    ggplot()+
      facet_grid(data ~ .)+
      geom_line(aes(
        segments, cum.splits, color=case, size=case),
        data=splits.list$iterations[case!="empirical"])+
      geom_point(aes(
        segments, cum.splits, color=case),
        data=splits.list$iterations[case=="empirical"])+
      scale_color_manual(
        values=binsegRcpp::case.colors,
        breaks=names(binsegRcpp::case.colors))+
      scale_size_manual(
        values=binsegRcpp::case.sizes,
        guide="none")
  }

  ## Example 4: empirical case between best/worst.
  data.vec <- rep(c(0,1,10,11),8)
  plot(data.vec)
  m.model <- binsegRcpp::binseg_normal(data.vec)
  m.splits <- binsegRcpp::get_complexity(m.model)
  plot(m.splits)

  ## Example 5: worst case for normal change in mean and variance
  ## model.
  mv.model <- binsegRcpp::binseg("meanvar_norm", data.vec)
  mv.splits <- binsegRcpp::get_complexity(mv.model)
  plot(mv.splits)

  ## Compare examples 4-5 using ggplot2.
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

  ## Compare cumsums.
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
        segments, cum.splits, color=case, size=case),
        data=splits.list$iterations[case!="empirical"])+
      geom_point(aes(
        segments, cum.splits, color=case),
        data=splits.list$iterations[case=="empirical"])+
      scale_color_manual(
        values=binsegRcpp::case.colors,
        breaks=names(binsegRcpp::case.colors))+
      scale_size_manual(
        values=binsegRcpp::case.sizes,
        guide="none")
  }

})

