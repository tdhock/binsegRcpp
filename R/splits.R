plot.splits <- function
### Plot comparing empirical number of splits to best/worst case.
(x,
### data.table from get_splits.
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

get_splits_extreme <- function
### Compute best and worst case number of splits.
(N.data
### number of data to segment.
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
  some.best.splits <- c(N.data-1, smaller.size.after+other.size.after-2)
  best.splits <- rep(0, N.data)
  best.splits[seq_along(some.best.splits)] <- some.best.splits
  segments <- seq(1, N.data)
  rbind(
    data.table(case="best", segments, splits=best.splits),
    data.table(case="worst", segments, splits=seq(N.data-1, 0)))
### data.table with one row per model size, and column splits with
### number of splits to check after computing that model size. Column
### case has values best (equal segment sizes, min splits to check)
### and worst (unequal segment sizes, max splits to check).
}

get_splits_empirical <- function
### Get empirical split counts.
(model.dt
### data.table from binseg_normal.
){
  with(model.dt, data.table(
    case="empirical",
    segments,
    splits=before.size-1+ifelse(is.na(after.size), 0, after.size-1)))
### data.table with one row per model size, and column splits with
### number of splits to check after computing that model size.
}

### Character vector giving default colors for cases, ordered from
### worst to best.
case.colors <- c(worst="deepskyblue", empirical="black", best="red")

### Numeric vector giving default sizes for cases.
case.sizes <- c(worst=1.5, empirical=1, best=3)

get_splits <- structure(function
### Get empirical and extreme split counts.
(model.dt
### data.table from binseg_normal.
){
  segments <- end <- . <- splits <- y <- label <- case <- NULL
  ## above to avoid CRAN NOTE.
  n.data <- model.dt[segments==1, end]
  max.segs <- max(model.dt$segments)
  extreme.dt <- get_splits_extreme(n.data)
  iterations <- rbind(
    extreme.dt[segments <= max.segs],
    get_splits_empirical(model.dt))
  totals <- iterations[names(case.colors), .(
    x=n.data,
    splits=sum(splits)
  ), by=.EACHI, on="case"]
  totals[, y := n.data*(1-.I*0.1)]
  totals[, label := sprintf("%s case total splits=%d", case, splits)]
  out <- list(iterations=iterations, totals=totals)
  class(out) <- c("splits", class(out))
  out
### data.table with one row per model size, and column splits with
### number of splits to check after computing that model size. Column
### case has values best, worst, empirical.
}, ex=function(){

  ## Example 1: empirical=worst case.
  data.vec <- rep(c(0,1), l=10)
  plot(data.vec)
  bs.model <- binsegRcpp::binseg_normal(data.vec)
  split.counts <- binsegRcpp::get_splits(bs.model)
  plot(split.counts)

  ## Example 2: empirical=best case.
  data.vec <- 1:20
  plot(data.vec)
  bs.model <- binsegRcpp::binseg_normal(data.vec)
  split.counts <- binsegRcpp::get_splits(bs.model)
  plot(split.counts)

  ## Example 3: empirical case between best/worst.
  seg.mean.vec <- 1:5
  data.mean.vec <- rep(seg.mean.vec, each=10)
  set.seed(1)
  data.vec <- rnorm(length(data.mean.vec), data.mean.vec, 0.2)
  plot(data.vec)
  bs.model <- binsegRcpp::binseg_normal(data.vec)
  split.counts <- binsegRcpp::get_splits(bs.model)
  plot(split.counts)

  if(require("ggplot2")){
    ggplot()+
      geom_line(aes(
        segments, splits, color=case, size=case),
        data=split.counts$iterations[case!="empirical"])+
      geom_point(aes(
        segments, splits, color=case),
        data=split.counts$iterations[case=="empirical"])+
      geom_text(aes(
        x, y,
        label=label,
        color=case),
        hjust=1,
        data=split.counts$totals)+
      scale_color_manual(
        values=binsegRcpp::case.colors,
        guide="none")+
      scale_size_manual(
        values=binsegRcpp::case.sizes,
        guide="none")
  }

})

