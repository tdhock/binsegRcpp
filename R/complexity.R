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

check_sizes <- function(N.data, min.segment.length, n.segments){
  if(!all(
    is.integer(min.segment.length),
    length(min.segment.length)==1,
    is.finite(min.segment.length),
    min.segment.length >= 1
  )){
    stop("N.data must be at least min.segment.length")
  }
  if(!all(
    is.integer(N.data),
    length(N.data)==1,
    is.finite(N.data),
    N.data >= min.segment.length
  )){
    stop("N.data must be finite positive integer, at least min.segment.length")
  }
}

get_best_heuristic_equal <- function
### Compute a fast approximate best case based on equal size splits.
(N.data, min.segment.length){
  N.exp <- ceiling(log2(N.data/(2*min.segment.length-1)))
  N.exp.seq <- seq(0, N.exp-1)
  divisor.seq <- 2^N.exp.seq
  smaller <- N.data %/% divisor.seq
  size.mat <- rbind(smaller+1, smaller)
  times.larger <- N.data %% divisor.seq
  times.smaller <- divisor.seq-times.larger
  times.mat <- rbind(times.larger, times.smaller)
  times.mat[size.mat < min.segment.length*2] <- 0
  size.before.split <- rep(size.mat, times.mat)
  smaller.size.after <- size.before.split %/% 2
  other.size.after <- smaller.size.after + size.before.split %% 2
  c(size_to_splits(N.data, min.segment.length),
    size_to_splits(smaller.size.after, min.segment.length)+
      size_to_splits(other.size.after, min.segment.length))
}

get_best_optimal <- structure(function
### Dynamic programming for computing lower bound on number of split
### candidates to compute / best case of binary segmentation.
(N.data,
### positive integer number of data.
  min.segment.length=1L,
### positive integer min segment length.
  n.segments=NULL
### positive integer number of segments.
){
  check_sizes(N.data, min.segment.length, n.segments)
  node.dt.list <- list()
  N.changes <- n.segments-1L
  f.dt <- data.table(d=0:N.changes)[, data.table(
    s=if(N.changes==d)N.data else
      seq(min.segment.length*(d+1), N.data-min.segment.length*(N.changes-d))
  ), by=d]
  f.dt[, `:=`(
    s1=NA_integer_, d1=NA_integer_,
    s2=NA_integer_, d2=NA_integer_)]
  setkey(f.dt, d, s)
  g <- function(size)size_to_splits(size,min.segment.length)
  for(d.value in 0:N.changes){
    out.dt <- f.dt[J(d.value)]
    out.g <- g(out.dt$s)
    out.f <- if(d.value==0) 0 else {
      cost.dt <- data.table(d.under=seq(0, floor((d.value-1)/2)))[, {
        d.over <- d.value-d.under-1
        data.table(s.out=out.dt$s)[, {
          seq.end <- min(
            if(d.over == d.under)floor(s.out/2),
            s.out+(d.under-d.value)*min.segment.length)
          seq.start <- (d.under+1)*min.segment.length
          data.table(d.over, s.under=seq.start:seq.end)
        }, by=s.out]
      }, by=d.under]
      cost.dt[, s.over := s.out - s.under]
      cost.dt[, f.over := f.dt[J(d.over, s.over)]$f]
      cost.dt[, f.under := f.dt[J(d.under, s.under)]$f]
      cost.dt[, f := f.under+f.over]
      cost.dt[, .(s.under,d.under,f,s.over,d.over)]
      best.cost <- cost.dt[out.dt, {
        min.rows <- .SD[f==min(f)]
        min.rows[.N]#more balanced at end.
      }, keyby=.EACHI, on=.(s.out=s)]
      f.dt[out.dt, `:=`(
        s1=best.cost$s.under, d1=best.cost$d.under,
        s2=best.cost$s.over, d2=best.cost$d.over)]
      best.cost$f
    }
    f.dt[out.dt, f := out.f+out.g]
  }
  ##decoding.
  level <- 0
  new.id <- 1
  new.nodes <- data.table(
    f.dt[.N, .(d,s,parent.x=NA,parent.y=NA,parent=0)])
  while(nrow(new.nodes)){
    node.info <- f.dt[new.nodes][order(parent.x,s)]
    node.info[, id := seq(new.id, new.id+.N-1)]
    node.info[, ord := 1:.N]
    node.info[, y := -level]
    node.info[, x := ord-(.N+1)/2]
    new.id <- node.info[.N, new.id+1]
    node.dt.list[[paste(n.segments, level)]] <- 
      node.info[, data.table(
        n.segments, level,ord,id,d,s,x,y,
        parent.x,parent.y,parent)]
    level <- level+1
    new.nodes <- node.info[, data.table(
      d=c(d1,d2),s=c(s1,s2),parent.x=x,parent.y=y,parent=id
    )][!is.na(d)][order(parent.x)]
  }
  do.call(rbind, node.dt.list)
### Data table with one row for each node in the tree.
}, ex=function(){

  N.data <- 29L
  min.seg.len <- 3L
  (heuristic.df <- binsegRcpp::depth_first_interface(N.data, min.seg.len))
  node.dt <- binsegRcpp::get_best_optimal(N.data, min.seg.len, nrow(heuristic.df))
  (opt.dt <- node.dt[, .(
    candidates=sum(binsegRcpp::size_to_splits(s, min.seg.len))
  ), by=.(parent,level)])

  ## Taking the first few steps of depth first search is not good
  ## enough to get the optimal number of splits. Here is an example
  ## where the depth first search gets four more splits in the first 3
  ## steps.
  node.dt <- binsegRcpp::get_best_optimal(N.data, min.seg.len, 3L)
  (opt.dt <- node.dt[, .(
    candidates=sum(binsegRcpp::size_to_splits(s, min.seg.len))
  ), by=.(parent,level)])


})

get_complexity_extreme <- function
### Compute best and worst case number of splits.
(N.data,
### number of data to segment, positive integer.
  min.segment.length=1L,
### minimum segment length, positive integer.
  n.segments=NULL
### number of segments, positive integer.
){
  node.dt <- get_best_optimal(N.data, min.segment.length, n.segments)
  best.dt <- node.dt[, .(
    candidates=sum(size_to_splits(s, min.segment.length))
  ), by=.(parent,level)]
  first.splits <- size_to_splits(N.data, min.segment.length)
  worst.splits <- c(
    if(1 <= first.splits)seq(first.splits, 1, by=-min.segment.length),
    0)[1:n.segments]
  rbind(
    data.table(
      case="best", 
      segments=1:nrow(best.dt),
      best.dt[, .(splits=candidates, depth=level)]),
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
  ifelse(splits < 0, 0L, splits)
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
  extreme.dt <- get_complexity_extreme(
    n.data, models$min.segment.length, max.segs)
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
    for(data.type in names(worst.counts)){
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

