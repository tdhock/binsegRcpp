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

### Checks types and values of size inputs.
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

get_complexity_worst <- function
### Get full sequence of splits which results in worst case time
### complexity.
(N.data, min.segment.length){
  first.splits <- size_to_splits(N.data, min.segment.length)
  c(if(1 <= first.splits)seq(first.splits, 1, by=-min.segment.length), 0)
}

get_complexity_best_heuristic_equal_depth_full <- function
### Heuristic depth first.
(N.data, min.segment.length){
  depth_first_interface(N.data, min.segment.length)
}

get_complexity_best_heuristic_equal_breadth_full <- function
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

get_tree_empirical <- function
### Compute tree for empirical binary segmentation model.
(fit){
  before.size <- . <- segments <- invalidates.index <-
    invalidates.after <- after.size <- parent <- J <- NULL
  child.dt.list <- list(
    data.table(id=0L,size=fit$splits[1,before.size],parent=NA))
  get_id <- function(idx,after)as.integer(2*idx-after-2)
  children.wide <- fit$splits[-1, .(
    segments,
    parent=get_id(invalidates.index,invalidates.after),
    before.size, after.size)]
  for(child in c("before", "after")){
    child.name <- paste0(child, ".size")
    size <- children.wide[[child.name]]
    child.dt.list[[child]] <- children.wide[, data.table(
      id=get_id(segments,child=="after"),
      size,
      parent)]
  }
  child.dt <- do.call(rbind, child.dt.list)
  depth <- 0
  setkey(child.dt, parent)
  new.parents <- child.dt[1]
  tree.dt.list <- list()
  while(nrow(new.parents)){
    new.children <- child.dt[J(new.parents$id), nomatch=0L]
    tree.dt.list[[paste(depth)]] <- data.table(depth, new.parents)
    new.parents <- new.children
    depth <- depth+1
  }
  do.call(rbind, tree.dt.list)
}

qp.x <- function
### Solve quadratic program to find x positions.
(target, y.up, y.lo){
  k <- length(target)
  D <- diag(rep(1, k))
  Ik <- diag(rep(1, k - 1))
  A <- rbind(0, Ik) - rbind(Ik, 0)
  b0 <- (y.up - target)[-k] + (target - y.lo)[-1]
  sol <- quadprog::solve.QP(D, target, A, b0)
  sol$solution
}

tree_layout <- structure(function
### Compute x,y coordinates for graphing a tree.
(node.dt, space=0.5){
  x <- id <- depth <- J <- parent.x <- size <-
    parent.depth <- parent <- NULL
  stopifnot(identical(names(node.dt), c("depth","id","size","parent")))
  id.tab <- table(node.dt$id)
  stopifnot(all(id.tab==1))
  tree.dt <- data.table(node.dt)
  tree.dt[, x := NA_real_]
  setkey(tree.dt, id)
  for(d in unique(tree.dt$depth)){
    if(d==0)tree.dt[depth==0, x := 0] else{
      d.nodes <- tree.dt[depth==d]
      px <- tree.dt[J(d.nodes$parent), x, nomatch=0L]
      d.nodes[, parent.x := px]
      ord.nodes <- d.nodes[order(parent.x, size)]
      new.x <- ord.nodes[, qp.x(parent.x,parent.x+space,parent.x-space)]
      tree.dt[J(ord.nodes$id), x := new.x]
    }
  }
  px <- tree.dt[J(tree.dt$parent), x]
  tree.dt[, parent.x := px]
  tree.dt[, parent.depth := ifelse(is.na(parent), NA, depth-1)]
  tree.dt
}, ex=function(){

  N.data <- 29L
  min.seg.len <- 3L
  max.segments <- 5L
  cost.dt <- binsegRcpp::get_complexity_best_optimal_cost(
    N.data, min.seg.len, max.segments)
  set.seed(1)
  data.vec <- rnorm(N.data)
  fit <- binsegRcpp::binseg_normal(data.vec, max.segments)
  tree.list <- list(
    best=binsegRcpp::get_complexity_best_optimal_tree(cost.dt),
    empirical=binsegRcpp::get_tree_empirical(fit))
  library(data.table)
  tree.dt <- data.table(type=names(tree.list))[, {
    binsegRcpp::tree_layout(tree.list[[type]])
  }, by=type]
  total.dt <- tree.dt[, .(
    candidate.splits=sum(binsegRcpp::size_to_splits(size, min.seg.len))
  ), by=type]
  join.dt <- total.dt[tree.dt, on="type"]
  if(require(ggplot2)){
    ggplot()+
      facet_grid(. ~ type + candidate.splits, labeller=label_both)+
      geom_segment(aes(
        x, depth, 
        xend=parent.x, yend=parent.depth),
        data=join.dt)+
      geom_label(aes(
        x, depth, label=size),
        data=join.dt)+
      scale_y_reverse()
  }

})

get_complexity_best_optimal_tree <- structure(function
### decoding.
(f.dt){
  . <- d <- s <- id <- ord <- y <- x <- parent <- d1 <-
    d2 <- s1 <- s2 <- NULL
  level <- 0
  new.id <- 0
  new.nodes <- data.table(
    f.dt[.N, .(d,s,parent=NA)])
  node.dt.list <- list()
  while(nrow(new.nodes)){
    node.info <- f.dt[new.nodes]
    node.info[, id := seq(new.id, new.id+.N-1)]
    node.info[, ord := 1:.N]
    node.info[, y := -level]
    node.info[, x := ord-(.N+1)/2]
    new.id <- node.info[.N, id+1]
    node.dt.list[[paste(level)]] <- 
      node.info[, data.table(
        depth=level,id,size=s,parent)]
    level <- level+1
    new.nodes <- node.info[, data.table(
      d=c(d1,d2),s=c(s1,s2),parent=id
    )][!is.na(d)]
  }
  do.call(rbind, node.dt.list)
### Data table with one row for each node in the tree.
}, ex=function(){

  N.data <- 19L
  min.seg.len <- 3L
  max.segments <- 4L
  cost.dt <- binsegRcpp::get_complexity_best_optimal_cost(
    N.data, min.seg.len, max.segments)
  binsegRcpp::get_complexity_best_optimal_tree(cost.dt)

})

get_complexity_best_optimal_splits <- function
### Convert output of get_complexity_best_optimal_tree to counts of
### candidate splits that need to be considered at each iteration.
(node.dt, min.segment.length){
  . <- size <- parent <- depth <- NULL
  node.dt[, .(
    splits=sum(size_to_splits(size, min.segment.length))
  ), by=.(parent,depth)]
### Data table with one row for each segment.
}

get_complexity_best_optimal_cost <- structure(function
### Dynamic programming for computing lower bound on number of split
### candidates to compute / best case of binary segmentation. The
### dynamic programming recursion is on f(d,s) = best number of splits
### for segment of size s which is split d times. Need to optimize
### f(d,s) = g(s) + min f(d1,s1) + f(d2,s2) over s1,d1 given that
### s1+s2=s, d1+d2+1=d, and g(s) is the number of splits for segment
### of size s.
(N.data,
### positive integer number of data.
  min.segment.length=1L,
### positive integer min segment length.
  n.segments=NULL
### positive integer number of segments.
){
  d <- s <- J <- d.under <- s.out <- s.over <- s.under <- 
    f.over <- f.under <- f <- . <- parent.x <- id <- 
      ord <- y <- x <- parent.y <- parent <- d1 <- d2 <- 
        s1 <- s2 <- NULL
  check_sizes(N.data, min.segment.length, n.segments)
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
  f.dt
### data table with one row for each f(d,s) value computed.
}, ex=function(){

  binsegRcpp::get_complexity_best_optimal_cost(
    N.data = 19L, 
    min.segment.length = 3L, 
    n.segments = 4L)

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
  . <- s <- parent <- splits <- depth <- NULL
  cost.dt <- get_complexity_best_optimal_cost(
    N.data, min.segment.length, n.segments)
  tree.dt <- get_complexity_best_optimal_tree(cost.dt)
  best.dt <- get_complexity_best_optimal_splits(tree.dt, min.segment.length)
  worst.splits <- get_complexity_worst(
    N.data, min.segment.length)[1:n.segments]
  rbind(
    data.table(
      case="best", 
      segments=1:nrow(best.dt),
      best.dt[, .(splits, depth)]),
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
  segments <- end <- . <- splits <- y <- label <- case <- 
    cum.splits <- NULL
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

