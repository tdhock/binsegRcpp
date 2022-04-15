\name{tree_layout}
\alias{tree_layout}
\title{tree layout}

\usage{tree_layout(node.dt)}
\arguments{
  \item{node.dt}{
}
}



\author{Toby Dylan Hocking}




\examples{

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

}
