\name{get_complexity_best_optimal_tree}
\alias{get_complexity_best_optimal_tree}
\title{get complexity best optimal tree}
\description{decoding.}
\usage{get_complexity_best_optimal_tree(f.dt)}
\arguments{
  \item{f.dt}{
}
}

\value{Data table with one row for each node in the tree.}

\author{Toby Dylan Hocking}




\examples{

N.data <- 19L
min.seg.len <- 3L
max.segments <- 4L
(heuristic.df <- binsegRcpp::depth_first_interface(N.data, min.seg.len))
cost.dt <- binsegRcpp::get_complexity_best_optimal_cost(
  N.data, min.seg.len, nrow(heuristic.df))
node.dt <- binsegRcpp::get_complexity_best_optimal_tree(cost.dt)
set.seed(1)
data.vec <- rnorm(N.data)
fit <- binsegRcpp::binseg_normal(data.vec, max.segments)
child.dt.list <- list(
  data.table(id=0,size=fit$splits[1,before.size],parent=NA))
get_id <- function(idx,after)2*idx-after-2
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
level <- 0
setkey(child.dt, parent)
new.parents <- child.dt[1]
tree.dt.list <- list()
while(nrow(new.parents)){
  new.children <- child.dt[J(new.parents$id), nomatch=0L]
  tree.dt.list[[paste(level)]] <- data.table(level, new.parents)
  new.parents <- new.children
  level <- level+1
}
tree.dt <- do.call(rbind, tree.dt.list)

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
}
