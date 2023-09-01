\name{get_complexity_best_optimal_tree}
\alias{get_complexity_best_optimal_tree}
\title{get complexity best optimal tree}
\description{decoding.}
\usage{get_complexity_best_optimal_tree(f.dt)}
\arguments{
  \item{f.dt}{f.dt }
}

\value{Data table with one row for each node in the tree.}

\author{Toby Dylan Hocking}




\examples{

N.data <- 19L
min.seg.len <- 3L
max.segments <- 4L
cost.dt <- binsegRcpp::get_complexity_best_optimal_cost(
  N.data, min.seg.len, max.segments)
binsegRcpp::get_complexity_best_optimal_tree(cost.dt)

}
