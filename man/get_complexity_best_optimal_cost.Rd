\name{get_complexity_best_optimal_cost}
\alias{get_complexity_best_optimal_cost}
\title{get complexity best optimal cost}
\description{Dynamic programming for computing lower bound on number of split
candidates to compute / best case of binary segmentation. The
dynamic programming recursion is on f(d,s) = best number of splits
for segment of size s which is split d times. Need to optimize
f(d,s) = g(s) + min f(d1,s1) + f(d2,s2) over s1,d1 given that
s1+s2=s, d1+d2+1=d, and g(s) is the number of splits for segment
of size s.}
\usage{get_complexity_best_optimal_cost(N.data, 
    min.segment.length = 1L, 
    n.segments = NULL)}
\arguments{
  \item{N.data}{positive integer number of data.}
  \item{min.segment.length}{positive integer min segment length.}
  \item{n.segments}{positive integer number of segments.}
}

\value{data table with one row for each f(d,s) value computed.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre]}




\examples{

binsegRcpp::get_complexity_best_optimal_cost(
  N.data = 19L, 
  min.segment.length = 3L, 
  n.segments = 4L)

}
