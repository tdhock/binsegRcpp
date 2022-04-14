\name{get_best_optimal}
\alias{get_best_optimal}
\title{get best optimal}
\description{Dynamic programming for computing lower bound on number of split
candidates to compute / best case of binary segmentation.}
\usage{get_best_optimal(N.data, 
    min.segment.length = 1L, 
    n.segments = NULL)}
\arguments{
  \item{N.data}{positive integer number of data.}
  \item{min.segment.length}{positive integer min segment length.}
  \item{n.segments}{positive integer number of segments.}
}

\value{Data table with one row for each node in the tree.}

\author{Toby Dylan Hocking}





