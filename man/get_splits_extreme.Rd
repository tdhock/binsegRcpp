\name{get_splits_extreme}
\alias{get_splits_extreme}
\title{get splits extreme}
\description{Compute best and worst case number of splits.}
\usage{get_splits_extreme(N.data)}
\arguments{
  \item{N.data}{number of data to segment.}
}

\value{data.table with one row per model size, and column splits with
number of splits to check after computing that model size. Column
case has values best (equal segment sizes, min splits to check)
and worst (unequal segment sizes, max splits to check).}

\author{Toby Dylan Hocking}





