\name{get_complexity_empirical}
\alias{get_complexity_empirical}
\title{get complexity empirical}
\description{Get empirical split counts.}
\usage{get_complexity_empirical(model.dt, 
    min.segment.length = 1L)}
\arguments{
  \item{model.dt}{data.table from \code{\link{binseg_normal}}.}
  \item{min.segment.length}{Minimum segment length, positive integer.}
}

\value{data.table with one row per model size, and column splits with
number of splits to check after computing that model size.}

\author{Toby Dylan Hocking}





