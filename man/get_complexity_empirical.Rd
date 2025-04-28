\name{get_complexity_empirical}
\alias{get_complexity_empirical}
\title{get complexity empirical}
\description{Get empirical split counts. This is a sub-routine of
\code{\link{get_complexity}}, which should typically be used instead.}
\usage{get_complexity_empirical(model.dt, 
    min.segment.length = 1L)}
\arguments{
  \item{model.dt}{splits data table from \code{\link{binseg}} result list.}
  \item{min.segment.length}{Minimum segment length, positive integer.}
}

\value{data.table with one row per model size, and column splits with
number of splits to check after computing that model size.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre]}





