\name{get_complexity_empirical}
\alias{get_complexity_empirical}
\title{get complexity empirical}
\description{Get empirical split counts.}
\usage{get_complexity_empirical(model.dt)}
\arguments{
  \item{model.dt}{data.table from \code{\link{binseg_normal}}.}
}

\value{data.table with one row per model size, and column splits with
number of splits to check after computing that model size.}

\author{Toby Dylan Hocking}




