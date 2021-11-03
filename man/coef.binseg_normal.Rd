\name{coef.binseg_normal}
\alias{coef.binseg_normal}
\title{coef binseg normal}
\description{Compute a data table of segment start/end/mean values for all
models given by \code{segments}.}
\usage{\method{coef}{binseg_normal}(object, 
    segments = 1:min(nrow(object$splits), 
        10), ...)}
\arguments{
  \item{object}{data.table from \code{\link{binseg_normal}}.}
  \item{segments}{integer vector, model sizes in number of \code{segments}.}
  \item{\dots}{ignored.}
}

\value{data.table with one row for each segment.}

\author{Toby Dylan Hocking}





