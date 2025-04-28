\name{coef.binseg_normal_cv}
\alias{coef.binseg_normal_cv}
\title{coef binseg normal cv}
\description{Compute a data table of segment start/end/mean values for all
models given by \code{segments}.}
\usage{\method{coef}{binseg_normal_cv}(object, 
    segments = max(nrow(object$splits)), 
    ...)}
\arguments{
  \item{object}{data.table from \code{\link{binseg_normal_cv}}.}
  \item{segments}{integer vector, model sizes in number of \code{segments}. default=number
of selected \code{segments}.}
  \item{\dots}{ignored.}
}

\value{data.table with one row for each segment.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre]}





