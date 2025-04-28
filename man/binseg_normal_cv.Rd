\name{binseg_normal_cv}
\alias{binseg_normal_cv}
\title{Binary segmentation, normal change in mean, cross-validation for model selection}
\description{Efficient implementation of binary segmentation for change in
mean, with automatic model selection via cross-validation.}
\usage{binseg_normal_cv(data.vec, 
    max.segments = length(data.vec), 
    position.vec = seq_along(data.vec), 
    n.validation.sets = 100L, 
    prop.validation = 0.5)}
\arguments{
  \item{data.vec}{Vector of numeric data to segment.}
  \item{max.segments}{Maximum number of segments to compute, default=length(\code{data.vec}).}
  \item{position.vec}{integer vector of positions at which data are measured,
default=1:length(\code{data.vec}).}
  \item{n.validation.sets}{Number of validation sets.}
  \item{prop.validation}{Proportion of validation set.}
}



\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre]}




\examples{

data.table::setDTthreads(1)

seg.mean.vec <- 1:5
data.mean.vec <- rep(seg.mean.vec, each=20)
set.seed(1)
n.data <- length(data.mean.vec)
data.vec <- rnorm(n.data, data.mean.vec, 0.2)
plot(data.vec)
(fit <- binsegRcpp::binseg_normal_cv(data.vec))
seg.dt <- coef(fit)
model.color <- "red"
seg.dt[, segments(start.pos, mean, end.pos, mean, col=model.color)]
seg.dt[start>1, abline(v=start.pos, col=model.color)]

## plot method shows number of times selected.
plot(fit)

if(requireNamespace("neuroblastoma")){
  data(neuroblastoma, package="neuroblastoma", envir=environment())
  library(data.table)
  profiles.dt <- data.table(neuroblastoma$profiles)
  one.chrom <- profiles.dt[profile.id=="4" & chromosome=="2"]
  fit <- one.chrom[, binsegRcpp::binseg_normal_cv(
    logratio, position.vec=position)]
  selected.segs <- coef(fit)
  if(require(ggplot2)){
    ggplot()+
      geom_point(aes(
        position, logratio),
        data=one.chrom)+
      geom_segment(aes(
        start.pos, mean,
        xend=end.pos, yend=mean),
        data=selected.segs,
        color=model.color)+
      geom_vline(aes(
        xintercept=start.pos),
        data=selected.segs[start>1],
        color=model.color)
  }
}

}
