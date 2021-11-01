\name{binseg_normal}
\alias{binseg_normal}
\title{Binary segmentation, normal change in mean}
\description{Efficient implementation of binary segmentation for change in
mean, max normal likelihood = min square loss. Output includes
columns which can be used to compute parameters for a single model
in log-linear time.}
\usage{binseg_normal(data.vec, 
    max.segments = length(data.vec), 
    is.validation.vec = rep(FALSE, 
        length(data.vec)), 
    position.vec = seq_along(data.vec))}
\arguments{
  \item{data.vec}{Vector of numeric data to segment.}
  \item{max.segments}{Maximum number of segments to compute, default=length(\code{data.vec}).}
  \item{is.validation.vec}{logical vector indicating which data are to be used in validation
set, default=all FALSE (no validation set).}
  \item{position.vec}{integer vector of positions at which data are measured,
default=1:length(\code{data.vec}).}
}

\value{data.table with a row for each model and columns
\item{segments}{number of parameters}
\item{loss}{square loss}
\item{end}{index of last data point per segment}
\item{before.mean}{mean before changepoint}
\item{after.mean}{mean after changepoint}
\item{before.size}{number of data before changepoint}
\item{after.size}{number of data after changepoint}
\item{invalidates.index}{index of model parameter no longer used after this changepoint is used}
\item{invalidates.after}{idem}}

\author{Toby Dylan Hocking}




\examples{

x <- c(0.1, 0, 1, 1.1, 0.1, 0)
## Compute full path of binary segmentation models from 1 to 6
## segments.
(models.dt <- binsegRcpp::binseg_normal(x))

## Plot loss values using base graphics.
plot(models.dt)

## Same loss values using ggplot2.
if(require("ggplot2")){
  ggplot()+
    geom_point(aes(
      segments, loss),
      data=models.dt)
}

## Compute data table of segments to plot.
(segs.dt <- coef(models.dt, 2:4))

## Plot data, segments, changepoints.
if(require("ggplot2")){
  ggplot()+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    facet_grid(segments ~ ., labeller=label_both)+
    geom_vline(aes(
      xintercept=start-0.5),
      color="green",
      data=segs.dt[1<start])+
    geom_segment(aes(
      start-0.5, mean,
      xend=end+0.5, yend=mean),
      data=segs.dt,
      color="green")+
    xlab("Position/index")+
    ylab("Data/mean value")+
    geom_point(aes(
      pos, x),
      data=data.frame(x, pos=seq_along(x)))
}

}
