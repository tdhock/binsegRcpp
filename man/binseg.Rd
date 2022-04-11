\name{binseg}
\alias{binseg}
\alias{binsegRcpp}
\title{Binary segmentation}
\description{Efficient C++ implementation of the classic binary segmentation
algorithm for finding changepoints in a sequence of N data. Output
includes columns which can be used to compute parameters for a
single model in log-linear time, using coef method.}
\usage{binseg(distribution.str, 
    data.vec, max.segments = NULL, 
    is.validation.vec = rep(FALSE, 
        length(data.vec)), 
    position.vec = seq_along(data.vec), 
    weight.vec = rep(1, 
        length(data.vec)), 
    min.segment.length = NULL, 
    container.str = "multiset")}
\arguments{
  \item{distribution.str}{String indicating distribution, use \code{\link{get_distribution_info}} to see
possible values.}
  \item{data.vec}{Vector of numeric data to segment.}
  \item{max.segments}{Maximum number of segments to compute, default=NULL which means to
compute the largest number possible, given \code{is.validation.vec} and
\code{min.segment.length}. Note that the returned number of segments may
be less than this, if there are min segment length constraints.}
  \item{is.validation.vec}{logical vector indicating which data are to be used in validation
set, default=all FALSE (no validation set).}
  \item{position.vec}{integer vector of positions at which data are measured,
default=1:length(\code{data.vec}).}
  \item{weight.vec}{Numeric vector of non-negative weights for each data point.}
  \item{min.segment.length}{Positive integer, minimum number of data points per
segment. Default NULL means to use min given \code{distribution.str}.}
  \item{container.str}{C++ container to use for storing breakpoints/cost. Most users
should leave this at the default "multiset" for efficiency but you
could use "list" if you want to study the time complexity of a
slower implementation of binary segmentation.}
}
\details{Each iteration involves first computing and storing the
best split point on one or two segments, then looking up the
segment with the best split so far. The best case time complexity
occurs when splits are equal (N data split into two segments of
size N/2), and the worst case is when splits are unequal (N data
split into one big segment with N-1 data and one small segment
with 1 data point). Looking up the segment with the best split so
far is a constant O(1) time operation using C++ multimap, so O(K)
overall for K iterations/segments. Storage of a new best split
point/cost involves the multimap insert method which is
logarithmic time in the size of the multimap, overall O(K log K)
for equal splits and O(K) for unequal splits. Computing the cost
values, and overall time complexity, depends on the loss. For
normal and poisson distributions the best case O(N log K) time
for equal splits and worst case O(N K) time for unequal
splits. For l1/laplace distributions the best case is O(N log N
log K) time for equal splits and worst case is O(N log N K) time
for unequal splits.}
\value{list of class binsegRcpp with elements \code{min.segment.length},
\code{distribution.str}, param.names, subtrain.borders and splits, which
is a data.table with columns:
\item{segments}{number of segments}
\item{loss}{subtrain loss}
\item{validation.loss}{validation loss}
\item{end}{index of last data point per segment}
\item{depth}{number of splits to reach segment}
\item{before}{params before changepoint}
\item{after}{params after changepoint}
\item{before.size}{number of data before changepoint}
\item{after.size}{number of data after changepoint}
\item{invalidates.index}{index of param invalidated by this split.}
\item{invalidates.after}{indicates if before/after params invalidated by this split.}}

\author{Toby Dylan Hocking}




\examples{

x <- c(0.1, 0, 1, 1.1, 0.1, 0)
## Compute full path of binary segmentation models from 1 to 6
## segments.
(models <- binsegRcpp::binseg("mean_norm", x))

## Plot loss values using base graphics.
plot(models)

## Same loss values using ggplot2.
if(require("ggplot2")){
  ggplot()+
    geom_point(aes(
      segments, loss),
      data=models$splits)
}

## Compute data table of segments to plot.
(segs.dt <- coef(models, 2:4))

## Plot data, segments, changepoints.
if(require("ggplot2")){
  ggplot()+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    facet_grid(segments ~ ., labeller=label_both)+
    geom_vline(aes(
      xintercept=start.pos),
      color="green",
      data=segs.dt[1<start])+
    geom_segment(aes(
      start.pos, mean,
      xend=end.pos, yend=mean),
      data=segs.dt,
      color="green")+
    xlab("Position/index")+
    ylab("Data/mean value")+
    geom_point(aes(
      pos, x),
      data=data.frame(x, pos=seq_along(x)))
}

## Use min.segment.length to constrain segment sizes.
(constrained.models <- binsegRcpp::binseg("mean_norm", x, min.segment.length = 2L))

## Demonstration of model selection using cross-validation in
## simulated data.
seg.mean.vec <- 1:5
data.mean.vec <- rep(seg.mean.vec, each=20)
set.seed(1)
n.data <- length(data.mean.vec)
data.vec <- rnorm(n.data, data.mean.vec, 0.2)
plot(data.vec)

library(data.table)
loss.dt <- data.table(seed=1:10)[, {
  set.seed(seed)
  is.valid <- sample(rep(c(TRUE,FALSE), l=n.data))
  bs.model <- binsegRcpp::binseg("mean_norm", data.vec, is.validation.vec=is.valid)
  bs.model$splits[, data.table(
    segments,
    validation.loss)]
}, by=seed]
loss.stats <- loss.dt[, .(
  mean.valid.loss=mean(validation.loss)
), by=segments]
plot(
  mean.valid.loss ~ segments, loss.stats,
  col=ifelse(
    mean.valid.loss==min(mean.valid.loss),
    "black",
    "red"))

selected.segments <- loss.stats[which.min(mean.valid.loss), segments]
full.model <- binsegRcpp::binseg("mean_norm", data.vec, selected.segments)
(segs.dt <- coef(full.model, selected.segments))
if(require("ggplot2")){
  ggplot()+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    geom_vline(aes(
      xintercept=start.pos),
      color="green",
      data=segs.dt[1<start])+
    geom_segment(aes(
      start.pos, mean,
      xend=end.pos, yend=mean),
      data=segs.dt,
      color="green")+
    xlab("Position/index")+
    ylab("Data/mean value")+
    geom_point(aes(
      pos, data.vec),
      data=data.frame(data.vec, pos=seq_along(data.vec)))
}

## Demo of poisson loss, weights.
data.vec <- c(3,4,10,20)
(fit1 <- binsegRcpp::binseg("poisson", data.vec, weight.vec=c(1,1,1,10)))
coef(fit1, 2L)
(fit2 <- binsegRcpp::binseg("poisson", data.vec, weight.vec=c(1,1,10,1)))
coef(fit2, 2L)

}
