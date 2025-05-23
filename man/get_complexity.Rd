\name{get_complexity}
\alias{get_complexity}
\title{get complexity}
\description{Get empirical and extreme split counts, in order to compare the
empirical and theoretical time complexity of the binary
segmentation algorithm.}
\usage{get_complexity(models, 
    y.increment = 0.1)}
\arguments{
  \item{models}{result of \code{\link{binseg}}.}
  \item{y.increment}{Offset for y column values of totals output table.}
}

\value{List of class "complexity" which has a plot method. Elements
include "iterations" which is a data table with one row per model
size, and column splits with number of splits to check after
computing that model size; "totals" which is a data table with
total number of splits for each case. }

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre]}




\examples{

## Example 1: empirical=worst case.
data.vec <- rep(0:1, l=8)
plot(data.vec)
worst.model <- binsegRcpp::binseg_normal(data.vec)
worst.counts <- binsegRcpp::get_complexity(worst.model)
plot(worst.counts)

## Example 2: empirical=best case for full path.
data.vec <- 1:8
plot(data.vec)
full.model <- binsegRcpp::binseg_normal(data.vec)
full.counts <- binsegRcpp::get_complexity(full.model)
plot(full.counts)

## Example 3: empirical=best case for all partial paths.
data.vec <- c(0,3,6,10,21,22,23,24)
plot(data.vec)
best.model <- binsegRcpp::binseg_normal(data.vec)
best.counts <- binsegRcpp::get_complexity(best.model)
plot(best.counts)

## ggplot comparing examples 1-3.
if(require("ggplot2")){
  library(data.table)
  splits.list <- list()
  for(data.type in names(worst.counts)){
    splits.list[[data.type]] <- rbind(
      data.table(data="worst", worst.counts[[data.type]]),
      data.table(data="best always", best.counts[[data.type]]),
      data.table(data="best full", full.counts[[data.type]]))
  }
  ggplot()+
    facet_grid(data ~ .)+
    geom_line(aes(
      segments, cum.splits, color=case, size=case),
      data=splits.list$iterations[case!="empirical"])+
    geom_point(aes(
      segments, cum.splits, color=case),
      data=splits.list$iterations[case=="empirical"])+
    scale_color_manual(
      values=binsegRcpp::case.colors,
      breaks=names(binsegRcpp::case.colors))+
    scale_size_manual(
      values=binsegRcpp::case.sizes,
      guide="none")
}

## Example 4: empirical case between best/worst.
data.vec <- rep(c(0,1,10,11),8)
plot(data.vec)
m.model <- binsegRcpp::binseg_normal(data.vec)
m.splits <- binsegRcpp::get_complexity(m.model)
plot(m.splits)

## Example 5: worst case for normal change in mean and variance
## model.
mv.model <- binsegRcpp::binseg("meanvar_norm", data.vec)
mv.splits <- binsegRcpp::get_complexity(mv.model)
plot(mv.splits)

## Compare examples 4-5 using ggplot2.
if(require("ggplot2")){
  library(data.table)
  splits.list <- list()
  for(data.type in names(m.splits)){
    splits.list[[data.type]] <- rbind(
      data.table(model="mean and variance", mv.splits[[data.type]]),
      data.table(model="mean only", m.splits[[data.type]]))
  }
  ggplot()+
    facet_grid(model ~ .)+
    geom_line(aes(
      segments, splits, color=case, size=case),
      data=splits.list$iterations[case!="empirical"])+
    geom_point(aes(
      segments, splits, color=case),
      data=splits.list$iterations[case=="empirical"])+
    geom_text(aes(
      x, y,
      label=label,
      color=case),
      hjust=1,
      data=splits.list$totals)+
    scale_color_manual(
      values=binsegRcpp::case.colors,
      guide="none")+
    scale_size_manual(
      values=binsegRcpp::case.sizes,
      guide="none")
}

## Compare cumsums.
if(require("ggplot2")){
  library(data.table)
  splits.list <- list()
  for(data.type in names(m.splits)){
    splits.list[[data.type]] <- rbind(
      data.table(model="mean and variance", mv.splits[[data.type]]),
      data.table(model="mean only", m.splits[[data.type]]))
  }
  ggplot()+
    facet_grid(model ~ .)+
    geom_line(aes(
      segments, cum.splits, color=case, size=case),
      data=splits.list$iterations[case!="empirical"])+
    geom_point(aes(
      segments, cum.splits, color=case),
      data=splits.list$iterations[case=="empirical"])+
    scale_color_manual(
      values=binsegRcpp::case.colors,
      breaks=names(binsegRcpp::case.colors))+
    scale_size_manual(
      values=binsegRcpp::case.sizes,
      guide="none")
}

}
