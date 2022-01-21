\name{get_complexity}
\alias{get_complexity}
\title{get complexity}
\description{Get empirical and extreme split counts.}
\usage{get_complexity(models)}
\arguments{
  \item{models}{data.table from \code{\link{binseg_normal}}.}
}

\value{data.table with one row per model size, and column splits with
number of splits to check after computing that model size. Column
case has values best, worst, empirical.}

\author{Toby Dylan Hocking}




\examples{

## Example 1: empirical=worst case.
data.vec <- rep(c(0,1), l=10)
plot(data.vec)
bs.model <- binsegRcpp::binseg_normal(data.vec)
split.counts <- binsegRcpp::get_complexity(bs.model)
plot(split.counts)

## Example 2: empirical=best case.
data.vec <- 1:20
plot(data.vec)
bs.model <- binsegRcpp::binseg_normal(data.vec)
split.counts <- binsegRcpp::get_complexity(bs.model)
plot(split.counts)

## Example 3: empirical case between best/worst.
seg.mean.vec <- 1:5
data.mean.vec <- rep(seg.mean.vec, each=10)
set.seed(1)
data.vec <- rnorm(length(data.mean.vec), data.mean.vec, 0.2)
plot(data.vec)
bs.model <- binsegRcpp::binseg_normal(data.vec)
split.counts <- binsegRcpp::get_complexity(bs.model)
plot(split.counts)

if(require("ggplot2")){
  ggplot()+
    geom_line(aes(
      segments, splits, color=case, size=case),
      data=split.counts$iterations[case!="empirical"])+
    geom_point(aes(
      segments, splits, color=case),
      data=split.counts$iterations[case=="empirical"])+
    geom_text(aes(
      x, y,
      label=label,
      color=case),
      hjust=1,
      data=split.counts$totals)+
    scale_color_manual(
      values=binsegRcpp::case.colors,
      guide="none")+
    scale_size_manual(
      values=binsegRcpp::case.sizes,
      guide="none")
}

}