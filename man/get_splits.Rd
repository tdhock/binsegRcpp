\name{get_splits}
\alias{get_splits}
\title{get splits}
\description{Get empirical and extreme split counts.}
\usage{get_splits(model.dt)}
\arguments{
  \item{model.dt}{data.table from \code{\link{binseg_normal}}.}
}

\value{data.table with one row per model size, and column splits with
number of splits to check after computing that model size. Column
case has values best, worst, empirical.}

\author{Toby Dylan Hocking}




\examples{

## Example 1: empirical=worst case.
data.vec <- rep(c(0,1), l=10)
bs.model <- binsegRcpp::binseg_normal(data.vec)
split.counts <- binsegRcpp::get_splits(bs.model)
plot(split.counts)

## Example 2: empirical=best case.
data.vec <- 1:20
bs.model <- binsegRcpp::binseg_normal(data.vec)
split.counts <- binsegRcpp::get_splits(bs.model)
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
