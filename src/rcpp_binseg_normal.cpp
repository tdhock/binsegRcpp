/* -*- compile-command: "R -e \"Rcpp::compileAttributes('..')\" && R CMD INSTALL .. && R --vanilla < ../tests/testthat/test_binseg_normal.R" -*- */

#include <Rcpp.h>
#include "binseg_normal.h"
 
// [[Rcpp::export]]
Rcpp::List rcpp_binseg_normal
(const Rcpp::NumericVector &data_vec,
 const Rcpp::IntegerVector &max_segments) {
  int kmax = max_segments[0];
  Rcpp::IntegerVector end(kmax);
  Rcpp::NumericVector loss(kmax);
  Rcpp::NumericVector before_mean(kmax);
  Rcpp::NumericVector after_mean(kmax);
  Rcpp::IntegerVector before_size(kmax);
  Rcpp::IntegerVector after_size(kmax);
  Rcpp::IntegerVector invalidates_index(kmax);
  Rcpp::IntegerVector invalidates_after(kmax);
  int status = binseg_normal
    (&data_vec[0], data_vec.size(), kmax,
     //inputs above, outputs below.
     &end[0], &loss[0],
     &before_mean[0], &after_mean[0],
     &before_size[0], &after_size[0],
     &invalidates_index[0], &invalidates_after[0]);
  if(status == ERROR_NO_DATA){
    Rcpp::stop("no data"); 
  }
  if(status == ERROR_NO_SEGMENTS){
    Rcpp::stop("no segments"); 
  }
  if(status == ERROR_TOO_MANY_SEGMENTS){
    Rcpp::stop("too many segments"); 
  }
  return Rcpp::List::create
    (Rcpp::Named("loss", loss),
     Rcpp::Named("end", end),
     Rcpp::Named("before.mean", before_mean),
     Rcpp::Named("after.mean", after_mean),
     Rcpp::Named("before.size", before_size),
     Rcpp::Named("after.size", after_size),
     Rcpp::Named("invalidates.index", invalidates_index),
     Rcpp::Named("invalidates.after", invalidates_after)
     ) ;
}
