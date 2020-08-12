#include <Rcpp.h>
#include <R.h>
#include "binseg_normal.h"
#include "binseg_normal_cost.h"

// [[Rcpp::export]]
int read_vector(int i){
  Rcpp::IntegerVector x(0);
  return x[i];
}

// [[Rcpp::export]]
int read_new(int i){
  int* ptr = new int[0];
  int x = ptr[i];
  delete[] ptr;
  return x;
}

// [[Rcpp::export]]
int read_malloc(int i){
  int *ptr = (int*)malloc(0);
  int x = ptr[i];
  free(ptr);
  return x;
}

// [[Rcpp::export]]
Rcpp::List rcpp_binseg_normal
(Rcpp::NumericVector data_vec,
 Rcpp::IntegerVector max_segments) {
  int kmax = max_segments[0];
  max_segments[10000000] = 5;
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

// [[Rcpp::export]]
Rcpp::List rcpp_binseg_normal_cost
(const Rcpp::NumericVector data_vec,
 const Rcpp::IntegerVector max_segments) {
  int kmax = max_segments[0];
  Rcpp::NumericVector loss(kmax);
  binseg_normal_cost
    (&data_vec[0], data_vec.size(), kmax,
     //inputs above, outputs below.
     &loss[0]);
  return Rcpp::List::create(Rcpp::Named("loss", loss));
}
