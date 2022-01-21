#include <Rcpp.h>
#include <R.h>
#include "binseg_normal.h"
 
// [[Rcpp::export]]
Rcpp::List rcpp_binseg_normal
(const Rcpp::NumericVector data_vec,
 const int kmax,
 const Rcpp::LogicalVector is_validation_vec,
 const Rcpp::NumericVector position_vec
 ) {
  int n_data = data_vec.size();
  if(n_data < 1){
    Rcpp::stop("need at least one data point"); 
  }
  if(is_validation_vec.size() != n_data){
    Rcpp::stop("length of is_validation_vec must be same as data_vec");
  }
  if(position_vec.size() != n_data){
    Rcpp::stop("length of position_vec must be same as data_vec");
  }
  //TODO stop for non-finite data.
  int n_subtrain = get_n_subtrain(n_data, &is_validation_vec[0]);
  if(n_subtrain == 0){
    Rcpp::stop("need at least one subtrain data");
  }
  if(kmax < 1){
    Rcpp::stop("kmax must be positive"); 
  }
  Rcpp::IntegerVector end(kmax);
  Rcpp::NumericVector subtrain_borders(n_subtrain+1);
  Rcpp::NumericVector loss(kmax);
  Rcpp::NumericVector validation_loss(kmax);
  Rcpp::NumericVector before_mean(kmax);
  Rcpp::NumericVector after_mean(kmax);
  Rcpp::IntegerVector before_size(kmax);
  Rcpp::IntegerVector after_size(kmax);
  Rcpp::IntegerVector invalidates_index(kmax);
  Rcpp::IntegerVector invalidates_after(kmax);
  int status = binseg_normal 
    (&data_vec[0], n_data, kmax, &is_validation_vec[0], &position_vec[0],
     //inputs above, outputs below.
     &end[0], &subtrain_borders[0], &loss[0], &validation_loss[0],
     &before_mean[0], &after_mean[0],
     &before_size[0], &after_size[0],
     &invalidates_index[0], &invalidates_after[0]);
  //Rprintf("status=%d\n", status);
  if(status == ERROR_TOO_MANY_SEGMENTS){
    Rcpp::stop("too many segments"); 
  }
  if(status == ERROR_POSITIONS_MUST_INCREASE){
    Rcpp::stop("positions must increase");
  }
  return Rcpp::List::create
    (Rcpp::Named("loss", loss),
     Rcpp::Named("validation.loss", validation_loss),
     Rcpp::Named("end", end),
     Rcpp::Named("subtrain.borders", subtrain_borders),
     Rcpp::Named("before.mean", before_mean),
     Rcpp::Named("after.mean", after_mean),
     Rcpp::Named("before.size", before_size),
     Rcpp::Named("after.size", after_size),
     Rcpp::Named("invalidates.index", invalidates_index),
     Rcpp::Named("invalidates.after", invalidates_after)
     ) ;
}

