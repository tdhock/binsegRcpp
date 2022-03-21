#include <Rcpp.h>
#include <R.h>
#include "binseg.h"
 
//' Lookup the string values used to represent different distributions
//'
//' @return Character vector corresponding to supported distributions
// [[Rcpp::export]]
Rcpp::CharacterVector get_distribution_names(){
  map_type *dmap = get_dist_map();
  int n_items = dmap->size();
  Rcpp::CharacterVector names(n_items);
  int i=0;
  for(map_type::iterator it=dmap->begin(); it != dmap->end(); it++){
    names[i++] = it->first;
  }
  return names;
}

// [[Rcpp::export]]
Rcpp::List binseg_interface
(const Rcpp::NumericVector data_vec,
 const Rcpp::NumericVector weight_vec,
 const int kmax,
 const std::string distribution_str,
 const Rcpp::LogicalVector is_validation_vec,
 const Rcpp::NumericVector position_vec
 ) {
  int n_data = data_vec.size();
  if(n_data < 1){
    Rcpp::stop("need at least one data point"); 
  }
  if(weight_vec.size() != n_data){
    Rcpp::stop("length of weight_vec must be same as data_vec");
  }
  if(is_validation_vec.size() != n_data){
    Rcpp::stop("length of is_validation_vec must be same as data_vec");
  }
  if(position_vec.size() != n_data){
    Rcpp::stop("length of position_vec must be same as data_vec");
  }
  if(Rcpp::any(!Rcpp::is_finite(data_vec))){
    Rcpp::stop("data must be finite");
  }
  int n_subtrain = get_n_subtrain(n_data, &is_validation_vec[0]);
  if(n_subtrain == 0){
    Rcpp::stop("need at least one subtrain data");
  }
  if(kmax < 1){
    Rcpp::stop("kmax must be positive"); 
  }
  Rcpp::NumericVector subtrain_borders(n_subtrain+1);
  Rcpp::IntegerVector end(kmax);
  Rcpp::NumericVector loss(kmax);
  Rcpp::NumericVector validation_loss(kmax);
  Rcpp::NumericVector before_mean(kmax);
  Rcpp::NumericVector after_mean(kmax);
  Rcpp::IntegerVector before_size(kmax);
  Rcpp::IntegerVector after_size(kmax);
  Rcpp::IntegerVector invalidates_index(kmax);
  Rcpp::IntegerVector invalidates_after(kmax);
  int status = binseg 
    (&data_vec[0], &weight_vec[0],
     n_data, kmax, &is_validation_vec[0], &position_vec[0],
     distribution_str.c_str(),
     //inputs above, outputs below.
     &subtrain_borders[0],
     &end[0], &loss[0], &validation_loss[0],
     &before_mean[0], &after_mean[0],
     &before_size[0], &after_size[0],
     &invalidates_index[0], &invalidates_after[0]);
  if(status == ERROR_UNRECOGNIZED_DISTRIBUTION){
    Rcpp::stop("unrecognized distribution"); 
  }
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

