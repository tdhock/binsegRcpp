#include <Rcpp.h>
#include "binseg.h"
 
//' Lookup the string values used to represent different distributions
// [[Rcpp::export]]
Rcpp::CharacterVector get_distribution_names(){
  dist_map_type *dmap = get_dist_map();
  int n_items = dmap->size();
  Rcpp::CharacterVector names(n_items);
  int i=0;
  for(dist_map_type::iterator it=dmap->begin(); it != dmap->end(); it++){
    names[i++] = it->first;
  }
  return names;
}

template < typename T >
std::string unrecognized
(std::string what, T* (*get_map)(void)){
  std::string msg = "unrecognized ";
  msg += what;
  msg += ", try one of: ";
  T *map = get_map();
  typename T::iterator it=map->begin();
  while(1){
    msg += it->first;
    if(++it != map->end()){
      msg += ", ";
    }else break;
  }
  return msg;
}

//' Low-level interface to binary segmentation algorithm.
// [[Rcpp::export]]
Rcpp::List binseg_interface
(const Rcpp::NumericVector data_vec,
 const Rcpp::NumericVector weight_vec,
 const int max_segments,
 const int min_segment_length,
 const std::string distribution_str,
 const std::string container_str,
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
  if(max_segments < 1){//keep below one subtrain error msg.
    Rcpp::stop("max_segments must be positive"); 
  }
  param_names_type *param_names_ptr;
  try{
    param_names_ptr = get_param_names(distribution_str.c_str());
  }
  catch(const std::out_of_range& err){
    Rcpp::stop(unrecognized<dist_map_type>("distribution", get_dist_map));
  }
  int n_params = param_names_ptr->size();
  Rcpp::CharacterVector param_names_vec(n_params);
  int param_i = 0;
  for(param_names_type::iterator it=param_names_ptr->begin();
      it != param_names_ptr->end(); it++){
    param_names_vec[param_i++] = *it;
  }
  Rcpp::NumericVector subtrain_borders(n_subtrain+1);
  Rcpp::IntegerVector end(max_segments);
  Rcpp::NumericVector loss(max_segments);
  Rcpp::NumericVector validation_loss(max_segments);
  Rcpp::NumericMatrix before_param_mat(max_segments, n_params);
  Rcpp::colnames(before_param_mat) = param_names_vec;
  Rcpp::NumericMatrix after_param_mat(max_segments, n_params);
  Rcpp::colnames(after_param_mat) = param_names_vec;
  Rcpp::IntegerVector before_size(max_segments);
  Rcpp::IntegerVector after_size(max_segments);
  Rcpp::IntegerVector invalidates_index(max_segments);
  Rcpp::IntegerVector invalidates_after(max_segments);
  int status = binseg 
    (&data_vec[0], &weight_vec[0],
     n_data, max_segments, min_segment_length,
     &is_validation_vec[0], &position_vec[0],
     distribution_str.c_str(),
     container_str.c_str(),
     //inputs above, outputs below.
     &subtrain_borders[0],
     &end[0], &loss[0], &validation_loss[0],
     &before_param_mat[0], &after_param_mat[0],
     &before_size[0], &after_size[0],
     &invalidates_index[0], &invalidates_after[0]);
  if(status == ERROR_UNRECOGNIZED_CONTAINER){
    Rcpp::stop(unrecognized<factory_map_type>("container", get_factory_map));
  }
  if(status == ERROR_TOO_MANY_SEGMENTS){
    Rcpp::stop("too many segments, max_segments=%d and min_segment_length=%d which would require at least %d data but n_subtrain=%d", max_segments, min_segment_length, max_segments*min_segment_length, n_subtrain); 
  }
  if(status == ERROR_MIN_SEGMENT_LENGTH_MUST_BE_POSITIVE){
    Rcpp::stop("min segment length must be positive"); 
  }
  if(status == ERROR_POSITIONS_MUST_INCREASE){
    Rcpp::stop("positions must increase");
  }
  return Rcpp::List::create
    (Rcpp::Named("subtrain.borders", subtrain_borders),
     Rcpp::Named("end", end),
     Rcpp::Named("loss", loss),
     Rcpp::Named("validation.loss", validation_loss),
     Rcpp::Named("before.param.mat", before_param_mat),
     Rcpp::Named("after.param.mat", after_param_mat),
     Rcpp::Named("before.size", before_size),
     Rcpp::Named("after.size", after_size),
     Rcpp::Named("invalidates.index", invalidates_index),
     Rcpp::Named("invalidates.after", invalidates_after)
     ) ;
}

