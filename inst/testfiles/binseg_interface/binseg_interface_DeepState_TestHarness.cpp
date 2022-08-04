#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

RInside Rinstance;

/** FUNCTION SIGNATURE */
Rcpp::List binseg_interface(const Rcpp::NumericVector data_vec, const Rcpp::NumericVector weight_vec, 
const int max_segments, const int min_segment_length, const std::string distribution_str, 
const std::string container_str, const Rcpp::LogicalVector is_validation_vec, const Rcpp::NumericVector position_vec);

Rcpp::LogicalVector RcppDeepState_LogicalVector(int rand_size){
  Rcpp::LogicalVector rand_vec(rand_size);
    for(int i = 0 ; i < rand_size;i++){      
      rand_vec[i] = DeepState_IntInRange(0,1);  
    }

  return rand_vec;
}

#define INPUTS \
  int n_data = DeepState_IntInRange(1,1000); \
  Rcpp::NumericVector data_vec = RcppDeepState_NumericVector(n_data, 0, 1000); \
  Rcpp::NumericVector weight_vec = RcppDeepState_NumericVector(n_data, 0, 1000); \
  int max_segments = RcppDeepState_int(1, n_data); \
  int min_segment_length = RcppDeepState_int(1, 10); \
  std::string distribution_str = "test"; \
  std::string container_str = "test"; \
  Rcpp::LogicalVector is_validation_vec = RcppDeepState_LogicalVector(n_data); \
  Rcpp::NumericVector position_vec = RcppDeepState_NumericVector(n_data, 0, 1000); 

TEST(binsegRcpp, generator){
  INPUTS
}

TEST(binsegRcpp, runner){
  INPUTS

  /** INPUTS DUMP
  Skipped: the final analysis table will contain an empty inputs column	
  */

  try{
    binseg_interface(data_vec, weight_vec, max_segments,min_segment_length,
         distribution_str, container_str, is_validation_vec, position_vec);

  }catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
