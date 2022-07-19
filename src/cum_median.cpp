#include <math.h> //isfinite
#include "cum_median.h"
#include "PiecewiseFunction.h"

int cum_median
(int n_data, double *data_vec, double *weight_vec, double *med_vec){
  PiecewiseFunction function;
  for(int data_i=0; data_i < n_data; data_i++){
    double data_value = data_vec[data_i];
    if(!isfinite(data_value)){
      return ERROR_CUM_MEDIAN_DATA_NOT_FINITE;
    }
    double weight_value = weight_vec[data_i];
    function.insert_l1(data_value, weight_value);
    med_vec[data_i] = function.get_minimum_position();
  }
  return 0;
}
