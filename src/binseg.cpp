#include "binseg.h"

double Cumsum::get_sum(int first, int last){
  double total = cumsum_vec[last];
  if(0 < first){
    total -= cumsum_vec[first-1];
  }
  return total;
}

void Set::set_mean_var_loss
(int first, int last, double *mean, double *var, double *loss){
  total_weights = weights.get_sum(first, last);
  total_weighted_data = weighted_data.get_sum(first, last);
  total_weighted_squares = weighted_squares.get_sum(first, last);
  *mean = total_weighted_data/total_weights;
  *var  = total_weighted_squares/total_weights +
    (*mean)*((*mean)-2*total_weighted_data/total_weights);
  *loss = get_loss(first, last, *mean, *var);
}
void Set::set_mean_var_loss(int first, int last, MeanVarLoss *MVL){
  set_mean_var_loss(first, last, &(MVL->mean), &(MVL->var), &(MVL->loss));
}
double Set::get_loss
(int first, int last, MeanVarLoss& subtrain_mvl){
  return get_loss
    (first, last, subtrain_mvl.mean, subtrain_mvl.var);
}
double Set::get_loss
(int first, int last, double mean, double var){
  total_weights = weights.get_sum(first, last);
  total_weighted_data = weighted_data.get_sum(first, last);
  total_weighted_squares = weighted_squares.get_sum(first, last);
  return dist_ptr->compute_loss
    (total_weights, total_weighted_data, total_weighted_squares,
     mean, var, max_zero_var);
}
void Set::resize_cumsums(int vec_size){
  weights.cumsum_vec.resize(vec_size);
  weighted_data.cumsum_vec.resize(vec_size);
  weighted_squares.cumsum_vec.resize(vec_size);
}
void Set::write_cumsums(int write_index){
  weights.cumsum_vec[write_index] = total_weights;
  weighted_data.cumsum_vec[write_index] = total_weighted_data;
  weighted_squares.cumsum_vec[write_index] = total_weighted_squares;
}

static dist_map_type dist_map;
// we need get_dist_map in this file to query map contents from
// interface.cpp
dist_map_type* get_dist_map(void){
  return &dist_map;
}

#define CONCAT(x,y) x##y
#define DISTRIBUTION(NAME, DESC, COMPUTE, VARIANCE)                     \
  class CONCAT(NAME,Distribution) : public Distribution {               \
  public:                                                               \
    double compute_loss                                                 \
      (double N, double sum, double squares,				\
       double mean, double var, double max_zero_var){			\
      return COMPUTE;                                                   \
    }									\
    CONCAT(NAME,Distribution)                                           \
      (const char *name, std::string desc, bool var_changes){           \
      description = desc;                                               \
      param_names_vec.push_back("mean");                                \
      if(var_changes)param_names_vec.push_back("var");                  \
      dist_map.emplace(name, this);                                     \
    }                                                                   \
  };                                                                    \
  static CONCAT(NAME,Distribution) NAME( #NAME, DESC, VARIANCE );

#define RSS (mean*(N*mean-2*sum)+squares)

DISTRIBUTION(mean_norm,
             "change in normal mean with constant variance (L2/square loss)",
             RSS, 
	     false) 
/* Above we compute the square loss for a segment with sum of data = s
   and mean parameter m.

   If x_i is data point i, and s = sum_i x_i is the sum over all N
   data points on that segment, then total loss is

   sum_i (m - x_i)^2 = N*m^2 - 2*m*s + sum_i x_i^2

   The last term (sum of squares of data) can be ignored during
   optimization, because it does not depend on the optimization
   variable (segment mean). It is added back after optimization, at
   the end of binseg.

   Including weights it becomes

   sum_i w_i (M - x_i)^2 = M^2 (sum_i w_i) - 2*M*(sum_i w_i x_i)
   + sum_i w_i x_i^2

   Ignoring last term it becomes

   = W*M^2 - 2*M*S = M*(W*M - 2*S), where

   W = sum_i w_i (N in code), S = sum_i w_i x_i (sum), M = S/W (mean).
   
*/

DISTRIBUTION(poisson,
             "change in poisson rate parameter (loss is negative log likelihood minus constant term)",
             mean*N - log(mean)*sum, // neg log lik minus constant term.
             false) // dont add constant term to loss.
/* poisson likelihood:

prob_i = m^{x_i} * exp(-m) / (x_i !)

log(prob_i) = x_i*log(m) - m - x_i!

-log(prob_i) = w_i*x_i*log(M) - w_i*M - w_i*(x_i!)
   
poisson loss with weights:

     sum_i w_i [M - x_i log(M)] = M*(sum_i w_i) - log(M)*(sum_i x_i w_i)

  = M*W - log(M)*S.
  
 */

DISTRIBUTION(meanvar_norm,
             "change in normal mean and variance (loss is negative log likelihood)",
	     (var>max_zero_var) ? (RSS/var+N*log(2*M_PI*var))/2 : INFINITY,
	     true)
/*
meanvar_norm loss is negative log likelihood =

0.5 [ (sum_i x_i^2 + M(NM-2 sum_i x_i))/var + log(2*pi*var) ]
 */

double Split::set_mean_var_loss(Set &subtrain, int first, int end_i, int last){
  this_end = end_i;
  subtrain.set_mean_var_loss(first, end_i, &before);
  subtrain.set_mean_var_loss(end_i+1, last, &after);
  return before.loss + after.loss;
}

int Segment::n_changes() const {
  return last_i-first_i;
}
Segment::Segment
(Set &subtrain, Set &validation,
 int first_data, int last_data,
 int first_candidate, int last_candidate,
 int invalidates_after, int invalidates_index,
 double loss_no_split, double validation_loss_no_split
 ): first_i(first_data), last_i(last_data),
    invalidates_index(invalidates_index),
    invalidates_after(invalidates_after){
  Split candidate_split;
  double best_loss_split = INFINITY, loss_split;
  // for loop over all possible splits on this Segment.
  for(int candidate=first_candidate; candidate<=last_candidate; candidate++){
    loss_split = candidate_split.set_mean_var_loss
      (subtrain, first_data, candidate, last_data);
    if(loss_split < best_loss_split){
      best_loss_split = loss_split;
      best_split = candidate_split;
    }
  }
  best_decrease = best_loss_split - loss_no_split;
  if(best_decrease == INFINITY)return;
  before_validation_loss = validation.get_loss
    (first_data, best_split.this_end, best_split.before);
  after_validation_loss = validation.get_loss
    (best_split.this_end+1, last_data, best_split.after);
  double validation_loss_split =
    before_validation_loss + after_validation_loss;
  validation_decrease =
    validation_loss_split - validation_loss_no_split;
}

template <typename T>
class MyContainer : public Container {
public:
  T segment_container;
  typename T::iterator best;
  int get_size(void){
    return segment_container.size();
  }
  void remove_best(void){
    segment_container.erase(best);
  }
  virtual typename T::iterator get_best_it(void) = 0;  
  const Segment* set_best(void){
    best = get_best_it();
    return &(*best);
  }
};

typedef std::multiset<Segment> segment_set_type;
static factory_map_type factory_map;
ContainerFactory::ContainerFactory
(const char *name, construct_fun_type construct, destruct_fun_type destruct){
  construct_fun_ptr = construct;
  destruct_fun_ptr = destruct;
  factory_map.emplace(name, this);
}
factory_map_type* get_factory_map(void){
  return &factory_map;
}

#define CMAKER(CONTAINER, INSERT, BEST) \
  class CONCAT(CONTAINER,Wrapper) : public MyContainer< std::CONTAINER<Segment> > { \
  public:                                                               \
    void insert(Segment& new_seg){                                      \
      segment_container.INSERT(new_seg);                                \
    }                                                                   \
    std::CONTAINER<Segment>::iterator get_best_it(void){                \
      return BEST;                                                      \
    }                                                                   \
  };                                                                    \
  Container* CONCAT(CONTAINER,construct) (){                            \
    return new CONCAT(CONTAINER,Wrapper);                               \
  }                                                                     \
  void CONCAT(CONTAINER,destruct) (Container *c_ptr){                   \
    delete c_ptr;           \
  }                                                                     \
  static ContainerFactory CONCAT(CONTAINER,_instance)                   \
    ( #CONTAINER, CONCAT(CONTAINER,construct), CONCAT(CONTAINER,destruct) );

CMAKER(multiset, insert, segment_container.begin())

CMAKER(list, push_back, std::min_element(segment_container.begin(),segment_container.end()))

class Candidates {
public:
  ContainerFactory *factory_ptr;
  Container *container_ptr = 0;
  Set subtrain, validation;
  int min_segment_length;
  ~Candidates(){
    if(container_ptr != 0)factory_ptr->destruct_fun_ptr(container_ptr);
  }
  // computes the cumulative sum vectors in linear O(n_data) time.
  int init
  (const char *container_str,
   const double *data_vec, const double *weight_vec, const int n_data,
   const double *position_vec, const int *is_validation_vec,
   double *subtrain_borders, Distribution *dist_ptr,
   int min_segment_length_arg
   ){
    factory_ptr = factory_map.at(container_str);
    container_ptr = factory_ptr->construct_fun_ptr();
    min_segment_length = min_segment_length_arg;
    subtrain.dist_ptr = dist_ptr;
    validation.dist_ptr = dist_ptr;
    int n_validation = 0;
    for(int data_i=0; data_i<n_data; data_i++){
      if(is_validation_vec[data_i]){
	n_validation++;
      }
    }
    int n_subtrain = n_data - n_validation;
    subtrain.resize_cumsums(n_subtrain);
    validation.resize_cumsums(n_subtrain);
    int last_subtrain_i=-1;
    double pos_total, pos_change;
    int read_start=0, write_index=0;
    for(int data_i=0; data_i<=n_data; data_i++){
      bool is_subtrain = false;
      bool write_subtrain = false;
      bool write_end = data_i == n_data;
      if(!write_end){
	is_subtrain = !is_validation_vec[data_i];
	write_subtrain = last_subtrain_i >= 0 && is_subtrain;
      }
      if(write_subtrain || write_end){
	if(write_subtrain){
	  pos_total = position_vec[data_i]+position_vec[last_subtrain_i];
	  pos_change = pos_total/2;
	  if(write_index==0){
	    subtrain_borders[write_index] = position_vec[last_subtrain_i]-0.5;
	  }
	}else{
	  pos_change = position_vec[data_i-1]+0.5;//last.
	}
	subtrain_borders[write_index+1] = pos_change;
	int read_index=read_start;
	while(read_index < n_data && position_vec[read_index] <= pos_change){
	  double data_value = data_vec[read_index];
          double weight_value = weight_vec[read_index];
          Set *this_set;
	  if(is_validation_vec[read_index]){
            this_set = &validation;
          }else{
            this_set = &subtrain;
          }
          this_set->total_weights +=
	    weight_value;
          this_set->total_weighted_data +=
	    data_value * weight_value;
          this_set->total_weighted_squares +=
            data_value * data_value * weight_value;
	  read_index++;
	}
        validation.write_cumsums(write_index);
        subtrain.write_cumsums(write_index);
	read_start = read_index;
	write_index++;
      }
      if(is_subtrain){
	last_subtrain_i = data_i;
      }
    }
    // analyze subtrain data to find what is the largest value that
    // should be considered numerically zero for a variance estimate.
    MeanVarLoss mvl;
    subtrain.max_zero_var = 0;
    for(int subtrain_i=0; subtrain_i < n_subtrain; subtrain_i++){
      subtrain.set_mean_var_loss(subtrain_i, subtrain_i, &mvl);
      if(subtrain.max_zero_var < mvl.var){
	subtrain.max_zero_var = mvl.var;
      }
    }
    validation.max_zero_var = subtrain.max_zero_var;
    return n_subtrain;
  }
  // Add a new Segment to candidates if it is big enough to split.
  void maybe_add
  (int first_data, int last_data,
   int invalidates_after, int invalidates_index,
   double loss_no_split, double validation_loss_no_split
   ){
    int first_candidate = first_data + min_segment_length-1;
    int last_candidate = last_data - min_segment_length;
    if(first_candidate <= last_candidate){
      // if it is possible to split then insert new segment into the
      // candidates set.
      Segment new_seg
	(subtrain, validation, 
	 first_data, last_data,
         first_candidate, last_candidate,
	 invalidates_after, invalidates_index,
	 loss_no_split, validation_loss_no_split);
      if(new_seg.best_decrease < INFINITY){
        container_ptr->insert(new_seg);
      }
    }
  }
};

class OutArrays {
public:
  int n_params, max_segments;
  int *seg_end, *before_size, *after_size,
    *invalidates_index, *invalidates_after;
  double *subtrain_loss, *validation_loss,
    *before_param_mat, *after_param_mat;
  OutArrays
  (Distribution *dist_ptr, int max_segments_,
   int *seg_end_, double *subtrain_loss_, double *validation_loss_,
   double *before_param_mat_, double *after_param_mat_,
   int *before_size_, int *after_size_,
   int *invalidates_index_, int *invalidates_after_){
    n_params = dist_ptr->param_names_vec.size();
    max_segments = max_segments_;
    seg_end = seg_end_;
    subtrain_loss = subtrain_loss_;
    validation_loss = validation_loss_;
    before_param_mat = before_param_mat_;
    after_param_mat = after_param_mat_;
    before_size = before_size_;
    after_size = after_size_;
    invalidates_index = invalidates_index_;
    invalidates_after = invalidates_after_;
  }
  void save
  (int seg_i,
   double subtrain_loss_value,
   double validation_loss_value,
   int seg_end_value,
   const MeanVarLoss &before_mvl,
   const MeanVarLoss &after_mvl, 
   int invalidates_index_value,
   int invalidates_after_value,
   int before_size_value,
   int after_size_value
   ){
    subtrain_loss[seg_i] = subtrain_loss_value;
    validation_loss[seg_i] = validation_loss_value;
    seg_end[seg_i] = seg_end_value;
    before_param_mat[seg_i] = before_mvl.mean;
    after_param_mat[seg_i] = after_mvl.mean;
    if(n_params == 2){
      before_param_mat[seg_i+max_segments] = before_mvl.var;
      after_param_mat[seg_i+max_segments] = after_mvl.var;
    }
    invalidates_index[seg_i] = invalidates_index_value;
    invalidates_after[seg_i] = invalidates_after_value;
    before_size[seg_i] = before_size_value;
    after_size[seg_i]  = after_size_value;
  }
};

/* Binary segmentation algorithm.

   This code assumes all array parameters are allocated by the code
   which calls this function, and the code which calls this function
   needs to have error checking for, the following:

   At least one data point (0 < n_data), arrays of this size:
   - data_vec: input data sequence to segment,
   - weight_vec: non-negative weight for each data point,
   - is_validation_vec: indicator for validation set,
   - position_vec: position where each data point is measured in time/space,

   distribution_str: string indicating loss function to use.
   subtrain_borders: array of size n_subtrain+1.

   Positive number of segments (0 < max_segments), arrays of this size:
   - seg_end: end of segment/split.
   - subtrain_loss: subtrain loss of each split.
   - validation_loss: validation loss of each split.
   - before_param_mat: params before split.
   - after_param_mat: params after split.
   - before_size: number of data before this split.
   - after_size: number of data after this split.
   - invalidates_index: index of params invalidated by this split.
   - invalidates_after: indicates if before/after params invalidated by this split.

   See coef method in R code for a procedure that uses these output
   arrays to efficiently compute the segment means for any model size.
 */
int binseg
(const double *data_vec, const double *weight_vec,
 const int n_data, const int max_segments, const int min_segment_length,
 const int *is_validation_vec, const double *position_vec,
 const char *distribution_str,
 const char *container_str,
 double *subtrain_borders, 
 int *seg_end, double *subtrain_loss, double *validation_loss,
 double *before_param_mat, double *after_param_mat, 
 int *before_size, int *after_size,  
 int *invalidates_index, int *invalidates_after
 ){
  if(min_segment_length < 1){
    return ERROR_MIN_SEGMENT_LENGTH_MUST_BE_POSITIVE;
  }
  for(int data_i=1; data_i<n_data; data_i++){
    if(position_vec[data_i] <= position_vec[data_i-1]){
      return ERROR_POSITIONS_MUST_INCREASE;
    }
  }
  Distribution* dist_ptr = dist_map.at(distribution_str);
  OutArrays out_arrays
    (dist_ptr, max_segments,
     seg_end, subtrain_loss, validation_loss,
     before_param_mat, after_param_mat,
     before_size, after_size,
     invalidates_index, invalidates_after);
  Candidates V;
  int n_subtrain;
  try{
    n_subtrain = V.init
    (container_str,
     data_vec, weight_vec, n_data, position_vec, is_validation_vec,
     subtrain_borders, dist_ptr, min_segment_length);
  }
  catch(const std::out_of_range& err){
    return ERROR_UNRECOGNIZED_CONTAINER;
  }
  if(n_subtrain < max_segments*min_segment_length){
    return ERROR_TOO_MANY_SEGMENTS;
  }
  // Then store the trivial segment mean/loss (which starts at the
  // first and ends at the last data point).
  MeanVarLoss full_mvl, missing_mvl;
  missing_mvl.mean = INFINITY;
  missing_mvl.var = INFINITY;
  V.subtrain.set_mean_var_loss(0, n_subtrain-1, &full_mvl);
  out_arrays.save
    (0,
     full_mvl.loss,
     V.validation.get_loss(0, n_subtrain-1, full_mvl),
     n_subtrain-1,
     full_mvl,
     missing_mvl,
     -2, -2, n_subtrain, -2);
  // Add a segment and split to the set of candidates.
  V.maybe_add(0, n_subtrain-1, 0, 0, subtrain_loss[0], validation_loss[0]);
  // initialize to infinite cost and missing values, which is
  // necessary when number of segments returned is less than
  // max_segments: infinite cost rows are removed from the resulting
  // splits table in the R code.
  for(int seg_i=1; seg_i < max_segments; seg_i++){
    out_arrays.save
      (seg_i, INFINITY, INFINITY, -2, full_mvl, full_mvl, -2, -2, -2, -2);
  }
  // Loop over splits. During each iteration we find the Segment/split
  // which results in the best loss decrease, store the resulting
  // model parameters, and add new Segment/split candidates if
  // necessary.
  int seg_i = 0;
  while(V.container_ptr->not_empty() && ++seg_i < max_segments){
    // Store loss and model parameters associated with this split.
    const Segment *seg_ptr = V.container_ptr->set_best();
    out_arrays.save
      (seg_i, 
       subtrain_loss[seg_i-1] + seg_ptr->best_decrease,
       validation_loss[seg_i-1] + seg_ptr->validation_decrease,
       seg_ptr->best_split.this_end,
       seg_ptr->best_split.before,
       seg_ptr->best_split.after,
       seg_ptr->invalidates_index,
       seg_ptr->invalidates_after,
       seg_ptr->best_split.this_end - seg_ptr->first_i + 1,
       seg_ptr->last_i - seg_ptr->best_split.this_end);
    // Finally add new split candidates if necessary.
    V.maybe_add
      (seg_ptr->first_i, seg_ptr->best_split.this_end,
       0,//invalidates_after=0 => before_mean invalidated.
       seg_i, seg_ptr->best_split.before.loss,
       seg_ptr->before_validation_loss);
    V.maybe_add
      (seg_ptr->best_split.this_end+1, seg_ptr->last_i,
       1,//invalidates_after=1 => after_mean invalidated.
       seg_i, seg_ptr->best_split.after.loss,
       seg_ptr->after_validation_loss);
    // Erase at end because we need seg_ptr->values during maybe_add
    // inserts above.
    V.container_ptr->remove_best();
  }
  return 0;//SUCCESS.
}

int get_n_subtrain
(const int n_data, const int *is_validation_vec){
  int n_subtrain = 0;
  for(int data_i=0; data_i<n_data; data_i++){
    if(!is_validation_vec[data_i]){
      n_subtrain++;
    }
  }
  return n_subtrain;
}
  
param_names_type* get_param_names
(const char *distribution_str){
  return &(dist_map.at(distribution_str)->param_names_vec);
}
