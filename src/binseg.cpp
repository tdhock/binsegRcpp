#include "binseg.h"
#include <cmath>
#include "PiecewiseFunction.h"

Split::Split() {}

void Split::set_end_dist(int first, int end, int last){
  // first=1, last=4, possible end values are 1-3, 2 is best.
  this_end = end;
  int dist_first = end-first;
  int dist_last = last-end-1;
  dist_from_edges = (dist_first < dist_last) ? dist_first : dist_last;
}

void Split::maybe_update(Split *candidate_ptr){
  double loss_diff = get_loss() - candidate_ptr->get_loss();
  bool bigger_dist = 
    loss_diff == 0 && dist_from_edges < candidate_ptr->dist_from_edges;
  if(loss_diff > 0 || bigger_dist){
    *this = *candidate_ptr;
  }
}

void Set::set_max_zero_var(void){
  max_zero_var = dist_ptr->get_max_zero_var(*this);
}

double Cumsum::get_sum(int first, int last){
  double total = cumsum_vec[last];
  if(0 < first){
    total -= cumsum_vec[first-1];
  }
  return total;
}

void Set::set_totals(int first, int last){
  total_weights = weights.get_sum(first, last);
  total_weighted_data = weighted_data.get_sum(first, last);
  total_weighted_squares = weighted_squares.get_sum(first, last);
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

static dist_umap_type dist_umap;
// we need get_dist_umap in this file to query map contents from
// interface.cpp
dist_umap_type* get_dist_umap(void){
  return &dist_umap;
}

class CumDistribution : public Distribution {
public:
  void set_best_split
  (Split *best_split_ptr, Set &subtrain,
   int first_data, int last_data,
   int first_candidate, int last_candidate,
   Split *candidate_split_ptr
   ){
    for(int candidate=first_candidate; candidate<=last_candidate; candidate++){
      candidate_split_ptr->set_end_dist(first_data, candidate, last_data);
      //CumDistribution::estimate_params is O(1) so set_best_split is
      //linear in the number of candidates.
      estimate_params
        (&candidate_split_ptr->before, subtrain, first_data, candidate);
      estimate_params
        (&candidate_split_ptr->after, subtrain, candidate+1, last_data);
      best_split_ptr->maybe_update(candidate_split_ptr);
    }
  }
  double get_max_zero_var(Set &subtrain){
    int n_subtrain = subtrain.weights.cumsum_vec.size();
    double max_zero_var = 0;
    ParamsLoss ploss;
    for(int subtrain_i=0; subtrain_i < n_subtrain; subtrain_i++){
      subtrain.dist_ptr->estimate_params
        (&ploss, subtrain, subtrain_i, subtrain_i);
      if(max_zero_var < ploss.spread){
	max_zero_var = ploss.spread;
      }
    }
    return max_zero_var;
  }    
  void estimate_params
  (ParamsLoss *ploss_ptr, Set &subtrain, int first, int last){
    subtrain.set_totals(first, last);
    ploss_ptr->center = subtrain.total_weighted_data/subtrain.total_weights;
    if(subtrain.dist_ptr->var_param){
      ploss_ptr->spread = 
        subtrain.total_weighted_squares/subtrain.total_weights +
        ploss_ptr->center*
        (ploss_ptr->center-2*
         subtrain.total_weighted_data/subtrain.total_weights);
    }
    ploss_ptr->loss = compute_loss
      (subtrain.total_weights,
       subtrain.total_weighted_data,
       subtrain.total_weighted_squares,
       ploss_ptr->center,
       ploss_ptr->spread,
       subtrain.max_zero_var);
  }
  virtual double compute_loss(double,double,double,double,double,double) = 0;
  double loss_for_params
  (Set &validation, ParamsLoss &ploss, int first, int last){
    validation.set_totals(first, last);
    return compute_loss
      (validation.total_weights,
       validation.total_weighted_data,
       validation.total_weighted_squares,
       ploss.center,
       ploss.spread,
       validation.max_zero_var);
  }
};  

#define CONCAT(x,y) x##y
#define CUM_DIST(NAME, DESC, COMPUTE, ERROR, VARIANCE)			\
  class CONCAT(NAME,Distribution) : public CumDistribution {            \
  public:                                                               \
    int check_data(double value){					\
      ERROR;								\
      return 0;								\
    }									\
    double compute_loss                                                 \
      (double N, double sum, double squares,				\
       double mean, double var, double max_zero_var){			\
      return COMPUTE;                                                   \
    }									\
    CONCAT(NAME,Distribution) () {					\
      var_param = VARIANCE;						\
      description = DESC;                                               \
      param_names_vec.push_back("mean");                                \
      if(var_param)param_names_vec.push_back("var");                  \
      dist_umap.emplace( #NAME, this );					\
    }                                                                   \
  };                                                                    \
  static CONCAT(NAME,Distribution) NAME;

class absDistribution : public Distribution {
  public:
  int check_data(double value){
    return 0;
  }
  double adjust(double sum_abs_dev, double N, double scale_est){
    if(var_param==false)return sum_abs_dev;
    if(scale_est==0)return INFINITY;
    return N*log(2*scale_est) + sum_abs_dev/scale_est;
  }
  double get_max_zero_var(Set &subtrain){
    return 0;
  }
  void set_best_split
  (Split *best_split_ptr, Set &subtrain,
   int first_data, int last_data,
   int first_candidate, int last_candidate,
   Split *candidate_split_ptr
   ){
    int n_candidates = last_candidate-first_candidate+1;
    int n_insertions = last_candidate-first_data+1;
    std::vector<double> before_median_vec(n_candidates);
    std::vector<double> before_loss_vec(n_candidates);
    std::vector<double> before_weight_vec(n_candidates);
    std::vector<double> after_median_vec(n_candidates);
    std::vector<double> after_loss_vec(n_candidates);
    std::vector<double> after_weight_vec(n_candidates);
    for(int direction=0; direction<2; direction++){
      PiecewiseFunction function;
      int start, increment, offset;
      std::vector<double> *loss_ptr, *median_ptr, *weight_ptr;
      if(direction==0){
        start = first_data;
        increment = 1;
        loss_ptr = &before_loss_vec;
        median_ptr = &before_median_vec;
        weight_ptr = &before_weight_vec;
        offset = 0;
      }else{
        start = last_data;
        increment = -1;
        loss_ptr = &after_loss_vec;
        median_ptr = &after_median_vec;
        weight_ptr = &after_weight_vec;
        offset = 1;
      }
      int out_i = 0;
      double total_weight = 0;
      for(int iteration=0; iteration < n_insertions; iteration++){
        int data_i = start + iteration*increment;
        double weight_value = subtrain.weights.get_sum(data_i,data_i);
        total_weight += weight_value;
        double weighted_data = subtrain.weighted_data.get_sum(data_i,data_i);
        double data_value = weighted_data/weight_value;
        function.insert_l1(data_value, weight_value);
        if
          (first_candidate+offset <= data_i &&
           data_i <= last_candidate+offset){
          (*median_ptr)[out_i] = function.get_minimum_position();
          (*loss_ptr)[out_i] = function.get_minimum_value();
          (*weight_ptr)[out_i] = total_weight;
          out_i++;
        }
      }//for(iteration
    }//for(direction
    //now that all *pred_vec and *loss_vec entries have been computed,
    //we can find the split with min loss.
    for(int before_i=0; before_i<n_candidates; before_i++){
      int after_i = n_candidates-1-before_i;
      candidate_split_ptr->set_end_dist
        (first_data, first_candidate+before_i, last_data);
      candidate_split_ptr->before.center = 
        before_median_vec[before_i];
      candidate_split_ptr->after.center = 
        after_median_vec[after_i];
      candidate_split_ptr->before.spread =
        before_loss_vec[before_i]/before_weight_vec[before_i];
      candidate_split_ptr->after.spread =
        after_loss_vec[after_i]/after_weight_vec[after_i];
      candidate_split_ptr->before.loss = adjust
        (before_loss_vec[before_i],
         before_weight_vec[before_i],
         candidate_split_ptr->before.spread);
      candidate_split_ptr->after.loss = adjust
        (after_loss_vec[after_i],
         after_weight_vec[after_i],
         candidate_split_ptr->after.spread);
      best_split_ptr->maybe_update(candidate_split_ptr);
    }
  }
  double loss_for_params
  (Set &validation, ParamsLoss &ploss, int first, int last){
    double total_loss=0, total_weight=0;
    double median = ploss.center;
    // linear O(last-first).
    for(int data_i=first; data_i <= last; data_i++){
      double weight_value = validation.weights.get_sum(data_i,data_i);
      if(0 < weight_value){
        total_weight += weight_value;
        double weighted_data = validation.weighted_data.get_sum(data_i,data_i);
        double data_value = weighted_data/weight_value;
        total_loss += abs(median - data_value)*weight_value;
      }
    }
    return adjust(total_loss, total_weight, ploss.spread);
  }
  void estimate_params
  (ParamsLoss *ploss_ptr, Set &subtrain, int first, int last){
    // log-linear in number of data between first and last.
    PiecewiseFunction function;
    double total_weight = 0;
    for(int data_i=first; data_i <= last; data_i++){
      double weight_value = subtrain.weights.get_sum(data_i,data_i);
      double weighted_data = subtrain.weighted_data.get_sum(data_i,data_i);
      double data_value = weighted_data/weight_value;
      function.insert_l1(data_value, weight_value);
      total_weight += weight_value;
    }
    ploss_ptr->center = function.get_minimum_position();
    double sum_abs_dev = function.get_minimum_value();
    ploss_ptr->spread = sum_abs_dev/total_weight;
    ploss_ptr->loss = adjust(sum_abs_dev, total_weight, ploss_ptr->spread);
  }
};

#define ABS_DIST(NAME, DESC, VARIANCE)                          \
  class CONCAT(NAME, Distribution) : public absDistribution {   \
  public:                                                       \
    CONCAT(NAME, Distribution) (){                              \
      var_param = VARIANCE;                                   \
      description = DESC;                                       \
      param_names_vec.push_back("median");                      \
      if(var_param)param_names_vec.push_back("scale");        \
      dist_umap.emplace( #NAME, this );                          \
    }                                                           \
  };                                                            \
static CONCAT(NAME, Distribution) NAME;

ABS_DIST(l1, "change in median (loss is total absolute deviation)", false)

ABS_DIST(laplace, "change in Laplace median and scale (loss is negative log likelihood)", true)

#define RSS (mean*(N*mean-2*sum)+squares)

CUM_DIST(mean_norm,
         "change in normal mean with constant variance (L2/square loss)",
         RSS,
	 ,
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

CUM_DIST(poisson,
         "change in Poisson rate parameter (loss is negative log likelihood minus constant term)",
         (mean>0) ? (mean*N - log(mean)*sum) : ( (sum==0) ? 0 : INFINITY ),
	 if(round(value)!=value)return ERROR_DATA_MUST_BE_INTEGER_FOR_POISSON_LOSS;if(value < 0)return ERROR_DATA_MUST_BE_NON_NEGATIVE_FOR_POISSON_LOSS;,
         false)
/* poisson likelihood:

prob_i = m^{x_i} * exp(-m) / (x_i !)

log(prob_i) = x_i*log(m) - m - x_i!

-log(prob_i) = w_i*x_i*log(M) - w_i*M - w_i*(x_i!)
   
poisson loss with weights:

     sum_i w_i [M - x_i log(M)] = M*(sum_i w_i) - log(M)*(sum_i x_i w_i)

  = M*W - log(M)*S.
  
 */

CUM_DIST(meanvar_norm,
         "change in normal mean and variance (loss is negative log likelihood)",
         (var>max_zero_var) ? (RSS/var+N*log(2*M_PI*var))/2 : INFINITY,
	 ,
         true)
/*
meanvar_norm loss is negative log likelihood =

0.5 [ (sum_i x_i^2 + M(NM-2 sum_i x_i))/var + log(2*pi*var) ]
 */

Segment::Segment
(Set &subtrain, Set &validation,
 int first_data, int last_data,
 int first_candidate, int last_candidate,
 int invalidates_after, int invalidates_index,
 double loss_no_split, double validation_loss_no_split,
 int depth, Split *candidate_split_ptr
 ): first_i(first_data), last_i(last_data),
  depth(depth),
    invalidates_index(invalidates_index),
    invalidates_after(invalidates_after){
  subtrain.dist_ptr->set_best_split
    (&best_split, subtrain, 
     first_data, last_data, 
     first_candidate, last_candidate,
     candidate_split_ptr);
  best_decrease = best_split.get_loss() - loss_no_split;
  if(best_decrease == INFINITY)return;
  before_validation_loss = validation.dist_ptr->loss_for_params
    (validation, best_split.before, first_data, best_split.this_end);
  after_validation_loss = validation.dist_ptr->loss_for_params
    (validation, best_split.after, best_split.this_end+1, last_data);
  double validation_loss_split =
    before_validation_loss + after_validation_loss;
  validation_decrease =
    validation_loss_split - validation_loss_no_split;
}

template <typename T>
class MyContainer : public Container {
public:
  T segment_container;
  int get_size(void){
    return segment_container.size();
  }
};

typedef std::multiset<Segment> segment_set_type;
static container_umap_type container_umap;
ContainerFactory::ContainerFactory
(const char *name, construct_fun_type construct, destruct_fun_type destruct){
  construct_fun_ptr = construct;
  destruct_fun_ptr = destruct;
  container_umap.emplace(name, this);
}
container_umap_type* get_container_umap(void){
  return &container_umap;
}

#define CMAKER(CONTAINER, STRUCT, INSERT, SET_IT, GET_SEG, ERASE)	\
  class CONCAT(CONTAINER,Wrapper) : public MyContainer< STRUCT > {	\
  public:                                                               \
    void insert(Segment& new_seg){                                      \
      segment_container.INSERT(new_seg);                                \
    }                                                                   \
    Segment get_best(void){						\
      SET_IT;								\
      Segment seg = GET_SEG;						\
      ERASE;								\
      return seg;							\
    }                                                                   \
  };                                                                    \
  Container* CONCAT(CONTAINER,construct) (){                            \
    return new CONCAT(CONTAINER,Wrapper);                               \
  }                                                                     \
  void CONCAT(CONTAINER,destruct) (Container *c_ptr){                   \
    delete c_ptr;                                                       \
  }                                                                     \
  static ContainerFactory CONCAT(CONTAINER,_instance)                   \
    ( #CONTAINER, CONCAT(CONTAINER,construct), CONCAT(CONTAINER,destruct) );

#define CIT(CONTAINER, INSERT, BEST)				\
  CMAKER(CONTAINER, std::CONTAINER<Segment>, INSERT, std::CONTAINER<Segment>::iterator it = BEST, *it, segment_container.erase(it))

CIT(multiset, insert, segment_container.begin())

CIT(list, push_back, std::min_element(segment_container.begin(),segment_container.end()))

class PQ_Compare {
public:
  bool operator()(Segment a, Segment b){
    return !(a < b);
  }
};
#define PQ_STRUCT std::priority_queue<Segment,std::vector<Segment>,PQ_Compare>
CMAKER(priority_queue, PQ_STRUCT, push, , segment_container.top(), segment_container.pop())
//template <class T, class Container = vector<T>,  class Compare = less<typename Container::value_type> >
//template < class T,                        // multiset::key_type/value_type           class Compare = less<T>,        // multiset::key_compare/value_compare           class Alloc = allocator<T> >    // multiset::allocator_type           > class multiset;

class Candidates {
public:
  ContainerFactory *factory_ptr;
  Container *container_ptr = 0;
  Set subtrain, validation;
  int min_segment_length;
  Split candidate_split;
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
    factory_ptr = container_umap.at(container_str);
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
    if(dist_ptr->var_param){
      subtrain.set_max_zero_var();
      validation.max_zero_var = subtrain.max_zero_var;
    }
    return n_subtrain;
  }
  // Add a new Segment to candidates if it is big enough to split.
  void maybe_add
  (int first_data, int last_data,
   int invalidates_after, int invalidates_index,
   double loss_no_split, double validation_loss_no_split,
   int depth
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
	 loss_no_split, validation_loss_no_split,
         depth+1, &candidate_split);
      if(new_seg.best_decrease < INFINITY){
        container_ptr->insert(new_seg);
      }
    }
  }
};

class OutArrays {
public:
  Distribution *dist_ptr;
  int max_segments;
  int *seg_end, *depth, *before_size, *after_size,
    *invalidates_index, *invalidates_after;
  double *subtrain_loss, *validation_loss,
    *before_param_mat, *after_param_mat;
  OutArrays
  (Distribution *dist_ptr_, int max_segments_,
   int *seg_end_, int *depth_, 
   double *subtrain_loss_, double *validation_loss_,
   double *before_param_mat_, double *after_param_mat_,
   int *before_size_, int *after_size_,
   int *invalidates_index_, int *invalidates_after_){
    dist_ptr = dist_ptr_;
    max_segments = max_segments_;
    seg_end = seg_end_;
    depth = depth_;
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
   int depth_value,
   const ParamsLoss &before_ploss,
   const ParamsLoss &after_ploss, 
   int invalidates_index_value,
   int invalidates_after_value,
   int before_size_value,
   int after_size_value
   ){
    subtrain_loss[seg_i] = subtrain_loss_value;
    validation_loss[seg_i] = validation_loss_value;
    seg_end[seg_i] = seg_end_value;
    depth[seg_i] = depth_value;
    int param_i=0;
    for
      (param_names_type::iterator it=dist_ptr->param_names_vec.begin();
       it != dist_ptr->param_names_vec.end();
       it++){
      int out_i = seg_i+max_segments*param_i;
      if(param_i==0){
        before_param_mat[out_i] = before_ploss.center;
        after_param_mat[out_i] = after_ploss.center;
      }else{
        before_param_mat[out_i] = before_ploss.spread;
        after_param_mat[out_i] = after_ploss.spread;
      }
      param_i++;
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
 int *seg_end, int *depth_vec,
 double *subtrain_loss, double *validation_loss,
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
  for(int data_i=0; data_i<n_data; data_i++){
    if(!std::isfinite(data_vec[data_i])){
      return ERROR_DATA_MUST_BE_FINITE;
    }
  }
  Distribution* dist_ptr = dist_umap.at(distribution_str);
  for(int data_i=0; data_i<n_data; data_i++){
    int status = dist_ptr->check_data(data_vec[data_i]);
    if(status != 0){
      return status;
    }
  }
  OutArrays out_arrays
    (dist_ptr, max_segments,
     seg_end, depth_vec, subtrain_loss, validation_loss,
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
  ParamsLoss full_ploss, missing_ploss;
  dist_ptr->estimate_params(&full_ploss, V.subtrain, 0, n_subtrain-1),
  out_arrays.save
    (0,
     full_ploss.loss,
     dist_ptr->loss_for_params(V.validation, full_ploss, 0, n_subtrain-1),
     n_subtrain-1,
     0,
     full_ploss,
     missing_ploss,
     -2, -2, n_subtrain, -2);
  // Add a segment and split to the set of candidates.
  V.maybe_add
    (0, n_subtrain-1, 0, 0, subtrain_loss[0], validation_loss[0], 0);
  // initialize to infinite cost and missing values, which is
  // necessary when number of segments returned is less than
  // max_segments: infinite cost rows are removed from the resulting
  // splits table in the R code.
  for(int seg_i=1; seg_i < max_segments; seg_i++){
    out_arrays.save
      (seg_i, INFINITY, INFINITY, -2, -2,
       missing_ploss, missing_ploss, -2, -2, -2, -2);
  }
  // Loop over splits. During each iteration we find the Segment/split
  // which results in the best loss decrease, store the resulting
  // model parameters, and add new Segment/split candidates if
  // necessary.
  int seg_i = 0;
  while(V.container_ptr->not_empty() && ++seg_i < max_segments){
    // Store loss and model parameters associated with this split.
    const Segment seg = V.container_ptr->get_best();
    out_arrays.save
      (seg_i, 
       subtrain_loss[seg_i-1] + seg.best_decrease,
       validation_loss[seg_i-1] + seg.validation_decrease,
       seg.best_split.this_end,
       seg.depth,
       seg.best_split.before,
       seg.best_split.after,
       seg.invalidates_index,
       seg.invalidates_after,
       seg.best_split.this_end - seg.first_i + 1,
       seg.last_i - seg.best_split.this_end);
    // Finally add new split candidates if necessary.
    V.maybe_add
      (seg.first_i, seg.best_split.this_end,
       0,//invalidates_after=0 => before_mean invalidated.
       seg_i, seg.best_split.before.loss,
       seg.before_validation_loss,
       seg.depth);
    V.maybe_add
      (seg.best_split.this_end+1, seg.last_i,
       1,//invalidates_after=1 => after_mean invalidated.
       seg_i, seg.best_split.after.loss,
       seg.after_validation_loss,
       seg.depth);
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
  return &(dist_umap.at(distribution_str)->param_names_vec);
}
