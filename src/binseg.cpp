#include "binseg.h"
#include <math.h>//INFINITY
#include <stdexcept>      // std::out_of_range

double Cumsum::get_sum(int first, int last){
  double total = cumsum_vec[last];
  if(0 < first){
    total -= cumsum_vec[first-1];
  }
  return total;
}

double Set::get_mean(int first, int last){
  total_weighted_data = weighted_data.get_sum(first, last);
  total_weights = weights.get_sum(first, last);
  return total_weighted_data/total_weights;
}
void Set::set_mean_loss(int first, int last, double *mean, double *loss){
  *mean = get_mean(first, last);
  *loss = get_loss(first, last, *mean);
}
void Set::set_mean_loss(int first, int last, MeanLoss *ML){
  set_mean_loss(first, last, &(ML->mean), &(ML->loss));
}
double Set::get_loss(int first, int last, double subtrain_mean){
  total_weighted_data = weighted_data.get_sum(first, last);
  total_weights = weights.get_sum(first, last);
  return instance_loss(total_weights, total_weighted_data, subtrain_mean);
}
void Set::resize_cumsums(int vec_size){
  weights.cumsum_vec.resize(vec_size);
  weighted_data.cumsum_vec.resize(vec_size);
}
void Set::write_cumsums(int write_index){
  weights.cumsum_vec[write_index] = total_weights;
  weighted_data.cumsum_vec[write_index] = total_weighted_data;
}

static dist_map_type dist_map;
Distribution::Distribution
(const char *name, compute_fun compute, update_fun update){
  compute_loss = compute;
  update_loss = update;
  dist_map.emplace(name, this);
}
// we need get_dist_map in this file to query map contents from
// interface.cpp
dist_map_type* get_dist_map(void){
  return &dist_map;
}

#define FUN(x,y) x##y
#define DISTRIBUTION(NAME, COMPUTE, UPDATE) \
  double FUN(NAME,compute) (double N, double sum, double mean){ \
  return COMPUTE;\
  } \
  void FUN(NAME,update) (double *loss, Set &set){ \
  UPDATE; \
  } \
  static Distribution NAME( #NAME, FUN(NAME,compute), FUN(NAME,update) ); 

DISTRIBUTION(mean_norm,
             mean*(N*mean-2*sum), // square loss minus constant term.
             *loss += set.total_weighted_squares) //add back constant.
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
             mean*N - log(mean)*sum, // neg log lik minus constant term.
             ) // dont add constant term to loss.
/* poisson likelihood:

prob_i = m^{x_i} * exp(-m) / (x_i !)

log(prob_i) = x_i*log(m) - m - x_i!

-log(prob_i) = w_i*x_i*log(M) - w_i*M - w_i*(x_i!)
   
poisson loss with weights:

     sum_i w_i [M - x_i log(M)] = M*(sum_i w_i) - log(M)*(sum_i x_i w_i)

  = M*W - log(M)*S.
  
 */

// Split class stores info for a single candidate split to consider.
class Split {
public:
  int this_end;//index of last data point on the first/before segment.
  MeanLoss before, after;
  double set_mean_loss
  (Set &subtrain, int first, int end_i, int last){
    this_end = end_i;
    subtrain.set_mean_loss(first, end_i, &before);
    subtrain.set_mean_loss(end_i+1, last, &after);
    return before.loss + after.loss;
  }
};

// Segment class stores information about a segment that could be
// split.
class Segment {
public:
  //indices of first/last data points on this segment.
  int first_i, last_i;
  //index and after indicator of a previous segment whose mean
  //parameter is invalidated if this segment in included. For example
  //invalidates_index=2 and invalidates_after=0 implies before_mean
  //parameter at index 2 is invalidated; invalidates_index=1 and
  //invalidates_after=1 implies after_mean parameter at index 1 is
  //invalidated. It is necessary to store this information in order to
  //recover the optimal mean parameters of any model size.
  int invalidates_index, invalidates_after;
  double best_decrease, validation_decrease;
  double before_validation_loss, after_validation_loss;
  Split best_split;
  int n_changes() const {
    return last_i-first_i;
  }
  // Segments are kept sorted by best_decrease value, so that we can
  // find the best segment to split in constant O(1) time.
  friend bool operator<(const Segment& l, const Segment& r)
  {
    if(l.best_decrease == r.best_decrease){
      // if two segments are equally good to split in terms of the
      // loss, then to save time we should split the larger.
      return l.n_changes() > r.n_changes();
    }else{
      return l.best_decrease < r.best_decrease;
    }
  }
  // constructor which considers all possible splits from first to
  // last, and then stores the best split and loss decrease.
  Segment
  (Set &subtrain, Set &validation,
   int first_data, int last_data,
   int first_candidate, int last_candidate,
   int invalidates_after, int invalidates_index,
   double loss_no_split, double validation_loss_no_split
   ): first_i(first_data), last_i(last_data),
      invalidates_index(invalidates_index),
      invalidates_after(invalidates_after)
  {
    Split candidate_split;
    double best_loss_split = INFINITY, loss_split;
    // for loop over all possible splits on this Segment.
    for(int candidate=first_candidate; candidate<=last_candidate; candidate++){
      loss_split = candidate_split.set_mean_loss
	(subtrain, first_data, candidate, last_data);
      if(loss_split < best_loss_split){
	best_loss_split = loss_split;
	best_split = candidate_split;
      }
    }
    best_decrease = best_loss_split - loss_no_split;
    // TODO combine this logic with Split class?
    before_validation_loss = validation.get_loss
      (first_data, best_split.this_end, best_split.before.mean);
    after_validation_loss = validation.get_loss
      (best_split.this_end+1, last_data, best_split.after.mean);
    double validation_loss_split =
      before_validation_loss + after_validation_loss;
    validation_decrease =
      validation_loss_split - validation_loss_no_split;
  }
};

typedef std::multiset<Segment> segment_set_type;

class Candidates {
public:
  segment_set_type segment_set;
  Set subtrain, validation;
  int min_segment_length;
  // computes the cumulative sum vectors in linear O(n_data) time.
  int init
  (const double *data_vec, const double *weight_vec, const int n_data,
   const double *position_vec, const int *is_validation_vec,
   double *subtrain_borders, compute_fun distribution_loss,
   int min_segment_length_arg
   ){
    min_segment_length = min_segment_length_arg;
    subtrain.instance_loss = distribution_loss;
    validation.instance_loss = distribution_loss;
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
          this_set->total_weighted_data += data_value * weight_value;
          this_set->total_weighted_squares +=
            data_value * data_value * weight_value;
          this_set->total_weights += weight_value;
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
      // if it is possible to split (more than one data point on this
      // segment) then insert new segment into the candidates set.
      segment_set.emplace
	(subtrain, validation, 
	 first_data, last_data,
         first_candidate, last_candidate,
	 invalidates_after, invalidates_index,
	 loss_no_split, validation_loss_no_split);
    }
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
   - before_mean: mean before split.
   - after_mean: mean after split.
   - before_size: number of data before this split.
   - after_size: number of data after this split.
   - invalidates_index: index of mean invalidated by this split.
   - invalidates_after: indicates if after mean invalidated by this split.

   See coef method in R code for a procedure that uses these output
   arrays to efficiently compute the segment means for any model size.
 */
int binseg
(const double *data_vec, const double *weight_vec,
 const int n_data, const int max_segments, const int min_segment_length,
 const int *is_validation_vec, const double *position_vec,
 const char *distribution_str,
 double *subtrain_borders, 
 int *seg_end, double *subtrain_loss, double *validation_loss,
 double *before_mean, double *after_mean, 
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
  Distribution* dist_ptr;
  try{
    dist_ptr = dist_map.at(distribution_str);
  }
  catch(const std::out_of_range& err){
    return ERROR_UNRECOGNIZED_DISTRIBUTION;
  }
  Candidates V;
  // Begin by initializing cumulative sum vectors.
  int n_subtrain = V.init
    (data_vec, weight_vec, n_data, position_vec, is_validation_vec,
     subtrain_borders, dist_ptr->compute_loss, min_segment_length);
  if(n_subtrain < max_segments*min_segment_length){
    return ERROR_TOO_MANY_SEGMENTS;
  }
  // Then store the trivial segment mean/loss (which starts at the
  // first and ends at the last data point).
  V.subtrain.set_mean_loss
    (0, n_subtrain-1, before_mean, subtrain_loss);
  validation_loss[0] = V.validation.get_loss
    (0, n_subtrain-1, *before_mean);
  before_size[0] = n_subtrain;
  seg_end[0] = n_subtrain-1;
  after_mean[0] = INFINITY;
  after_size[0] = -2; // unused/missing indicator.
  invalidates_index[0]=-2;
  invalidates_after[0]=-2;
  // Add a segment and split to the set of candidates.
  V.maybe_add(0, n_subtrain-1, 0, 0, subtrain_loss[0], validation_loss[0]);
  // initialize to infinite cost and missing values, which is
  // necessary when number of segments returned is less than
  // max_segments: infinite cost rows are removed from the resulting
  // splits table in the R code.
  for(int seg_i=1; seg_i < max_segments; seg_i++){
    subtrain_loss[seg_i] = INFINITY;
    validation_loss[seg_i] = INFINITY;
    seg_end[seg_i] = -2;
    before_mean[seg_i] = INFINITY;
    after_mean[seg_i] = INFINITY;
    invalidates_index[seg_i] = -2;
    invalidates_after[seg_i] = -2;
    before_size[seg_i] = -2;
    after_size[seg_i]  = -2;
  }
  // Loop over splits. During each iteration we find the Segment/split
  // which results in the best loss decrease, store the resulting
  // model parameters, and add new Segment/split candidates if
  // necessary.
  segment_set_type::iterator it;
  // The multiset is a red-black tree which is sorted in increasing
  // order (by best_decrease values), so the first element is always
  // the segment which results in the best loss decrease.
  int seg_i = 0;
  while
    ((it = V.segment_set.begin()) != V.segment_set.end() &&
     ++seg_i < max_segments
     ){
    // Store loss and model parameters associated with this split.
    subtrain_loss[seg_i] = subtrain_loss[seg_i-1] + it->best_decrease;
    validation_loss[seg_i] =
      validation_loss[seg_i-1] + it->validation_decrease;
    seg_end[seg_i] = it->best_split.this_end;
    before_mean[seg_i] = it->best_split.before.mean;
    after_mean[seg_i] = it->best_split.after.mean;
    // Also store invalidates index/after so we know which previous
    // model parameters are no longer used because of this split.
    invalidates_index[seg_i] = it->invalidates_index;
    invalidates_after[seg_i] = it->invalidates_after;
    // The sizes below are not strictly necessary to store (they can
    // be derived from start/end) but it makes it easier to analyze
    // the time complexity of the algorithm. Splits which result in
    // equal sizes before/after make the algorithm run fastest: first
    // split sizes N/2,N/2, second split sizes N/4,N/4, etc. Worst
    // case is when first split sizes 1,N-1 second split sizes 1,N-2,
    // etc.
    before_size[seg_i] = it->best_split.this_end - it->first_i + 1;
    after_size[seg_i]  = it->last_i - it->best_split.this_end;
    // Finally add new split candidates if necessary.
    V.maybe_add
      (it->first_i, it->best_split.this_end,
       0,//invalidates_after=0 => before_mean invalidated.
       seg_i, it->best_split.before.loss,
       it->before_validation_loss);
    V.maybe_add
      (it->best_split.this_end+1, it->last_i,
       1,//invalidates_after=1 => after_mean invalidated.
       seg_i, it->best_split.after.loss,
       it->after_validation_loss);
    // Erase at end because we need it->values during maybe_add
    // inserts above.
    V.segment_set.erase(it);
  }
  for(int seg_i=0; seg_i < max_segments; seg_i++){
    dist_ptr->update_loss(subtrain_loss+seg_i, V.subtrain);
    dist_ptr->update_loss(validation_loss+seg_i, V.validation);
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
  
