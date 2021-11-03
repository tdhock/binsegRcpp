#include "binseg_normal.h"
#include <math.h>//INFINITY
#include <set>//multiset
#include <vector>
#include <R.h>//Rprintf

// This class saves the optimal parameters/loss value for each segment
// (before and after) resulting from a split. Currently there is only
// one parameter (mean) because the model is normal change in mean
// with constant variance but there could be more parameters for other
// models (e.g., normal change in mean and variance).
class MeanLoss {
public:
  double mean, loss;
};

double square_loss(double N, double sum, double mean){
  return mean*(N*mean-2*sum);
}
/* Above we compute the square loss for a segment with sum of data = s
   and mean parameter m.

   If x_i is data point i, and s = sum_i x_i is the sum over all N
   data points on that segment, then total loss is

   sum_i (m - x_i)^2 = N*m^2 - 2*m*s + sum_i x_i^2

   The last term (sum of squares of data) can be ignored during
   optimization, because it does not depend on the optimization
   variable (segment mean). It is added back after optimization, at
   the end of binseg_normal.
*/

// This class computes and stores the statistics that we need to
// compute the optimal loss/parameters of a segment from first to
// last. In the case of normal change in mean with constant variance
// the only statistic we need is the cumulative sum.
class SetCumsum {
public:
  std::vector<double> cumsum_vec;
  double get_sum(int first, int last){
    double total = cumsum_vec[last];
    if(0 < first){
      total -= cumsum_vec[first-1];
    }
    return total;
  }
  // Compute/store optimal mean and loss values for a single segment
  // from first to last in constant O(1) time.
  void first_last_mean_loss
  (int first, int last, double *mean, double *loss){
    double s = get_sum(first, last);
    int N = last-first+1;
    *mean = s/N;
    *loss = square_loss(N, s, *mean);
  }
  void first_last_mean_loss(int first, int last, MeanLoss *mc){
    first_last_mean_loss(first, last, &(mc->mean), &(mc->loss));
  }
};

double get_validation_loss
(int first, int last, double subtrain_mean,
 SetCumsum &validation, SetCumsum &validation_count
 ){
  double validation_sum = validation.get_sum(first, last);
  double validation_N = validation_count.get_sum(first, last);
  return square_loss(validation_N, validation_sum, subtrain_mean);
}


// Split class stores info for a single candidate split to consider.
class Split {
public:
  int this_end;//index of last data point on the first/before segment.
  MeanLoss before, after;
  double set_mean_loss
  (SetCumsum &cumsums, int first, int end_i, int last){
    this_end = end_i;
    cumsums.first_last_mean_loss(first, end_i, &before);
    cumsums.first_last_mean_loss(end_i+1, last, &after);
    return before.loss + after.loss;
  }
};

// Segment class stores information about a segment that could be
// split.
class Segment {
public:
  //indices of first/last data points on this segment.
  int first, last;
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
    return last-first;
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
  (SetCumsum &subtrain, SetCumsum &validation, SetCumsum &validation_count,
   int first, int last,
   int invalidates_after, int invalidates_index,
   double loss_no_split, double validation_loss_no_split
   ): first(first), last(last),
      invalidates_index(invalidates_index),
      invalidates_after(invalidates_after)
  {
    Split candidate_split;
    int n_candidates = last-first;
    double best_loss_split = INFINITY, loss_split;
    // for loop over all possible splits on this Segment.
    for(int ci=0; ci<n_candidates; ci++){
      loss_split = candidate_split.set_mean_loss
	(subtrain, first, first+ci, last);
      if(loss_split < best_loss_split){
	best_loss_split = loss_split;
	best_split = candidate_split;
      }
    }
    best_decrease = best_loss_split - loss_no_split;
    before_validation_loss = get_validation_loss
      (first, best_split.this_end, best_split.before.mean,
       validation, validation_count);
    after_validation_loss = get_validation_loss
      (best_split.this_end+1, last, best_split.after.mean,
       validation, validation_count);
    double validation_loss_split =
      before_validation_loss + after_validation_loss;
    validation_decrease =
      validation_loss_split - validation_loss_no_split;
  }
};

class Candidates {
public:
  std::multiset<Segment> candidates;
  std::vector<double> change_pos_vec;
  std::vector<int> change_index_vec;
  SetCumsum subtrain, validation, validation_count;
  double subtrain_squares, validation_squares;
  // computes the cumulative sum vectors in linear O(n_data) time.
  int init
  (const double *data_vec, const int n_data,
   const double *position_vec, const int *is_validation_vec
   ){
    int n_validation = 0;
    for(int data_i=0; data_i<n_data; data_i++){
      if(is_validation_vec[data_i]){
	n_validation++;
      }
    }
    int n_subtrain = n_data - n_validation;
    if(n_subtrain == 0)return(n_subtrain);
    subtrain.cumsum_vec.resize(n_subtrain);
    validation.cumsum_vec.resize(n_subtrain);
    validation_count.cumsum_vec.resize(n_subtrain);
    change_pos_vec.resize(n_subtrain+1);
    int last_subtrain_i=-1;
    double pos_total, pos_change, subtrain_total=0, validation_total=0;
    double validation_count_total = 0;
    subtrain_squares=0, validation_squares=0;
    int read_start=0, write_index=0;
    for(int data_i=0; data_i<=n_data; data_i++){
      bool is_subtrain = false;
      bool write_subtrain = false;
      bool write_end = data_i == n_data;
      if(!write_end){
	is_subtrain = !is_validation_vec[data_i];
	write_subtrain = last_subtrain_i >= 0 && is_subtrain;
      }
      //Rprintf("data_i=%d write_subtrain=%d write_end=%d\n", data_i, write_subtrain, write_end);
      if(write_subtrain || write_end){
	if(write_subtrain){
	  pos_total = position_vec[data_i]+position_vec[last_subtrain_i];
	  pos_change = pos_total/2;
	  if(write_index==0){
	    change_pos_vec[write_index] = position_vec[last_subtrain_i]-1;
	    change_index_vec[write_index] = last_subtrain_i;
	  }
	}else{
	  pos_change = position_vec[data_i-1]+1;//last.
	}
	change_pos_vec[write_index+1] = pos_change;
	change_index_vec[write_index+1] = data_i;
	int read_index=read_start;
	while(read_index < n_data && position_vec[read_index] <= pos_change){
	  double data_value = data_vec[read_index];
	  if(is_validation_vec[read_index]){
	    validation_total += data_value;
	    validation_squares += data_value * data_value;
	    validation_count_total += 1;
	  }else{
	    subtrain_total += data_value;
	    subtrain_squares += data_value * data_value;
	  }
	  read_index++;
	}
	//Rprintf("write_index=%d validation_total=%f subtrain_total=%f\n", write_index, validation_total, subtrain_total);
	validation_count.cumsum_vec[write_index] = validation_count_total;
	validation.cumsum_vec[write_index] = validation_total;
	subtrain.cumsum_vec[write_index] = subtrain_total;
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
  (int first, int last,
   int invalidates_after, int invalidates_index,
   double loss_no_split, double validation_loss_no_split
   ){
    if(first < last){
      // if it is possible to split (more than one data point on this
      // segment) then insert new segment into the candidates set.
      candidates.emplace
	(subtrain, validation, validation_count,
	 first, last,
	 invalidates_after, invalidates_index,
	 loss_no_split, validation_loss_no_split);
    }
  }
};

/* Binary segmentation algorithm for change in mean in the normal
   distribution (square loss function).

   This code assumes, and the code which calls this function needs to
   have error checking for, the following:

   At least one data point (0 < n_data), data_vec is an array
   of input data, size n_data.

   Positive number of segments (0 < max_segments), all of the other
   pointers are output arrays of size max_segments (need to be
   allocated by the code which calls this function).

   See coef.binseg_normal in R code for a procedure that uses these
   output arrays to efficiently compute the segment means for any
   model size.
 */
int binseg_normal
(const double *data_vec, const int n_data, const int max_segments,
 const int *is_validation_vec, const double *position_vec, 
 int *seg_end, double *pos_end, double *loss, double *validation_loss,
 double *before_mean, double *after_mean, 
 int *before_size, int *after_size, 
 int *invalidates_index, int *invalidates_after
 ){
  for(int data_i=1; data_i<n_data; data_i++){
    if(position_vec[data_i] <= position_vec[data_i-1]){
      return ERROR_POSITIONS_MUST_INCREASE;
    }
  }
  Candidates V;
  // Begin by initializing cumulative sum vectors.
  int n_subtrain = V.init(data_vec, n_data, position_vec, is_validation_vec);
  if(n_subtrain == 0){
    return ERROR_NEED_AT_LEAST_ONE_SUBTRAIN_DATA;
  }
  if(n_subtrain < max_segments){
    return ERROR_TOO_MANY_SEGMENTS;
  }
  // Then store the trivial segment mean/loss (which starts at the
  // first and ends at the last data point).
  V.subtrain.first_last_mean_loss
    (0, n_subtrain-1, before_mean, loss);
  validation_loss[0] = get_validation_loss
    (0, n_subtrain-1, *before_mean, V.validation, V.validation_count);
  before_size[0] = n_subtrain;
  seg_end[0] = V.change_index_vec[n_subtrain];
  pos_end[0] = V.change_pos_vec[n_subtrain];
  after_mean[0] = INFINITY;
  after_size[0] = -2; // unused/missing indicator.
  invalidates_index[0]=-2;
  invalidates_after[0]=-2;
  // Add a segment and split to the set of candidates.
  V.maybe_add(0, n_subtrain-1, 0, 0, loss[0], validation_loss[0]);
  // For loop over number of segments. During each iteration we find
  // the Segment/split which results in the best loss decrease, store
  // the resulting model parameters, and add new Segment/split
  // candidates if necessary.
  for(int seg_i=1; seg_i < max_segments; seg_i++){
    // The multiset is a red-black tree which is sorted in increasing
    // order (by best_decrease values), so the first element is always
    // the segment which results in the best loss decrease.
    std::multiset<Segment>::iterator it = V.candidates.begin();
    // Store loss and model parameters associated with this split.
    loss[seg_i] = loss[seg_i-1] + it->best_decrease;
    validation_loss[seg_i] = validation_loss[seg_i-1] + it->validation_decrease;
    seg_end[seg_i] = V.change_index_vec[it->best_split.this_end+1];
    pos_end[seg_i] = V.change_pos_vec[it->best_split.this_end+1];
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
    before_size[seg_i] = it->best_split.this_end - it->first + 1;
    after_size[seg_i]  = it->last - it->best_split.this_end;
    // Finally add new split candidates if necessary.
    V.maybe_add
      (it->first, it->best_split.this_end,
       0,//invalidates_after=0 => before_mean invalidated.
       seg_i, it->best_split.before.loss,
       it->before_validation_loss);
    V.maybe_add
      (it->best_split.this_end+1, it->last,
       1,//invalidates_after=1 => after_mean invalidated.
       seg_i, it->best_split.after.loss,
       it->after_validation_loss);
    // Erase at end because we need it->values during maybe_add
    // inserts above.
    V.candidates.erase(it);
  }
  for(int seg_i=0; seg_i < max_segments; seg_i++){
    loss[seg_i] += V.subtrain_squares;
    validation_loss[seg_i] += V.validation_squares;
  }
  return 0;//SUCCESS.
}

