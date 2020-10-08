#include "binseg_normal.h"
#include <math.h>//INFINITY
#include <set>//multiset
#include <vector>
#include <array>

class MeanCost {
public:
  double mean, cost;
};

class Cumsums {
public:
  std::vector<double> cumsum_vec;
  void init(const double *data_vec, const int n_data){
    cumsum_vec.resize(n_data);
    double total = 0.0;
    for(int data_i=0; data_i<n_data; data_i++){
      total += data_vec[data_i];
      cumsum_vec[data_i] = total;
    }
  }
  double get_sum(int first, int last){
    double total = cumsum_vec[last];
    if(0 < first){
      total -= cumsum_vec[first-1];
    }
    return total;
  }
  void first_last_mean_cost
  (int first, int last, double *mean, double *cost){
    double s = get_sum(first, last);
    int N = last-first+1;
    *cost = -s*s/N;
    *mean = s/N;
  }
  void first_last_mean_cost(int first, int last, MeanCost *mc){
    first_last_mean_cost(first, last, &(mc->mean), &(mc->cost));
  }
};
/* Above we compute the optimal square loss for a segment with sum of
   data = s and number of data = N.

   If x_i is data point i, and s = sum_{i=1}^N x_i is the sum over N
   data points, then the optimal mean is s/n and the optimal square
   loss is

   sum_{i=1}^N (x_i - s/N)^2 =

   sum [x_i^2] - 2*(s/N)*s + N*(s/N)^2 =

   sum [x_i^2] - s^2 / N

   The first term (sum of squares of data) can be ignored during
   optimization, because it does not depend on the optimization
   variable (segment mean). It is added back after optimization, in
   the R code.
 */

class Split {
public:
  int this_end;
  std::array<MeanCost, 2> segs;
  double set_mean_cost
  (Cumsums &cumsums, int first, int end_i, int last){
    this_end = end_i;
    cumsums.first_last_mean_cost(first, end_i, &segs[0]);
    cumsums.first_last_mean_cost(end_i+1, last, &segs[1]);
    return segs[0].cost + segs[1].cost;
  }
};

class Segment {
public:
  int first, last, invalidates_after, invalidates_index;
  double best_decrease;
  Split best_split;
  friend bool operator<(const Segment& l, const Segment& r)
  {
    return l.best_decrease < r.best_decrease;
  }  
  Segment
  (Cumsums &cumsums, int first, int last,
   int invalidates_after, int invalidates_index,
   double cost_no_split
   ): first(first), last(last),
      invalidates_after(invalidates_after),
      invalidates_index(invalidates_index)
  {
    Split candidate_split;
    int n_candidates = last-first;
    double best_cost_split = INFINITY, cost_split;
    for(int ci=0; ci<n_candidates; ci++){
      cost_split = candidate_split.set_mean_cost
	(cumsums, first, first+ci, last);
      if(cost_split < best_cost_split){
	best_cost_split = cost_split;
	best_split = candidate_split;
      }
    }
    best_decrease = best_cost_split - cost_no_split;
  }
};

class Candidates {
public:
  std::multiset<Segment> candidates;
  Cumsums cumsums;
  void maybe_add
  (int first, int last,
   int invalidates_after, int invalidates_index,
   double cost_no_split)
  {
    if(first < last){
      // if it is possible to split (more than one data point on this
      // segment) then insert into the tree with key=best decrease in
      // cost, value=index of this segment for retrieval of parameters
      // / inserting new segments / deleting this segment.
      candidates.emplace
	(cumsums, first, last,
	 invalidates_after, invalidates_index, cost_no_split);
    }
  }
};

/* Binary segmentation algorithm for change in mean in the normal
   distribution, cost function is square loss.

   This code assumes, and the code which calls this function needs to
   have error checking for, the following:

   At least two data points (1 < n_data), data_vec is an array
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
 int *seg_end, double *cost,
 double *before_mean, double *after_mean, 
 int *before_size, int *after_size, 
 int *invalidates_index, int *invalidates_after){
  if(n_data < max_segments){
    return ERROR_TOO_MANY_SEGMENTS;
  }
  Candidates V;
  V.cumsums.init(data_vec, n_data);
  seg_end[0] = n_data-1;
  after_mean[0] = INFINITY;
  after_size[0] = -2; // unused/missing indicator.
  invalidates_index[0]=-2;
  invalidates_after[0]=-2;
  V.cumsums.first_last_mean_cost
    (0, n_data-1, before_mean, cost);
  before_size[0] = n_data;
  V.maybe_add(0, n_data-1, 0, 0, cost[0]);
  for(int seg_i=1; seg_i < max_segments; seg_i++){
    // The multiset is a red-black tree which is sorted in increasing
    // order (by best_decrease values), so the first element is always
    // the segment which results in the best cost decrease.
    std::multiset<Segment>::iterator it = V.candidates.begin();
    const Segment* s = &(*(it));
    seg_end[seg_i] = s->best_split.this_end;
    cost[seg_i] = cost[seg_i-1] + s->best_decrease;
    before_mean[seg_i] = s->best_split.segs[0].mean;
    after_mean[seg_i] = s->best_split.segs[1].mean;
    invalidates_index[seg_i] = s->invalidates_index;
    invalidates_after[seg_i] = s->invalidates_after;
    before_size[seg_i] = s->best_split.this_end - s->first + 1;
    after_size[seg_i]  = s->last - s->best_split.this_end;
    // it is not necessary to store sizes (they can be derived from
    // start/end) but it makes it easier to analyze the time
    // complexity of the algorithm. Splits which result in equal sizes
    // before/after make the algorithm run fastest: first split sizes
    // N/2,N/2, second split sizes N/4,N/4, etc. Worst case is when
    // first split sizes 1,N-1 second split sizes 1,N-2, etc.
    V.maybe_add(s->first, s->best_split.this_end,
		0, seg_i, s->best_split.segs[0].cost);
    V.maybe_add(s->best_split.this_end+1, s->last,
		1, seg_i, s->best_split.segs[1].cost);
    V.candidates.erase(it);
  }
  return 0;//SUCCESS.
}

