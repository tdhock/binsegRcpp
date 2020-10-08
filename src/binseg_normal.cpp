#include "binseg_normal.h"
#include <math.h>//INFINITY
#include <set>//multiset

/* Compute optimal square loss for a segment with sum of data = s and
   number of data = N.

   If x_i is data point i, and s = sum_{i=1}^N x_i is the sum over N
   data points, then the optimal mean is s/n and the optimal square
   loss is

   sum_{i=1}^N (x_i - s/N)^2 =

   sum [x_i^2 - 2*(s/N)*s + N*(s/N)^2] =

   sum [x_i^2] - s^2 / N

   The first term (sum of squares of data) can be ignored during
   optimization, because it does not depend on the optimization
   variable (segment mean). It is added back after optimization, in
   the R code.
 */
double optimal_cost(double s, int N){
  return -s*s/N;
}

class Segment {
public:
  int first, last, best_end, invalidates_after, invalidates_index;
  double cost, mean;
  double best_mean_before, best_mean_after;
  double best_cost, best_decrease;
  friend bool operator<(const Segment& l, const Segment& r)
  {
    return l.best_decrease < r.best_decrease;
  }  
  Segment
  (const double *data_vec, int first, int last,
   int invalidates_after, int invalidates_index
   ): first(first), last(last),
      invalidates_after(invalidates_after),
      invalidates_index(invalidates_index)
  {
    double sum_before = 0.0, sum_after = 0.0;
    double cost_before, cost_after, cost_split, data_value;
    for(int i=first; i <= last; i++){
      sum_after += data_vec[i];
    }
    int n_total = last-first+1;
    mean = sum_after/n_total;
    cost = optimal_cost(sum_after, n_total);
    int n_before, n_after, end_i;
    int n_candidates = last-first;
    best_cost = INFINITY;
    for(int ci=0; ci<n_candidates; ci++){
      end_i = first+ci;
      data_value = data_vec[end_i];
      sum_before += data_value;
      sum_after  -= data_value;
      n_before = ci+1;
      n_after = n_candidates-ci;
      cost_before = optimal_cost(sum_before, n_before);
      cost_after  = optimal_cost(sum_after, n_after);
      cost_split = cost_before + cost_after;
      if(cost_split < best_cost){
	best_cost = cost_split;
	best_mean_before = sum_before/n_before;
	best_mean_after = sum_after/n_after;
	best_end = end_i;
      }
    }
    best_decrease = best_cost - cost;
  }
};

class Candidates {
public:
  std::multiset<Segment> candidates;
  const double *data_vec;
  Candidates(const double *data_vec): data_vec(data_vec) {};
  void maybe_add
  (int first, int last, int invalidates_after, int invalidates_index){
    if(first < last){
      // if it is possible to split (more than one data point on this
      // segment) then insert into the tree with key=best decrease in
      // cost, value=index of this segment for retrieval of parameters
      // / inserting new segments / deleting this segment.
      candidates.emplace
	(data_vec, first, last, invalidates_after, invalidates_index);
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
  Candidates V(data_vec);
  V.maybe_add(0, n_data-1, 0, 0);
  const Segment *first = &(*(V.candidates.begin()));
  seg_end[0] = n_data-1;
  cost[0] = first->cost;
  before_mean[0] = first->mean;
  after_mean[0] = INFINITY;
  before_size[0] = n_data;
  after_size[0] = -2; // unused/missing indicator.
  invalidates_index[0]=-2;
  invalidates_after[0]=-2;
  for(int seg_i=1; seg_i < max_segments; seg_i++){
    // The multiset is a red-black tree which is sorted in increasing
    // order (by best_decrease values), so the first element is always
    // the segment which results in the best cost decrease.
    std::multiset<Segment>::iterator it = V.candidates.begin();
    const Segment* s = &(*(it));
    seg_end[seg_i] = s->best_end;
    cost[seg_i] = cost[seg_i-1] + s->best_decrease;
    before_mean[seg_i] = s->best_mean_before;
    after_mean[seg_i] = s->best_mean_after;
    invalidates_index[seg_i]=s->invalidates_index;
    invalidates_after[seg_i]=s->invalidates_after;
    before_size[seg_i]=s->best_end - s->first + 1;
    after_size[seg_i]=s->last - s->best_end;
    // it is not necessary to store sizes (they can be derived from
    // start/end) but it makes it easier to analyze the time
    // complexity of the algorithm. Splits which result in equal sizes
    // before/after make the algorithm run fastest: first split sizes
    // N/2,N/2, second split sizes N/4,N/4, etc. Worst case is when
    // first split sizes 1,N-1 second split sizes 1,N-2, etc.
    V.maybe_add(s->first, s->best_end, 0, seg_i);
    V.maybe_add(s->best_end+1, s->last, 1, seg_i);
    V.candidates.erase(it);
  }
  return 0;
}

