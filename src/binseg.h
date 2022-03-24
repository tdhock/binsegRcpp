#include <math.h>//INFINITY
#include <stdexcept>      // std::out_of_range
#include <string>
#include <string.h>
#include <list>
#include <algorithm>
#include <set>//multiset
#include <unordered_map>
#include <vector>
#define ERROR_TOO_MANY_SEGMENTS 2
#define ERROR_UNRECOGNIZED_DISTRIBUTION 3
#define ERROR_MIN_SEGMENT_LENGTH_MUST_BE_POSITIVE 5
#define ERROR_UNRECOGNIZED_CONTAINER 6
#define ERROR_POSITIONS_MUST_INCREASE -4

// This class saves the optimal parameters/loss value for each segment
// (before and after) resulting from a split. Currently there is only
// one parameter (mean) because the model is normal change in mean
// with constant variance but there could be more parameters for other
// models (e.g., normal change in mean and variance).
class MeanLoss {
public:
  double mean, loss;
};

typedef double(*compute_fun)(double,double,double);

// This class computes and stores the statistics that we need to
// compute the optimal loss/parameters of a segment from first to
// last. In the case of normal change in mean with constant variance
// the only statistic we need is the cumulative sum.
class Cumsum {
public:
  compute_fun instance_loss;
  std::vector<double> cumsum_vec;
  double get_sum(int first, int last);
};

class Set {// either subtrain or validation.
public:
  Cumsum weights, weighted_data;
  compute_fun instance_loss;
  double total_weighted_data=0, total_weights=0, total_weighted_squares=0;
  double get_mean(int first, int last);
  void set_mean_loss(int first, int last, double *mean, double *loss);
  void set_mean_loss(int first, int last, MeanLoss *ML);
  double get_loss(int first, int last, double subtrain_mean);
  void resize_cumsums(int vec_size);
  void write_cumsums(int write_index);
};

typedef void (*update_fun)(double*, Set&);

class Distribution {
public:
  compute_fun compute_loss;
  update_fun update_loss;
  Distribution();
  Distribution(const char *name, compute_fun compute, update_fun update);
};

typedef std::unordered_map<std::string, Distribution*> dist_map_type;
dist_map_type* get_dist_map(void);  

int get_n_subtrain(const int, const int*);

int binseg 
(const double *data_vec, const double *weight_vec,
 const int n_data, const int max_segments, const int min_segment_length,
 const int *is_validation_vec, const double *position_vec,
 const char *distribution_str,
 const char *container_str,
 double *pos_end,
 int *seg_end, double *loss, double *validation_loss,
 double *before_mean, double *after_mean,
 int *, int *,
 int *invalidates_index, int *invalidates_before);

// Split class stores info for a single candidate split to consider.
class Split {
public:
  int this_end;//index of last data point on the first/before segment.
  MeanLoss before, after;
  double set_mean_loss(Set &subtrain, int first, int end_i, int last);
};

class Segment {
public:
  int first_i, last_i;
  int invalidates_index, invalidates_after;
  double best_decrease, validation_decrease;
  double before_validation_loss, after_validation_loss;
  Split best_split;
  int n_changes(void) const;
  friend bool operator<(const Segment& l, const Segment& r){
    if(l.best_decrease == r.best_decrease){
      // if two segments are equally good to split in terms of the
      // loss, then to save time we should split the larger.
      return l.n_changes() > r.n_changes();
    }else{
      return l.best_decrease < r.best_decrease;
    }
  }
  Segment
  (Set &subtrain, Set &validation,
   int first_data, int last_data,
   int first_candidate, int last_candidate,
   int invalidates_after, int invalidates_index,
   double loss_no_split, double validation_loss_no_split
   );
};

class Container {
public:
  virtual void insert(Segment&) = 0;
  virtual int get_size(void) = 0;
  virtual const Segment* set_best(void) = 0;
  virtual void remove_best(void) = 0;
  virtual ~Container() {};
  bool not_empty(void){
    return get_size() > 0;
  }
};

typedef Container* (*construct_fun_type)(void);
typedef void (*destruct_fun_type)(Container*);
class ContainerFactory {
public:
  construct_fun_type construct_fun_ptr;
  destruct_fun_type destruct_fun_ptr;
  ContainerFactory(const char *name, construct_fun_type construct, destruct_fun_type destruct);
};
typedef std::unordered_map<std::string, ContainerFactory*> factory_map_type;
factory_map_type* get_factory_map(void);
