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
#define ERROR_MIN_SEGMENT_LENGTH_MUST_BE_POSITIVE 5
#define ERROR_UNRECOGNIZED_CONTAINER 6
#define ERROR_POSITIONS_MUST_INCREASE -4

// This class saves the optimal parameters/loss value for each segment
// (before and after) resulting from a split. Currently there is only
// one parameter (mean) because the model is normal change in mean
// with constant variance but there could be more parameters for other
// models (e.g., normal change in mean and variance).
class MeanVarLoss {
public:
  double mean, var, loss;
};

typedef std::vector<std::string> param_names_type;
param_names_type* get_param_names(const char*);
class Distribution {
public:
  std::string description;
  param_names_type param_names_vec;
  virtual double compute_loss(double,double,double,double,double,double) = 0;
};

// This class computes and stores a cumsum that we need to compute the
// optimal loss/parameters of a segment from first to last.
class Cumsum {
public:
  Distribution *dist_ptr;
  std::vector<double> cumsum_vec;
  double get_sum(int first, int last);
};

class Set {// either subtrain or validation.
public:
  Cumsum weights, weighted_data, weighted_squares;
  double max_zero_var;
  Distribution *dist_ptr;
  double total_weighted_data=0, total_weights=0, total_weighted_squares=0;
  double get_mean(int first, int last);
  double get_var(int first, int last);
  void set_mean_var_loss(int first, int last, double *mean, double *var, double *loss);
  void set_mean_var_loss(int first, int last, MeanVarLoss*);
  double get_loss(int first, int last, MeanVarLoss&);
  double get_loss(int first, int last, double, double);
  void resize_cumsums(int vec_size);
  void write_cumsums(int write_index);
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
  MeanVarLoss before, after;
  double set_mean_var_loss(Set &subtrain, int first, int end_i, int last);
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
