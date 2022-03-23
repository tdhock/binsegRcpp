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

