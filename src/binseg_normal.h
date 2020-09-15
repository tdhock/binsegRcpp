/* -*- compile-command: "R CMD INSTALL .." -*- */

#define ERROR_TOO_MANY_SEGMENTS 2

int binseg_normal
(const double *data_vec, const int n_data, const int max_segments,
 int *seg_end, double *cost,
 double *before_mean, double *after_mean,
 int *, int *,
 int *invalidates_index, int *invalidates_before);

