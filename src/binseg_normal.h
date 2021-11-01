#define ERROR_TOO_MANY_SEGMENTS 2

int binseg_normal
(const double *data_vec, const int n_data, const int max_segments,
 const int *is_validation_vec, const int *position_vec,
 int *seg_end, double *loss, double *validation_loss,
 double *before_mean, double *after_mean,
 int *, int *,
 int *invalidates_index, int *invalidates_before);

