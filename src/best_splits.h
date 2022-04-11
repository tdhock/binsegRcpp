#define ERROR_BEST_SPLITS_N_DATA_MUST_BE_AT_LEAST_MIN_SEGMENT_LENGTH 1
#define ERROR_BEST_SPLITS_MIN_SEGMENT_LENGTH_MUST_BE_POSITIVE 2

class Splitter {
public:
  int n_data, min_segment_length, max_segments, out_index;
  int *out_splits, *out_depth;
  Splitter(int n_data_, int min_segment_length_);
  int best_splits(int *out_splits_, int *out_depth_);
  void write_splits_depth(int depth, int splits);
  void split_if_possible(int, int);
  void children(int,int,int);
};
    

