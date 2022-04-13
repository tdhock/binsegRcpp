#include <R.h>
#include "best_splits.h"
#include <string>
#include <unordered_map>
#include <cmath> //log2
#define SIZE2SPLITS(SIZE) ( ((SIZE) < min_segment_length*2) ? 0 : (1+(SIZE)-min_segment_length*2) )

Splitter::Splitter
(int n_data_, int min_segment_length_, int n_segments_){
  n_data = n_data_;
  min_segment_length = min_segment_length_;
  n_segments = n_segments_;
}

class Table {
  std::unordered_map< std::string, int > f;
  std::string get_key(int n_data, int changes){
    std::string key = "";
    key += std::to_string(n_data);
    key += ",";
    key += std::to_string(changes);
    return key;
  }
  int get_best(int n_data, int changes){
    return f[get_key(n_data, changes)];
  }
  void put_best(int n_data, int changes, int best){
    f[get_key(n_data, changes)] = best;
  }
};

int Splitter::best_splits(int *out_splits_, int *out_depth_){
  if(min_segment_length < 1){
    return ERROR_BEST_SPLITS_MIN_SEGMENT_LENGTH_MUST_BE_POSITIVE;
  }
  if(n_data < min_segment_length){
    return ERROR_BEST_SPLITS_N_DATA_MUST_BE_AT_LEAST_MIN_SEGMENT_LENGTH;
  }
  int max_segments = n_data / min_segment_length;
  if(max_segments < n_segments){
    return ERROR_BEST_SPLITS_N_SEGMENTS_TOO_LARGE;
  }
  out_splits = out_splits_;
  out_depth = out_depth_;
  out_index=0;
  Table f;
  int n_changes = n_segments-1;
  for(int d_below=0; d_below < n_changes; d_below++){
    int best = 0;
    if(0 < d_below){
      int segs_first_f = d_below+1;
      int segs_second_f = n_changes-segs_first_f;
      int smallest_size = min_segment_length*segs_first_f;
      int largest_size = n_data-segs_second_f*min_segment_length;
      for(int size=smallest_size; size<=largest_size; size++){
      }
    }
    best += SIZE2SPLITS(
  }
  return 0;
}

void Splitter::write_splits_depth(int splits, int depth){
  out_splits[out_index] = splits;
  out_depth[out_index] = depth;
  out_index++;
}

void Splitter::children
(int smaller_size, int larger_size, int depth){
  int smaller_splits = SIZE2SPLITS(smaller_size);
  int larger_splits = SIZE2SPLITS(larger_size);
  write_splits_depth(smaller_splits + larger_splits, depth); 
  split_if_possible(smaller_size, depth);
  split_if_possible(larger_size, depth);
}

void Splitter::split_if_possible
(int size_to_split, int depth){
  int changes_remaining = n_segments-out_index;
  if(0 == changes_remaining)return;
  double half_size = 0.5 * (double)size_to_split;
  int smaller_size;
  if(changes_remaining==1){
    smaller_size = half_size;
  }else{
    int multiple_middle = round(half_size/min_segment_length);
    int multiple = (multiple_middle < changes_remaining) ?
      multiple_middle : changes_remaining;
    int one_size = multiple * min_segment_length;
    smaller_size = (one_size < half_size) ? one_size : size_to_split-one_size;
  }
  int larger_size = size_to_split - smaller_size;
  children(smaller_size, larger_size, depth+1);
}

