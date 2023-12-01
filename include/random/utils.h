#ifndef _RANDOM_UTILS_H_
#define _RANDOM_UTILS_H_

#include <algorithm>

#include<config.h>
#include<utils/macros.h>
#include<random/rng.h>

SARLAC_START_NAMESPACE

//Fisher-Yates shuffle algorithm
std::vector<int> randomPermutation(const std::vector<int> &in, RNGstore &rng = RNG){
  int N = in.size();
  std::vector<int> out(in);
  for(int i=N-1; i>0; i--){
    std::uniform_int_distribution<> dis(0, i);
    int swap_with = dis(rng());
    std::swap(out[i], out[swap_with]);
  }
  return out;
}

template<typename T>
std::vector<T> randomPermutation(const std::vector<T> &in, RNGstore &rng = RNG){
  int N = in.size();
  std::vector<int> idx(N);
  for(int i=0;i<N;i++) idx[i] = i;
  
  idx = randomPermutation(idx, rng);
  
  std::vector<T> out(N);
  for(int i=0;i<N;i++) out[i] = in[idx[i]];
  return out;
}



SARLAC_END_NAMESPACE

#endif
