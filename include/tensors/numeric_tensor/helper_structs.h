#ifndef _SARLAC_NUMERIC_TENSOR_HELPER_STRUCTS_H_
#define _SARLAC_NUMERIC_TENSOR_HELPER_STRUCTS_H_

#include<iostream>

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

template<typename T,int Rank,int D>
struct _tensor_helper{
  static inline size_t vol(const int* sizes){
    return sizes[0]*_tensor_helper<T,Rank,D-1>::vol(sizes+1);
  }
  static inline size_t map(const size_t accum, const int *elem, const int *sizes){
    return _tensor_helper<T,Rank,D-1>::map( elem[0] + accum*sizes[0], elem+1,sizes+1);
  }
  static inline void unmap(int *into, const size_t off, const int *sizes, const size_t vol){
    size_t subvol = vol/sizes[0];
    *into = off / subvol;
    _tensor_helper<T,Rank,D-1>::unmap(into+1, off % subvol, sizes+1, subvol);
  }
};
template<typename T,int Rank>
struct _tensor_helper<T,Rank,-1>{
  static inline size_t vol(const int* lst){
    return 1;
  }
  static inline size_t map(const size_t accum, const int *elem, const int *sizes){
    return accum;
  }
  static inline void unmap(int *into, const size_t off, const int *sizes, const size_t vol){};
};

template<typename T, int R>
class NumericTensor;

template<typename T, int R>
struct _NTprinter{
  static void print(std::ostream &os, const NumericTensor<T,R> &t);
};

template<typename T, int N>
struct iterate<NumericTensor<T,N> >{
  static inline int size(const NumericTensor<T,N> &from){ 
    int sz = 1;
    for(int i=0;i<N;i++) sz *= from.size(i);
    return sz;
  }
  static inline std::vector<int> unmap(int i, const NumericTensor<T,N> &from){ 
    std::vector<int> coord(N);
    for(int d=N-1;d>=0;d--){
      coord[d] = i % from.size(d);
      i /= from.size(d);
    }
    return coord;
  }    
  static inline const T& at(const int i, const NumericTensor<T,N> &from){
    std::vector<int> coord = unmap(i,from);
    return from(coord.data());
  }
  static inline T & at(const int i, NumericTensor<T,N> &from){
    std::vector<int> coord = unmap(i,from);
    return from(coord.data());
  }
};  

SARLAC_END_NAMESPACE
#endif
