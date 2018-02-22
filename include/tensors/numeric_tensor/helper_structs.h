#ifndef _CPSFIT_NUMERIC_TENSOR_HELPER_STRUCTS_H_
#define _CPSFIT_NUMERIC_TENSOR_HELPER_STRUCTS_H_

#include<iostream>

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

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

CPSFIT_END_NAMESPACE
#endif
