#ifndef _CPSFIT_NUMERIC_VECTOR_MATH_H_
#define _CPSFIT_NUMERIC_VECTOR_MATH_H_

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_vector/class.h>

CPSFIT_START_NAMESPACE

//Vector dot product
template<typename T>
T dot(const NumericVector<T> &a, const NumericVector<T> &b){
  assert(a.size() == b.size());
  T out(a(0)); zeroit(out);
  for(int i=0;i<a.size();i++) out = out + a(i) * b(i);
  return out;
}

template<typename T>
T mod2(const NumericVector<T> &m){
  return dot(m,m);
}

CPSFIT_END_NAMESPACE

#endif
