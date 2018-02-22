#ifndef _CPSFIT_NUMERIC_VECTOR_ET_H_
#define _CPSFIT_NUMERIC_VECTOR_ET_H_

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_vector/class.h>

CPSFIT_START_NAMESPACE

template<typename Tag>
struct is_NumericVector_tag{
  enum {value = 0};
};
template<typename A>
struct is_NumericVector_tag<NumericVector<A> >{
  enum {value = 1};
};

//Disable * and / for NumericVector as the meaning is ambiguous
template<typename Numeric>
struct disableGenericETbinOp<ETtimes, NumericVector<Numeric> >{
  enum {value = 1};
};
template<typename Numeric>
struct disableGenericETbinOp<ETdivide, NumericVector<Numeric> >{
  enum {value = 1};
};

CPSFIT_END_NAMESPACE

#endif
