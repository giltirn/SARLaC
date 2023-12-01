#ifndef _GENERIC_ET_BINOP_VS_H_
#define _GENERIC_ET_BINOP_VS_H_

#include<utils/macros.h>
#include<ET/generic_ET/tagging.h>
#include<ET/generic_ET/getelem.h>
#include<ET/generic_ET/value_storage.h>

//Binary operations with vector-scalar form

SARLAC_START_NAMESPACE

template<template<typename,typename> class Op, typename T, typename U>
struct binaryScalarRightHelper{
  typedef typename rvalueStoreType<T>::type rT;
  typedef typename rvalueStoreType<U>::type rU;
  
  static inline Op<rT,ETscalarEval<U> > doit(T &&a, const U &b){
    return Op<rT,ETscalarEval<U> >(std::move(a),b);
  }
  static inline Op<ETeval<T>,ETscalarEval<U> > doit(const T &a, const U &b){
    return Op<ETeval<T>,ETscalarEval<U> >(a,b);
  }
};

#define ET_BINOP_SCALAR_RIGHT(NAME, OP, OPNAME)		 \
template<typename Leaf1, typename Leaf2> \
struct NAME{ \
  typedef ETleafTag ET_leaf_mark; \
  Leaf1 a; \
  Leaf2 b; \
  typedef ENABLE_IF_TWO_ET_LEAF(Leaf1,Leaf2,typename Leaf1::ET_tag) ET_tag; \
  \
  NAME(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){} \
  \
  inline decltype(a[0] OP b.value()) operator[](const int i) const{ return a[i] OP b.value(); } \
  inline decltype(a.common_properties()) common_properties() const{ return a.common_properties(); } \
};   \
template<typename T,typename U, typename std::enable_if<has_ET_tag<typename std::decay<T>::type>::value && is_scalar_value<typename std::decay<U>::type>::value, int>::type = 0> \
inline auto OPNAME(T &&a, U &&b)->decltype(binaryScalarRightHelper<NAME,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T >(a),std::forward<U>(b))){ \
  return binaryScalarRightHelper<NAME,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)); \
}

ET_BINOP_SCALAR_RIGHT(ETscalarPlusRight, +, operator+);
ET_BINOP_SCALAR_RIGHT(ETscalarMinusRight, -, operator-);
ET_BINOP_SCALAR_RIGHT(ETscalarTimesRight, *, operator*);
ET_BINOP_SCALAR_RIGHT(ETscalarDivideRight, /, operator/);

SARLAC_END_NAMESPACE

#endif
