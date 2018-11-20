#ifndef _GENERIC_ET_UNARY_H_
#define _GENERIC_ET_UNARY_H_

#include<utils/macros.h>
#include<ET/generic_ET/tagging.h>
#include<ET/generic_ET/getelem.h>
#include<ET/generic_ET/value_storage.h>

//Unary operations

CPSFIT_START_NAMESPACE

//Make sure the compiler is able to find the defaults
using ::sqrt;
using ::exp;
using ::log;

template<template<typename> class Op, typename T>
struct unaryHelper{
  typedef typename rvalueStoreType<T>::type rT;
  
  static inline Op<ETeval<T> > doit(const T &a){
    return Op<ETeval<T> >(a);
  }
  static inline Op<rT> doit(T &&a){
    return Op<rT>(std::move(a));
  }
};
#define ET_UNOP(NAME, OP, OPNAME)		 \
template<typename Leaf> \
struct NAME{ \
  typedef ETleafTag ET_leaf_mark; \
  Leaf a; \
  typedef ENABLE_IF_ET_LEAF(Leaf, typename Leaf::ET_tag) ET_tag;	\
  \
  NAME(Leaf &&aa): a(std::move(aa)){} \
  \
  inline decltype(OP(a[0])) operator[](const int i) const{ return OP(a[i]); } \
  inline decltype(a.common_properties()) common_properties() const{ return a.common_properties(); } \
};   \
template<typename T, typename std::enable_if<has_ET_tag<typename std::decay<T>::type>::value , int>::type = 0> \
inline auto OPNAME(T &&a)->decltype(unaryHelper<NAME,typename std::decay<T>::type>::doit(std::forward<T >(a))){ \
  return unaryHelper<NAME,typename std::decay<T>::type>::doit(std::forward<T>(a)); \
}

template<typename T>
inline T negate(const T &a){ return -a; }

ET_UNOP(ETnegate, negate, operator-);
ET_UNOP(ETexp, exp, exp);
ET_UNOP(ETsqrt, sqrt, sqrt);
ET_UNOP(ETlog, log, log);

CPSFIT_END_NAMESPACE

#endif
