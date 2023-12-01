#ifndef _GENERIC_ET_BINOP_VV_H_
#define _GENERIC_ET_BINOP_VV_H_

#include<utils/macros.h>
#include<utils/utils.h>
#include<ET/generic_ET/tagging.h>
#include<ET/generic_ET/getelem.h>
#include<ET/generic_ET/value_storage.h>

//Binary operations with both operands vector type

SARLAC_START_NAMESPACE

template<template<typename,typename> class Op, typename T, typename U>
struct binaryHelper{
  typedef typename rvalueStoreType<T>::type rT;
  typedef typename rvalueStoreType<U>::type rU;
  
  static inline Op<rT,rU> doit(T &&a, U &&b){ //could be rvalue of operand or of leaf
    return Op<rT,rU>(std::move(a),std::move(b));
  }
  static inline Op<ETeval<T>, ETeval<U> > doit(const T &a, const U &b){
    return Op<ETeval<T>, ETeval<U> >(a,b);
  }
  static inline Op<rT, ETeval<U> > doit(T &&a, const U &b){
    return Op<rT, ETeval<U> >(std::move(a),b);
  }
  static inline Op<ETeval<T>, rU> doit(const T &a, U &&b){
    return Op<ETeval<T>, rU>(a,std::move(b));
  }
};

template<template<typename,typename> class Op, typename Tag>
struct disableGenericETbinOp{
  enum {value = 0};
};

#define ET_BINOP(NAME, OP, OPNAME)		 \
template<typename Leaf1, typename Leaf2> \
struct NAME{ \
  typedef ETleafTag ET_leaf_mark; \
  Leaf1 a; \
  Leaf2 b; \
  typedef ENABLE_IF_TWO_ET_LEAF_EQUAL_TAG(Leaf1,Leaf2, typename Leaf1::ET_tag) ET_tag; \
  \
  NAME(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){ \
    if(a.common_properties() != b.common_properties()) error_exit(std::cout << "Error: binop operands have different common_properties: " << a.common_properties() << ":" << b.common_properties() << std::endl); \
  } \
  \
  inline decltype(a[0] OP b[0]) operator[](const int i) const{ return a[i] OP b[i]; } \
  inline decltype(a.common_properties()) common_properties() const{ return a.common_properties(); } \
};						 \
template<typename T,typename U, \
         typename std::enable_if< \
				  has_ET_tag<typename std::decay<T>::type>::value && \
				  has_ET_tag<typename std::decay<U>::type>::value && \
				  !disableGenericETbinOp<NAME, typename std::decay<T>::type::ET_tag>::value \
				    , int>::type = 0>			\
inline auto OPNAME(T &&a, U &&b)->decltype( binaryHelper<NAME,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)) ) { \
  return binaryHelper<NAME,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)); \
}

ET_BINOP(ETplus, +, operator+);
ET_BINOP(ETminus, -, operator-);
ET_BINOP(ETtimes, *, operator*);
ET_BINOP(ETdivide, /, operator/);

SARLAC_END_NAMESPACE

#endif
