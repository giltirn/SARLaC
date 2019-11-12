#ifndef _GENERIC_ET_BINOP_SV_H_
#define _GENERIC_ET_BINOP_SV_H_

#include<utils/macros.h>
#include<ET/generic_ET/tagging.h>
#include<ET/generic_ET/getelem.h>
#include<ET/generic_ET/value_storage.h>

//Binary operations with scalar-vector form

CPSFIT_START_NAMESPACE

template<template<typename,typename> class Op, typename T, typename U>
struct binaryScalarLeftHelper{
  typedef typename rvalueStoreType<T>::type rT;
  typedef typename rvalueStoreType<U>::type rU;
  
  static inline Op<ETscalarEval<T>,rU> doit(const T &a, U &&b){
    return Op<ETscalarEval<T>,rU>(a,std::move(b));
  }
  static inline Op<ETscalarEval<T>,ETeval<U> > doit(const T &a, const U &b){
    return Op<ETscalarEval<T>,ETeval<U> >(a,b);
  }
};

#define ET_BINOP_SCALAR_LEFT(NAME, OP, OPNAME)		 \
template<typename Leaf1, typename Leaf2> \
struct NAME{ \
  typedef ETleafTag ET_leaf_mark; \
  Leaf1 a; \
  Leaf2 b; \
  typedef ENABLE_IF_TWO_ET_LEAF(Leaf1,Leaf2,typename Leaf2::ET_tag) ET_tag; \
  \
  NAME(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){} \
  \
  inline decltype(a.value() OP b[0]) operator[](const int i) const{ return a.value() OP b[i]; } \
  inline decltype(b.common_properties()) common_properties() const{ return b.common_properties(); } \
};   \
template<typename T,typename U, typename std::enable_if<is_scalar_value<typename std::decay<T>::type>::value && has_ET_tag<typename std::decay<U>::type>::value , int>::type = 0> \
inline auto OPNAME(T &&a, U &&b)->decltype(binaryScalarLeftHelper<NAME,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T >(a),std::forward<U>(b))){ \
  return binaryScalarLeftHelper<NAME,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)); \
}

ET_BINOP_SCALAR_LEFT(ETscalarPlusLeft, +, operator+);
ET_BINOP_SCALAR_LEFT(ETscalarMinusLeft, -, operator-);
ET_BINOP_SCALAR_LEFT(ETscalarTimesLeft, *, operator*);

template<typename T,int has_ET_tag>
struct _get_ET_elem_type{};

template<typename T>
struct _get_ET_elem_type<T,1>{
  typedef typename std::decay<T>::type Tbare;
  typedef typename std::decay<decltype(  getElem<Tbare>::elem(*((Tbare const*)NULL),0) )>::type type;
};
template<typename T>
struct _get_ET_elem_type<T,0>{
  typedef empty_t type;
};

#define BASE(T) typename std::decay<T>::type
#define GET_ELEMENT_TYPE(T) typename _get_ET_elem_type<T, has_ET_tag<BASE(T)>::value>::type


#define ET_BINOP_SCALAR_LEFT_ELEM(NAME, OP, OPNAME)		 \
template<typename Leaf1, typename Leaf2> \
struct NAME{ \
  typedef ETleafTag ET_leaf_mark; \
  Leaf1 a; \
  Leaf2 b; \
  typedef ENABLE_IF_TWO_ET_LEAF(Leaf1,Leaf2,typename Leaf2::ET_tag) ET_tag; \
  \
  NAME(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){} \
  \
  inline decltype(a.value() OP b[0]) operator[](const int i) const{ return a.value() OP b[i]; } \
  inline decltype(b.common_properties()) common_properties() const{ return b.common_properties(); } \
};   \
 template<typename T,typename U, typename std::enable_if<has_ET_tag<BASE(U)>::value && std::is_same<GET_ELEMENT_TYPE(U),BASE(T)>::value && !is_scalar_value<BASE(T)>::value, int>::type = 0> \
   inline auto OPNAME(T &&a, U &&b)->decltype(binaryScalarLeftHelper<NAME,BASE(T),BASE(U)>::doit(std::forward<T>(a),std::forward<U>(b))){ \
   return binaryScalarLeftHelper<NAME,BASE(T), BASE(U)>::doit(std::forward<T>(a),std::forward<U>(b)); \
}

ET_BINOP_SCALAR_LEFT_ELEM(ETscalarTimesLeftElem, *, operator*);


CPSFIT_END_NAMESPACE

#endif
