#ifndef _NUMERIC_TENSORS_ET_H_
#define _NUMERIC_TENSORS_ET_H_

#include<template_wizardry.h>
#include<generic_ET.h>
#include<numeric_tensors.h>

template<typename Numeric>
struct disableGenericETbinOp<ETtimes, NumericVector<Numeric> >{
  enum {value = 1};
};
template<typename Numeric>
struct disableGenericETbinOp<ETdivide, NumericVector<Numeric> >{
  enum {value = 1};
};

template<typename A>
struct getElem<NumericSquareMatrix<A> >{
  static inline auto elem(const NumericSquareMatrix<A> &v, const int i)->decltype(v(i/v.size(), i%v.size())){ return v(i/v.size(), i%v.size()); }
  static inline auto elem(NumericSquareMatrix<A> &v, const int i)->decltype(v(i/v.size(), i%v.size())){ return v(i/v.size(), i%v.size()); }
  static inline int common_properties(const NumericSquareMatrix<A> &v){ return v.size(); }
};
template<typename Numeric>
struct disableGenericETbinOp<ETtimes, NumericSquareMatrix<Numeric> >{
  enum {value = 1};
};
template<typename Numeric>
struct disableGenericETbinOp<ETdivide, NumericSquareMatrix<Numeric> >{
  enum {value = 1};
};
template<typename Tag>
struct is_NumericSquareMatrix_tag{
  enum {value = 0};
};
template<typename A>
struct is_NumericSquareMatrix_tag<NumericSquareMatrix<A> >{
  enum {value = 1};
};

template<typename Leaf1, typename Leaf2>
struct ETnumericMatrixMult{
  typedef ETleafTag ET_leaf_mark;
  Leaf1 a;
  Leaf2 b;
  typedef ENABLE_IF_TWO_ET_LEAF_EQUAL_TAG(Leaf1,Leaf2, typename Leaf1::ET_tag) ET_tag;
  
  ETnumericMatrixMult(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){
    assert(aa.common_properties() == bb.common_properties());
  }
  template<typename T>
  static inline auto elem(const T &m, const int i, const int j)->decltype(m[0]){ return m[j + m.common_properties()*i]; }
  
  inline typename std::decay<decltype(a[0])>::type operator[](const int i) const{
    typedef typename std::decay<decltype(a[0])>::type type;
    const int size = a.common_properties();
    const int ai = i/size; const int bi = i%size;
    //type ret = 0.;
    type ret(a[0]);
    zeroit(ret);
    for(int ci=0;ci<size;ci++)
      ret = ret + elem(a,ai,ci) * elem(b,ci,bi);
    return ret;
  }
    
  inline decltype(a.common_properties()) common_properties() const{ return a.common_properties(); }
};
template<typename T,typename U,
         typename std::enable_if<
	   is_NumericSquareMatrix_tag<typename std::decay<T>::type::ET_tag>::value && is_NumericSquareMatrix_tag<typename std::decay<U>::type::ET_tag>::value
				    , int>::type = 0>
inline auto operator*(T &&a, U &&b)->decltype( binaryHelper<ETnumericMatrixMult,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)) ) {
  return binaryHelper<ETnumericMatrixMult,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b));
}

#endif
