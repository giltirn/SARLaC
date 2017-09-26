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


struct _NTcommonproperties{
  std::vector<int> sizes;
  _NTcommonproperties(const std::vector<int> &s): sizes(s){}
  _NTcommonproperties(std::initializer_list<int> s): sizes(s){}
  
  inline bool operator==(const _NTcommonproperties &r) const{
    if(sizes.size() != r.sizes.size()) return false;
    for(int i=0;i<sizes.size();i++) if(sizes.at(i) != r.sizes.at(i)) return false;
    return true;
  }
  inline bool operator!=(const _NTcommonproperties &r) const{
    return !(*this == r);
  }
};
std::ostream & operator<<(std::ostream &os, const _NTcommonproperties &p){
  os << "(";
  for(int i=0;i<p.sizes.size()-1;i++) os << p.sizes.at(i) << ",";
  os << p.sizes.back() << ")";
  return os;
}

template<typename T, int R>
struct getElem<NumericTensor<T,R> >{
  static inline auto elem(const NumericTensor<T,R> &v, const int i)->decltype(v.data[i]){ return v.data[i]; }
  static inline auto elem(NumericTensor<T,R> &v, const int i)->decltype(v.data[i]){ return v.data[i]; }
  static inline _NTcommonproperties common_properties(const NumericTensor<T,R> &v){ return _NTcommonproperties(v.dsizes); }
};

template<typename T,int R>
struct disableGenericETbinOp<ETtimes, NumericTensor<T,R> >{
  enum {value = 1};
};
template<typename T,int R>
struct disableGenericETbinOp<ETdivide, NumericTensor<T,R> >{
  enum {value = 1};
};

template<typename Tag>
struct is_Rank2NumericTensor_tag{
  enum {value = 0};
};
template<typename A>
struct is_Rank2NumericTensor_tag<NumericTensor<A,2> >{
  enum {value = 1};
};

template<typename Tag>
struct is_Rank1NumericTensor_tag{
  enum {value = 0};
};
template<typename A>
struct is_Rank1NumericTensor_tag<NumericTensor<A,1> >{
  enum {value = 1};
};

template<typename Leaf1, typename Leaf2>
struct ETrank2numericTensorMult{
  typedef ETleafTag ET_leaf_mark;
  Leaf1 a;
  Leaf2 b;
  typedef ENABLE_IF_TWO_ET_LEAF_EQUAL_TAG(Leaf1,Leaf2, typename Leaf1::ET_tag) ET_tag;

  std::vector<int> sizesA;
  std::vector<int> sizesB;
  
  ETrank2numericTensorMult(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){
    assert(a.common_properties().sizes[1] == b.common_properties().sizes[0]);
    sizesA = a.common_properties().sizes;
    sizesB = b.common_properties().sizes;
  }

  template<typename T>
  static inline auto elem(const T &m, const int i, const int j, const std::vector<int> &msize)->decltype(m[0]){ return m[j + msize[1]*i]; }
  
  inline typename std::decay<decltype(a[0])>::type operator[](const int off) const{
    typedef typename std::decay<decltype(a[0])>::type type;

    //sizesC = {sizesA[0],sizesB[1]}
    //mapping  off = j + sizesB[1]*i
    int i = off / sizesB[1];
    int j = off % sizesB[1];
        
    const int nk = sizesA[1];
    type ret(a[0]);
    zeroit(ret);
    for(int k=0;k<nk;k++)
      ret = ret + elem(a,i,k,sizesA) * elem(b,k,j,sizesB);
    return ret;
  }
    
  inline _NTcommonproperties common_properties() const{ return _NTcommonproperties({sizesA[0],sizesB[1]}); }

};
template<typename T,typename U,
         typename std::enable_if<
	   is_Rank2NumericTensor_tag<typename std::decay<T>::type::ET_tag>::value && is_Rank2NumericTensor_tag<typename std::decay<U>::type::ET_tag>::value
	   , int>::type = 0>
inline auto operator*(T &&a, U &&b)->decltype( binaryHelper<ETrank2numericTensorMult,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)) ) {
  return binaryHelper<ETrank2numericTensorMult,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b));
}




template<typename Leaf1, typename Leaf2>
struct ETNumericTensorMatrixVectorMult{
  typedef ETleafTag ET_leaf_mark;
  Leaf1 a;
  Leaf2 b;
  typedef typename Leaf2::ET_tag ET_tag;//vector type

  std::vector<int> sizesA;
  std::vector<int> sizesB;
  
  ETNumericTensorMatrixVectorMult(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){
    assert(a.common_properties().sizes[1] == b.common_properties().sizes[0]);
    sizesA = a.common_properties().sizes;
    sizesB = b.common_properties().sizes;
  }

  template<typename T>
  static inline auto elem(const T &m, const int i, const int j, const std::vector<int> &msize)->decltype(m[0]){ return m[j + msize[1]*i]; }
  
  inline typename std::decay<decltype(a[0])>::type operator[](const int off) const{
    typedef typename std::decay<decltype(a[0])>::type type;        
    const int nk = sizesA[1];
    type ret(a[0]);
    zeroit(ret);
    for(int k=0;k<nk;k++)
      ret = ret + elem(a,off,k,sizesA) *b[k];
    return ret;
  }
    
  inline _NTcommonproperties common_properties() const{ return _NTcommonproperties({sizesA[0]}); }

};
template<typename T,typename U,
         typename std::enable_if<
	   is_Rank2NumericTensor_tag<typename std::decay<T>::type::ET_tag>::value && is_Rank1NumericTensor_tag<typename std::decay<U>::type::ET_tag>::value
	   , int>::type = 0>
inline auto operator*(T &&a, U &&b)->decltype( binaryHelper<ETNumericTensorMatrixVectorMult,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)) ) {
  return binaryHelper<ETNumericTensorMatrixVectorMult,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b));
}


#endif
