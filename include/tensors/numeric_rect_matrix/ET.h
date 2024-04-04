#pragma once 

#include<tensors/numeric_square_matrix/class.h>
#include<tensors/numeric_square_matrix/ET.h>
#include<tensors/numeric_rect_matrix/class.h>
#include<tensors/numeric_vector.h>

SARLAC_START_NAMESPACE

template<typename A>
struct getElem<NumericRectMatrix<A> >{
  static inline auto elem(const NumericRectMatrix<A> &v, const int i)->decltype(v.linElem(i)){ return v.linElem(i); }
  static inline auto elem(NumericRectMatrix<A> &v, const int i)->decltype(v.linElem(i)){ return v.linElem(i); }
  static inline std::pair<int,int> common_properties(const NumericRectMatrix<A> &v){ return v.size(); }
};
//Disable / for NumericRectMatrix as the meaning is ambiguous. Disable generic * and specialize matrix multiplication
template<typename Numeric>
struct disableGenericETbinOp<ETtimes, NumericRectMatrix<Numeric> >{
  enum {value = 1};
};
template<typename Numeric>
struct disableGenericETbinOp<ETdivide, NumericRectMatrix<Numeric> >{
  enum {value = 1};
};
template<typename Tag>
struct is_NumericRectMatrix_tag{
  enum {value = 0};
};
template<typename A>
struct is_NumericRectMatrix_tag<NumericRectMatrix<A> >{
  enum {value = 1};
};

template<typename Leaf1, typename Leaf2>
struct ETnumericRectMatrixMult{
  typedef ETleafTag ET_leaf_mark;
  Leaf1 a;
  Leaf2 b;
  typedef ENABLE_IF_TWO_ET_LEAF_EQUAL_TAG(Leaf1,Leaf2, typename Leaf1::ET_tag) ET_tag;
  
  ETnumericRectMatrixMult(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){
    assert(a.common_properties().second == b.common_properties().first);
  }
  template<typename T>
  static inline auto elem(const T &m, const int i, const int j)->decltype(m[0]){ return m[j + m.common_properties().second*i]; }
  
  inline typename std::decay<decltype(a[0])>::type operator[](const int i) const{
    typedef typename std::decay<decltype(a[0])>::type type;
    int nj = a.common_properties().second;
    int nk = b.common_properties().second;
    //out size   ni * nk    and indexing  k + nk*i

    int bk = i % nk;
    int ai = i / nk;

    type ret(a[0]);
    zeroit(ret);
    for(int j=0;j<nj;j++)
      ret = ret + elem(a,ai,j) * elem(b,j,bk);
    return ret;
  }
    
  inline decltype(a.common_properties()) common_properties() const{ return {a.common_properties().first, b.common_properties().second}; }
};
template<typename T,typename U,
         typename std::enable_if<
	   is_NumericRectMatrix_tag<typename std::decay<T>::type::ET_tag>::value && is_NumericRectMatrix_tag<typename std::decay<U>::type::ET_tag>::value
				    , int>::type = 0>
inline auto operator*(T &&a, U &&b)->decltype( binaryHelper<ETnumericRectMatrixMult,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)) ) {
  return binaryHelper<ETnumericRectMatrixMult,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b));
}

//Matrix-vector multiplication
template<typename Leaf1, typename Leaf2>
struct ETnumericRectMatrixVectorMult{
  typedef ETleafTag ET_leaf_mark;
  Leaf1 a;
  Leaf2 b;
  typedef typename Leaf2::ET_tag ET_tag; //result is vector
  
  ETnumericRectMatrixVectorMult(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){
    assert(a.common_properties().second == b.common_properties());
  }
  template<typename T>
  static inline auto melem(const T &m, const int i, const int j)->decltype(m[0]){ return m[j + m.common_properties().second*i]; }
  
  inline typename std::decay<decltype(b[0])>::type operator[](const int i) const{
    typedef typename std::decay<decltype(b[0])>::type type;
    const int size = b.common_properties();
    type ret(b[0]);
    zeroit(ret);
    for(int ci=0;ci<size;ci++)
      ret = ret + melem(a,i,ci) * b[ci];
    return ret;
  }
    
  inline decltype(b.common_properties()) common_properties() const{ return a.common_properties().first; }
};
template<typename T,typename U,
         typename std::enable_if<
	   is_NumericRectMatrix_tag<typename std::decay<T>::type::ET_tag>::value && is_NumericVector_tag<typename std::decay<U>::type::ET_tag>::value
				    , int>::type = 0>
inline auto operator*(T &&a, U &&b)->decltype( binaryHelper<ETnumericRectMatrixVectorMult,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)) ) {
  return binaryHelper<ETnumericRectMatrixVectorMult,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b));
}

template<typename Leaf1, typename Leaf2>
struct ETnumericRectSquareMatrixMult{
  typedef ETleafTag ET_leaf_mark;
  Leaf1 a;
  Leaf2 b;
  typedef typename Leaf1::ET_tag ET_tag;  //result is rect matrix

  ETnumericRectSquareMatrixMult(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){
    assert(a.common_properties().second == b.common_properties());
  }
  template<typename T>
  static inline auto elem(const T &m, const int i, const int j)->decltype(m[0]){ return m[j + m.common_properties().second*i]; }
  template<typename T>
  static inline auto belem(const T &m, const int i, const int j)->decltype(m[0]){ return m[j + m.common_properties()*i]; }

  inline typename std::decay<decltype(a[0])>::type operator[](const int i) const{
    typedef typename std::decay<decltype(a[0])>::type type;
    int nj = a.common_properties().second;
    int nk = nj;
    //out size   ni * nk    and indexing  k + nk*i

    int bk = i % nk;
    int ai = i / nk;

    type ret(a[0]);
    zeroit(ret);
    for(int j=0;j<nj;j++)
      ret = ret + elem(a,ai,j) * belem(b,j,bk);
    return ret;
  }
    
  inline decltype(a.common_properties()) common_properties() const{ return a.common_properties(); }
};
template<typename T,typename U,
         typename std::enable_if<
	   is_NumericRectMatrix_tag<typename std::decay<T>::type::ET_tag>::value && is_NumericSquareMatrix_tag<typename std::decay<U>::type::ET_tag>::value
				    , int>::type = 0>
inline auto operator*(T &&a, U &&b)->decltype( binaryHelper<ETnumericRectSquareMatrixMult,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)) ) {
  return binaryHelper<ETnumericRectSquareMatrixMult,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b));
}


SARLAC_END_NAMESPACE

