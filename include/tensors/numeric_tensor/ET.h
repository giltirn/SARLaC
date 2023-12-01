#ifndef _SARLAC_NUMERIC_TENSOR_ET_H_
#define _SARLAC_NUMERIC_TENSOR_ET_H_

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_tensor/class.h>

SARLAC_START_NAMESPACE

//Numeric tensor stuff
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

//Matrix-matrix multiplication
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



//Matrix-vector multiplication
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


SARLAC_END_NAMESPACE
#endif
