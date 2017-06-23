#ifndef _NUMERIC_TENSORS_H
#define _NUMERIC_TENSORS_H

#include<map>
#include<template_wizardry.h>
#include<generic_ET.h>


template<typename Numeric>
class NumericVector{
  std::vector<Numeric> v;
public:
  NumericVector():v(){}
  explicit NumericVector(const int n):v(n){}
  NumericVector(const int n, const Numeric &def):v(n,def){}
  NumericVector(const NumericVector &r) = default;
  NumericVector(NumericVector &&r) = default;

  NumericVector & operator=(const NumericVector &r) = default;
  NumericVector & operator=(NumericVector &&r) = default;
  
  ENABLE_GENERIC_ET(NumericVector, NumericVector<Numeric>);
  
  int size() const{ return v.size(); }
  
  void resize(const int n){
    v.resize(n);
  }
  void resize(const int n, const Numeric &init){
    v.resize(n,init);
  }
  void zero(){ for(int i=0;i<v.size();i++) v[i] = 0.; }

  std::string print() const{
    std::ostringstream os;
    os << v[0];
    for(int i=1;i<v.size();i++)
      os << " " << v[i];
    return os.str();
  }

  Numeric & operator()(const int i){ return v[i]; }
  const Numeric & operator()(const int i) const { return v[i]; }

  Numeric & operator[](const int i){ return v[i]; }
  const Numeric & operator[](const int i) const { return v[i]; }
  
  NumericVector &operator+=(const NumericVector &r){
    for(int i=0;i<v.size();i++) v[i] += r.v[i];
  }
};

template<typename Numeric> 
std::ostream & operator<<(std::ostream & stream, const NumericVector<Numeric> &vec){
  stream << "(";
  for(int i=0;i<vec.size();i++)
    stream << vec[i] << (i != vec.size()-1 ? " " : ")");
  return stream;
}


template<typename Numeric>
class SVDinvertPolicy;



template<typename Numeric, typename InvertPolicy = SVDinvertPolicy<Numeric> >
class NumericMatrix: public InvertPolicy{ //square matrix
  std::vector<std::vector<Numeric> > m;
public:
  NumericMatrix():m(){}
  explicit NumericMatrix(const int n): m(n, std::vector<Numeric>(n)){}
  NumericMatrix(const int n, const Numeric &init): m(n, std::vector<Numeric>(n,init)){}
  NumericMatrix(const NumericMatrix &r) = default;
  NumericMatrix(NumericMatrix &&r) = default;

  NumericMatrix & operator=(const NumericMatrix &r) = default;
  NumericMatrix & operator=(NumericMatrix &&r) = default;
  
  typedef NumericMatrix<Numeric,InvertPolicy> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,NumericMatrix<Numeric,InvertPolicy> >::value, int>::type = 0>
  NumericMatrix(U&& expr): NumericMatrix(expr.common_properties()){
    for(int i=0;i<this->size()*this->size();i++)
      getElem<NumericMatrix<Numeric,InvertPolicy> >::elem(*this, i) = expr[i];
  }
  
  inline int size() const{ return m.size(); }
  
  void resize(const int n){
    if(m.size() == n) return;
    m.resize(n);
    for(int i=0;i<n;i++)
      m[i].resize(n);      
  }
  void resize(const int n, const Numeric &init){
    if(m.size() == n) return;
    m.resize(n);
    for(int i=0;i<n;i++)
      m[i].resize(n,init);     
  }
  void zero(){
    const int n = m.size();
    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
	m[i][j] = 0.;
  }
  void invert(const NumericMatrix<Numeric> &what){
    this->InvertPolicy::invert(*this,what);
  }

  std::string print() const{
    std::ostringstream os;
    for(int i=0;i<m.size();i++){
      os << m[i][0];
      for(int j=1;j<m[i].size();j++)
  	os << " " << m[i][j];
      os << std::endl;
    }
    return os.str();
  }
  
  Numeric & operator()(const int i, const int j){ return m[i][j]; }
  const Numeric & operator()(const int i, const int j) const { return m[i][j]; }
};



template<typename Numeric, typename StreamType, typename std::enable_if< isStreamType<StreamType>::value, int>::type = 0> 
StreamType & operator<<(StreamType & stream, const NumericMatrix<Numeric> &mat){
  for(int i=0;i<mat.size();i++){
    for(int j=0;j<mat.size();j++){
      stream << mat(i,j) << " ";
    }
    stream << '\n';
  }
  return stream;
}
		     



template<typename Numeric>
class SVDinvertPolicy{
 protected:
  inline static void invert(NumericMatrix<Numeric> &inv_m, const NumericMatrix<Numeric> &m){
    inv_m.resize(m.size());
    svd_inverse(inv_m, m);
  }
};



template<typename NumericMatrixType>
class NumericMatrixSampleView{
  typedef typename _get_elem_type<NumericMatrixType>::type DistributionType;
  typedef typename std::remove_const<typename std::remove_reference<decltype( ((DistributionType*)(NULL))->sample(0) )>::type>::type SampleType;

  NumericMatrixType &M;
  int sample;
public:
  NumericMatrixSampleView(NumericMatrixType &_M, const int _sample): M(_M), sample(_sample){}
  
  inline int size() const{ return M.size(); }

  inline const SampleType& operator()(const int i, const int j) const{ return M(i,j).sample(sample); }
  
  template<typename U = NumericMatrixType>
  inline typename std::enable_if< !std::is_const<U>::value, SampleType& >::type operator()(const int i, const int j){ return M(i,j).sample(sample); }
};


//A vector type with an additional global mapping between some generic tag and it's elements
template<typename T, typename MapType>
class mappedVector: public NumericVector<T>{
  //MapType must conform to the following form
  // struct MapTypeExample{
  //   typedef std::string tagType;    //choose a tag type - here a string but can be anything
  //   inline int map(const tagType &tag) const; //mapping between tag and index
  //   inline tagType unmap(const int idx) const; //unmapping between index and tag
  //   inline int size() const; //number of elements
  // };

  MapType const* mapping;

public:
  typedef typename MapType::tagType tagType;

  mappedVector(const MapType & _mapping): mapping(&_mapping), NumericVector<T>(_mapping.size()){}
  
  typedef mappedVector<T,MapType> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,mappedVector<T,MapType> >::value, int>::type = 0>
  mappedVector(U&& expr){
    mapping = expr.common_properties();
    this->resize(mapping->size());
    for(int i=0;i<mapping->size();i++) this->operator()(i) = expr[i];
  }
  
  inline T & operator()(const tagType &tag){ return this->operator()(mapping->map(tag)); }
  inline const T & operator()(const tagType &tag) const{ return this->operator()(mapping->map(tag)); }
  inline T & operator()(const int idx){ return NumericVector<T>::operator()(idx); }
  inline const T & operator()(const int idx) const{ return NumericVector<T>::operator()(idx); }
  
  const MapType & getMapping() const{ return *mapping; }
};
template<typename T, typename MapType>
struct getElem<mappedVector<T,MapType> >{
  static inline auto elem(const mappedVector<T,MapType> &v, const int i)->decltype(v(i)){ return v(i); }
  // static inline auto elem(mappedVector<T,MapType> &v, const int i)->decltype(v(i)){ return v(i); }
  static inline MapType const* common_properties(const mappedVector<T,MapType> &v){ return &v.getMapping(); }
};

template<typename T, typename MapType>
inline void debug_print(const mappedVector<T,MapType> &v){ std::cout << &v.getMapping() << std::endl; std::cout.flush(); }


class stringTagMap{
  std::map<std::string,int> mp;
  std::map<int,std::string> ump;
public:
  typedef std::string tagType; 
  inline int map(const tagType &tag) const{
    auto it = mp.find(tag);
    assert(it != mp.end());
    return it->second;
  }
  inline tagType unmap(const int idx) const{
    auto it = ump.find(idx);
    assert(it != ump.end());
    return it->second;
  }
  void add(const std::string &tag, const int idx){
    mp[tag] = idx;
    ump[idx] = tag;
  }
  inline int size() const{ return mp.size(); }
};



#include<numeric_tensors_ET.h>

#endif
