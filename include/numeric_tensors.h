#ifndef _NUMERIC_TENSORS_H
#define _NUMERIC_TENSORS_H

#include<template_wizardry.h>
#include<numeric_tensors_ET.h>


template<typename Numeric>
class NumericVector{
  std::vector<Numeric> v;
public:
  NumericVector():v(){}
  explicit NumericVector(const int n):v(n){}
  NumericVector(const int n, const Numeric &def):v(n,def){}

  template<typename Expr, IS_EXPRESSION_WITH_VECTOR_BASE_TYPE(Expr, NumericVector<Numeric>)>
  NumericVector(const Expr &e): v(e.size()){
    for(int i=0;i<e.size();i++)
      v[i] = e[i];
  }
  
  typedef VectorType<Numeric> Tensor_ET_base_type;
  
  int size() const{ return v.size(); }
  
  void resize(const int n){
    v.resize(n);
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

template<typename Numeric, typename StreamType, typename std::enable_if< isStreamType<StreamType>::value, int>::type = 0> 
StreamType & operator<<(StreamType & stream, const NumericVector<Numeric> &vec){
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

  typedef MatrixType<Numeric> Tensor_ET_base_type;
  
  int size() const{ return m.size(); }
  
  void resize(const int n){
    m.resize(n);
    for(int i=0;i<n;i++)
      m[i].resize(n);      
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

#endif
