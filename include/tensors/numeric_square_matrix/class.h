#ifndef _CPSFIT_NUMERIC_SQUARE_MATRIX_CLASS_H_
#define _CPSFIT_NUMERIC_SQUARE_MATRIX_CLASS_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include<config.h>
#include<utils/macros.h>
#include<serialize/hdf5_serialize.h>
#include<ET/generic_ET.h>

CPSFIT_START_NAMESPACE

template<typename Numeric>
class NumericSquareMatrix{ //square matrix
  std::vector<std::vector<Numeric> > m;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & m;
  }
public:
  typedef Numeric ElementType;

  NumericSquareMatrix():m(){}
  explicit NumericSquareMatrix(const int n): m(n, std::vector<Numeric>(n)){}
  NumericSquareMatrix(const int n, const Numeric &init): m(n, std::vector<Numeric>(n,init)){}
  NumericSquareMatrix(const NumericSquareMatrix &r) = default;
  NumericSquareMatrix(NumericSquareMatrix &&r) = default;
  NumericSquareMatrix(std::initializer_list<Numeric> l): m(int(floor(sqrt(l.size()))), std::vector<Numeric>(int(floor(sqrt(l.size())))) ){
    assert(l.size() == this->size() * this->size());
    auto it=l.begin();
    for(int i=0;i<m.size();i++)
      for(int j=0;j<m.size();j++)
	m[i][j] = *it++;
  }  
  template<typename Initializer> //Initializer is a lambda-type with operator()(const int, const int)
  inline NumericSquareMatrix(const int n, const Initializer &initializer): m(n, std::vector<Numeric>(n)){
    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
	m[i][j] = initializer(i,j);
  }
  typedef NumericSquareMatrix<Numeric> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,NumericSquareMatrix<Numeric> >::value, int>::type = 0>
  NumericSquareMatrix(U&& expr): NumericSquareMatrix(expr.common_properties()){
#pragma omp parallel for
    for(int i=0;i<this->size()*this->size();i++)
      getElem<NumericSquareMatrix<Numeric> >::elem(*this, i) = expr[i];
  }

  NumericSquareMatrix & operator=(const NumericSquareMatrix &r) = default;
  NumericSquareMatrix & operator=(NumericSquareMatrix &&r) = default;
    
  NumericSquareMatrix & operator=(std::initializer_list<Numeric> l){
    int sz = int(floor(sqrt(l.size())));
    assert(l.size() == sz*sz);
    auto it=l.begin();
    m.resize(sz);
    for(int i=0;i<m.size();i++){
      m[i].resize(sz);
      for(int j=0;j<m.size();j++)
	m[i][j] = *it++;
    }
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
  template<typename Initializer>
  inline void resize(const int n, const Initializer &initializer){
    m.resize(n);
    for(int i=0;i<n;i++){
      m[i].resize(n);
      for(int j=0;j<n;j++) m[i][j] = initializer(i,j);
    }
  }
  
  void zero(){
    const int n = m.size();
    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
	m[i][j] = 0.;
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

  NumericSquareMatrix transpose() const{
    NumericSquareMatrix out(this->size());
    for(int i=0;i<this->size();i++)
      for(int j=0;j<this->size();j++)
	out(i,j) = (*this)(j,i);
    return out;
  }

  NumericSquareMatrix submatrix(const int istart, const int jstart, const int size) const{
    return NumericSquareMatrix(size, [&](const int i, const int j){ return m[i+istart][j+jstart]; });
  }

  GENERATE_HDF5_SERIALIZE_METHOD((m));
};


CPSFIT_END_NAMESPACE
#endif
