#pragma once

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include<config.h>
#include<utils/macros.h>
#include<serialize/hdf5_serialize.h>
#include<ET/generic_ET.h>
#include<tensors/numeric_vector.h>
#include<tensors/numeric_square_matrix.h>

SARLAC_START_NAMESPACE

template<typename Numeric>
class NumericRectMatrix{ //rectangular
  std::vector<Numeric> m;
  int n1;
  int n2;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & m & n1 & n2;
  }
public:
  typedef Numeric ElementType;

  NumericRectMatrix():m(){}
  explicit NumericRectMatrix(const int n1, const int n2): m(n1*n2), n1(n1), n2(n2){}
  explicit NumericRectMatrix(const std::pair<int,int> n): NumericRectMatrix(n.first,n.second){}

  NumericRectMatrix(const int n1, const int n2, const Numeric &init): m(n1*n2, init), n1(n1), n2(n2){}
  NumericRectMatrix(const NumericRectMatrix &r) = default;
  NumericRectMatrix(NumericRectMatrix &&r) = default;
  NumericRectMatrix(const int n1, const int n2, std::initializer_list<Numeric> l): NumericRectMatrix(n1,n2){ //j + n2*i
    assert(l.size() == n1*n2);
    auto it=l.begin();
    for(int i=0;i<n1*n2;i++) m[i] = *it++;
  }  
  template<typename Initializer> //Initializer is a lambda-type with operator()(const int, const int)
  inline NumericRectMatrix(const int n1, const int n2, const Initializer &initializer): NumericRectMatrix(n1,n2){
    for(int i=0;i<n1;i++)
      for(int j=0;j<n2;j++)
	m[j+n2*i] = initializer(i,j);
  }
  typedef NumericRectMatrix<Numeric> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,NumericRectMatrix<Numeric> >::value, int>::type = 0>
  NumericRectMatrix(U&& expr): NumericRectMatrix(expr.common_properties()){
#pragma omp parallel for
    for(int i=0;i<n1*n2;i++)
      m[i] = expr[i];
  }

  NumericRectMatrix & operator=(const NumericRectMatrix &r) = default;
  NumericRectMatrix & operator=(NumericRectMatrix &&r) = default;
    
  inline std::pair<int,int> size() const{ return {n1,n2}; }

  inline int rows() const{ return n1; }
  inline int cols() const{ return n2; }
  
  void resize(const int m1, const int m2){
    if(n1 == m1 && n2 == m2) return;
    this->n1 = m1; 
    this->n2 = m2;
    m.resize(m1*m2);
  }
  void resize(const int m1, const int m2, const Numeric &init){
    if(n1 == m1 && n2 == m2) return;    
    this->n1 = m1; 
    this->n2 = m2;
    m.resize(m1*m2, init);
  }
  template<typename Initializer>
  inline void resize(const int m1, const int m2, const Initializer &initializer){
    this->resize(m1,m2);
    for(int i=0;i<m1;i++){
      for(int j=0;j<m2;j++) 
	m[j+m2*i] = initializer(i,j);
    }
  }
  
  void zero(){
    for(int i=0;i<n1*n2;i++) m[i] = 0.;
  }

  std::string print() const{
    std::ostringstream os;
    for(int i=0;i<n1;i++){
      os << (*this)(i,0);
      for(int j=1;j<n2;j++)
  	os << " " << (*this)(i,j);
      os << std::endl;
    }
    return os.str();
  }
  
  Numeric & operator()(const int i, const int j){ return m[j+n2*i]; }
  const Numeric & operator()(const int i, const int j) const { return m[j+n2*i]; }

  Numeric & linElem(const int i){ return m[i]; }
  const Numeric & linElem(const int i) const { return m[i]; }

  NumericRectMatrix transpose() const{
    NumericRectMatrix out(this->n2, this->n1);
    for(int i=0;i<n1;i++)
      for(int j=0;j<n2;j++)
	out(j,i) = (*this)(i,j);
    return out;
  }

  template<typename T=Numeric, typename std::enable_if<is_std_complex<T>::value, int>::type = 0>
  NumericRectMatrix dagger() const{
    NumericRectMatrix out(this->n2, this->n1);
    for(int i=0;i<n1;i++)
      for(int j=0;j<n2;j++)
	out(j,i) = std::conj((*this)(i,j));
    return out;
  }

  NumericRectMatrix submatrix(const int istart, const int jstart, const int size1, const int size2) const{
    return NumericRectMatrix(size1, size2, [&](const int i, const int j){ return (*this)(i+istart,j+jstart); });
  }

  void extractRow(NumericVector<Numeric> &row, int row_idx) const{
    row.resize(this->n2);
    for(int j=0;j<this->n2;j++)
      row(j) = (*this)(row_idx,j);
  }

  void insertRow(const NumericVector<double> &row, int row_idx){
    assert(row.size() == this->n2);
    for(int j=0;j<this->n2;j++)
      (*this)(row_idx, j) = row(j);
  }

  GENERATE_HDF5_SERIALIZE_METHOD((m)(n1)(n2));
};


SARLAC_END_NAMESPACE

