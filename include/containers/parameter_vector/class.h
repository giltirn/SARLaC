#ifndef _CPSFIT_PARAMETER_VECTOR_CLASS_H_
#define _CPSFIT_PARAMETER_VECTOR_CLASS_H_

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_vector.h>

CPSFIT_START_NAMESPACE

template<typename T>
class parameterVector: public NumericVector<T>{
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){ static_cast<NumericVector<T>* >(this)->serialize(ar,version); }

public:
  parameterVector() = default;
  explicit parameterVector(const int n): NumericVector<T>(n){}
  parameterVector(const int n, const T &def): NumericVector<T>(n,def){}
  parameterVector(const parameterVector &r) = default;
  parameterVector(parameterVector &&r) = default;
  
  template<typename Initializer> //Initializer is a lambda-type with operator()(const int)
  inline parameterVector(const int n, const Initializer &initializer): NumericVector<T>(n,initializer){ }

  parameterVector(std::initializer_list<T> l): NumericVector<T>(l){}
  
  parameterVector & operator=(const parameterVector &r) = default;
  parameterVector & operator=(parameterVector &&r) = default;
  parameterVector & operator=(std::initializer_list<T> l){ this->NumericVector<T>::operator=(l); return *this; }      

  ENABLE_GENERIC_ET(parameterVector, parameterVector<T>, parameterVector<T>);
  
  parameterVector &operator+=(const parameterVector &r){ return static_cast<NumericVector<T>* >(this)->operator+=(r); }

  void write(HDF5writer &writer, const std::string &tag) const{ static_cast<NumericVector<T> const* >(this)->write(writer,tag); }

  void read(HDF5reader &reader, const std::string &tag){ static_cast<NumericVector<T>* >(this)->read(reader,tag); }
};


CPSFIT_END_NAMESPACE
#endif
