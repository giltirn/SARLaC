#ifndef _CPSFIT_NUMERIC_VECTOR_CLASS_H_
#define _CPSFIT_NUMERIC_VECTOR_CLASS_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include<config.h>
#include<utils/macros.h>
#include<serialize/hdf5_serialize.h>
#include<ET/generic_ET.h>

CPSFIT_START_NAMESPACE

template<typename Numeric>
class NumericVector{
  std::vector<Numeric> v;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & v;
  }
public:
  NumericVector():v(){}
  explicit NumericVector(const int n):v(n){}
  NumericVector(const int n, const Numeric &def):v(n,def){}
  NumericVector(const NumericVector &r) = default;
  NumericVector(NumericVector &&r) = default;
  
  template<typename Initializer> //Initializer is a lambda-type with operator()(const int)
  inline NumericVector(const int n, const Initializer &initializer): v(n){
    for(int i=0;i<n;i++) v[i] = initializer(i);
  }
  
  NumericVector & operator=(const NumericVector &r) = default;
  NumericVector & operator=(NumericVector &&r) = default;

  ENABLE_GENERIC_ET(NumericVector, NumericVector<Numeric>, NumericVector<Numeric>);
  
  int size() const{ return v.size(); }
  
  void resize(const int n){
    v.resize(n);
  }
  void resize(const int n, const Numeric &init){
    v.resize(n,init);
  }
  template<typename Initializer>
  inline void resize(const int n, const Initializer &initializer){
    v.resize(n);
    for(int i=0;i<n;i++) v[i] = initializer(i);
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

  inline void push_back(const Numeric &b){ v.push_back(b); }

  NumericVector subvector(const int start, const int size) const{
    return NumericVector(size, [&](const int i){ return v[i+start]; });
  }

  GENERATE_HDF5_SERIALIZE_METHOD((v));
};

CPSFIT_END_NAMESPACE

#endif
