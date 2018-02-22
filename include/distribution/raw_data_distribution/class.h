#ifndef _RAW_DATA_DISTRIBUTION_CLASS_H_
#define _RAW_DATA_DISTRIBUTION_CLASS_H_

#include<config.h>
#include<distribution/distribution.h>

CPSFIT_START_NAMESPACE

template<typename _DataType, template<typename> class _VectorType = basic_vector>
class rawDataDistribution: public distribution<_DataType, _VectorType>{
  friend class boost::serialization::access;
  typedef distribution<_DataType,_VectorType> baseType;
  typedef rawDataDistribution<_DataType,_VectorType> myType;
  
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<baseType>(*this);
  }
public:
  typedef _DataType DataType;
  
  template<typename T>
  using rebase = rawDataDistribution<T,_VectorType>;

  inline DataType best() const{ return this->mean(); } //"central value" of distribution used for printing/plotting
  DataType standardError() const{ return this->standardDeviation()/sqrt(double(this->size()-1)); }
  
  rawDataDistribution(): distribution<DataType>(){}

  rawDataDistribution(const rawDataDistribution &r): baseType(r){}
  
  template<template<typename> class U>
  rawDataDistribution(const rawDataDistribution<DataType,U> &r): baseType(r){}

  
  explicit rawDataDistribution(const int nsample): baseType(nsample){}
  rawDataDistribution(const int nsample, const DataType &init): baseType(nsample,init){}
  template<typename Initializer>
  rawDataDistribution(const int nsample, const Initializer &init): baseType(nsample,init){}
  rawDataDistribution(rawDataDistribution&& o) noexcept : baseType(std::forward<baseType>(o)){}

  ENABLE_GENERIC_ET(rawDataDistribution, myType, myType);
  
  rawDataDistribution & operator=(const rawDataDistribution &r){ static_cast<baseType*>(this)->operator=(r); return *this; }

  template<typename U=DataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  rawDataDistribution<typename U::value_type,_VectorType> real() const{
    rawDataDistribution<typename U::value_type,_VectorType> out(this->size());
#pragma omp parallel for
    for(int i=0;i<this->size();i++) out.sample(i) = this->sample(i).real();
    return out;
  }

  inline bool operator==(const rawDataDistribution<DataType,_VectorType> &r) const{ return this->baseType::operator==(r); }
  inline bool operator!=(const rawDataDistribution<DataType,_VectorType> &r) const{ return !( *this == r ); }

  rawDataDistribution<DataType,_VectorType> bin(const int bin_size) const{
    const int nsample = this->size();
    const int nbins = nsample / bin_size;
    if(nsample % bin_size != 0) 
      error_exit(std::cout << "rawDataDistribution::bin(const int) distribution size " << nsample << " is not divisible by bin size " << bin_size << std::endl);

    DataType zro; zeroit(zro);
    rawDataDistribution<DataType,_VectorType> out(nbins, zro);

    for(int b=0;b<nbins;b++){
      for(int e=0;e<bin_size;e++)
	out.sample(b) = out.sample(b) + this->sample(e + bin_size*b);
      out.sample(b) = out.sample(b) / double(bin_size);
    }
    return out;
  }

};

template<typename T, template<typename> class _VectorType = basic_vector>
std::ostream & operator<<(std::ostream &os, const rawDataDistribution<T,_VectorType> &d){
  typedef distributionPrint<rawDataDistribution<T,_VectorType> > printClass;
  assert(printClass::printer() != NULL); printClass::printer()->print(os, d);
  return os;
}

CPSFIT_END_NAMESPACE
#endif
