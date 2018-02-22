#ifndef _JACKKNIFE_C__CLASS_H_
#define _JACKKNIFE_C__CLASS_H_

#include<distribution/jackknife.h>

CPSFIT_START_NAMESPACE

template<typename _DataType, template<typename> class _VectorType = basic_vector>
class jackknifeCdistribution: public jackknifeDistribution<_DataType,_VectorType>{
  _DataType cen;
  typedef jackknifeDistribution<_DataType,_VectorType> baseType;
  typedef jackknifeCdistribution<_DataType,_VectorType> myType;
  
#ifdef HAVE_HDF5
  template<typename T, template<typename> class V>
  friend void write(HDF5writer &writer, const jackknifeCdistribution<T,V> &value, const std::string &tag);
  template<typename T, template<typename> class V>
  friend void read(HDF5reader &reader, jackknifeCdistribution<T,V> &value, const std::string &tag);
#endif
public:
  typedef _DataType DataType;
  
  template<typename T>
  using rebase = jackknifeCdistribution<T,_VectorType>;
  
  jackknifeCdistribution() = default;
  jackknifeCdistribution(const jackknifeCdistribution &r) = default;
  
  template<template<typename> class U>
  jackknifeCdistribution(const jackknifeCdistribution<DataType,U> &r): baseType(r), cen(r.cen){}
  
  explicit jackknifeCdistribution(const int nsample): baseType(nsample){}
  jackknifeCdistribution(const int nsample, const DataType &init): baseType(nsample,init), cen(init){}
  template<typename Initializer>
  jackknifeCdistribution(const int nsample, const Initializer &init): baseType(nsample,init){
    cen = this->mean();
  }
  jackknifeCdistribution(jackknifeCdistribution&& o) noexcept : baseType(std::forward<baseType>(o)){}

  template<template<typename> class U>
  jackknifeCdistribution(const rawDataDistribution<DataType,U> &raw): jackknifeCdistribution(raw.size()){ this->resample(raw); }
  
  typedef myType ET_tag;
  
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,myType >::value, int>::type = 0>
  jackknifeCdistribution(U&& expr): jackknifeCdistribution(expr.common_properties()){
#ifdef PARALLELIZE_DISTRIBUTION_ET
#pragma omp parallel for
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
#else
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
#endif
    cen = expr[-1];
  }
  
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,myType >::value, int>::type = 0>
  myType & operator=(U&& expr){
    this->resize(expr.common_properties());
#ifdef PARALLELIZE_DISTRIBUTION_ET
#pragma omp parallel for
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
#else
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
#endif
    cen = expr[-1];
    return *this;
  }

  jackknifeCdistribution & operator=(const jackknifeCdistribution &r){ static_cast<baseType*>(this)->operator=(r); cen=r.cen; return *this; }

  inline const DataType & propagatedCentral() const{ return cen; }
  inline DataType & propagatedCentral(){ return cen; }
  
  inline const DataType & best() const{ return propagatedCentral(); } //"central value" of distribution used for printing/plotting
  inline DataType & best(){ return propagatedCentral(); }

  inline const DataType & sample(const int idx) const{ return idx == -1 ? this->best() : this->baseType::sample(idx); }
  inline DataType & sample(const int idx){ return idx == -1 ? this->best() : this->baseType::sample(idx); }

  template<typename DistributionType>
  void resample(const DistributionType &in){
    this->baseType::resample(in);
    cen = in.mean();
  }
  
#ifdef UKVALENCE_COMPAT
  //UKvalence uses the propagated central value when computing the standard error
  DataType standardError() const{
    const int N = this->size();

    struct Op{
      const DataType &avg;
      const _VectorType<DataType> &data;
      Op(const DataType &_avg, const _VectorType<DataType> &_data): data(_data), avg(_avg){}
      inline int size() const{ return data.size(); }
      inline DataType operator()(const int i) const{ DataType tmp = data[i] - avg; return tmp * tmp;  }
    };
    Op op(this->cen,this->_data);
    return sqrt( threadedSum(op) * (double(N-1)/N) );
  }
#endif

  template<template<typename> class U>
  void import(const jackknifeDistribution<DataType,U> &jack){
    this->resize(jack.size());
    for(int i=0;i<this->size();i++) this->sample(i) = jack.sample(i);
    cen = jack.mean();
  }
  template<typename U=DataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  jackknifeCdistribution<typename U::value_type,_VectorType> real() const{
    jackknifeCdistribution<typename U::value_type,_VectorType> out(this->size());
#pragma omp parallel for
    for(int i=0;i<this->size();i++) out.sample(i) = this->sample(i).real();
    return out;
  }

  inline bool operator==(const myType &r) const{ return (this->propagatedCentral() == r.propagatedCentral()) && this->baseType::operator==(r); }
  inline bool operator!=(const myType &r) const{ return !( *this == r ); }
};

template<typename T, template<typename> class V>
std::ostream & operator<<(std::ostream &os, const jackknifeCdistribution<T,V> &d){
  typedef distributionPrint<jackknifeCdistribution<T,V> > printClass;
  assert(printClass::printer() != NULL); printClass::printer()->print(os, d);
  return os;
}


CPSFIT_END_NAMESPACE

#endif
