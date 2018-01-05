#ifndef _DATA_CONTAINER_H
#define _DATA_CONTAINER_H

#include<config.h>
#include<vector>
#include<omp.h>
#include<cassert>
#include<cmath>
#include<type_traits>
#include<template_wizardry.h>
#include<hdf5_serialize.h>
#include<xml_serialize.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include<utils.h>
#include<generic_ET.h>
#include<distribution_print.h>
#include<distribution_vectors.h>


template<typename _DataType, template<typename> class _VectorType = basic_vector>
class distribution{
public:
  typedef _DataType DataType;
  template<typename T>
  using VectorType = _VectorType<T>;
  
  template<typename T>
  using rebase = distribution<T, VectorType>;
protected:  
  VectorType<DataType> _data;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & _data;
  }
#ifdef HAVE_HDF5
  template<typename DistributionType, typename std::enable_if<hasSampleMethod<DistributionType>::value, int>::type>
  friend void write(HDF5writer &writer, const DistributionType &value, const std::string &tag);
  template<typename DistributionType, typename std::enable_if<hasSampleMethod<DistributionType>::value, int>::type>
  friend void read(HDF5reader &reader, DistributionType &value, const std::string &tag);
#endif

  typedef distribution<_DataType,_VectorType> myType;
public:
  distribution(){}

  template<template<typename> class U>
  distribution(const distribution<DataType,U> &r): _data(r._data){}

  distribution(const distribution &r): _data(r._data){}
  
  explicit distribution(const int nsample): _data(nsample){}
  distribution(const int nsample, const DataType &init): _data(nsample,init){}
  template<typename Initializer>
  distribution(const int nsample, const Initializer &init): _data(nsample){
    for(int i=0;i<nsample;i++) _data[i] = init(i);
  }	       
  distribution(distribution&& o) noexcept : _data(std::move(o._data)){}

  ENABLE_GENERIC_ET(distribution, myType);
  
  distribution & operator=(const distribution &r){ _data = r._data; return *this; }
  
  int size() const{ return _data.size(); }

  void resize(const int sz){ _data.resize(sz); } //need default constructor for DataType
  void resize(const int sz, const DataType &init_val){ _data.resize(sz,init_val); }
  
  const VectorType<DataType> &sampleVector() const{ return _data; }
  
  inline const DataType & sample(const int idx) const{ return _data[idx]; }
  inline DataType & sample(const int idx){ return _data[idx]; }

  DataType mean() const{ 
    const int N = _data.size() == 0 ? 1 : _data.size();
    return threadedSum(_data)/double(N); 
  }
  
  DataType variance() const{
    const int N = _data.size() == 0 ? 1 : _data.size();
    const DataType avg = mean();

    const int nthread = omp_get_max_threads();

    struct Op{
      DataType avg;
      const VectorType<DataType> &data;
      Op(const DataType &_avg, const std::vector<DataType> &_data): data(_data), avg(_avg){}
      inline int size() const{ return data.size(); }
      inline DataType operator()(const int i) const{ DataType tmp = data[i] - avg; return tmp * tmp;  }
    };
    Op op(avg,_data);
    
    return threadedSum(op)/double(N);
  }
  DataType standardDeviation() const{
    return sqrt(variance()); //need sqrt defined, obviously
  }
  
  static DataType covariance(const distribution<DataType,VectorType> &a, const distribution<DataType,VectorType> &b){
    assert(a.size() == b.size());
    const int N = a.size() == 0 ? 1 : a.size();
    DataType avg_a = a.mean();
    DataType avg_b = b.mean();

    struct Op{
      DataType avg_a, avg_b;
      const VectorType<DataType> &data_a, &data_b;
      Op(const DataType &_avg_a, const DataType &_avg_b, const VectorType<DataType> &_data_a, const VectorType<DataType> &_data_b): data_a(_data_a), avg_a(_avg_a), data_b(_data_b), avg_b(_avg_b){}
      inline int size() const{ return data_a.size(); }
      inline DataType operator()(const int i) const{ return (data_a[i]-avg_a)*(data_b[i]-avg_b); }
    };
    Op op(avg_a,avg_b,a.sampleVector(),b.sampleVector());
    return threadedSum(op)/double(N);
  }

  template<typename U=DataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  distribution<typename U::value_type,VectorType> real() const{
    distribution<typename U::value_type,VectorType> out(this->size());
#pragma omp parallel for
    for(int i=0;i<this->size();i++) out.sample(i) = this->sample(i).real();
    return out;
  }

  inline void zero(){
    for(int i=0;i<this->size();i++) zeroit(this->sample(i));
  }

  inline bool operator==(const distribution<DataType,VectorType> &r) const{
    if(r.size() != this->size()) return false;
    for(int s=0;s<this->size();s++) if(r.sample(s) != this->sample(s)) return false;
    return true;
  }
  inline bool operator!=(const distribution<DataType,VectorType> &r) const{
    return !(*this == r);
  }

  void range(DataType &min, DataType &max) const{
    max = min = _data[0];
    for(int s=1;s<_data.size();s++){
      if(_data[s] > max) max = _data[s];
      if(_data[s] < min) min = _data[s];
    }    
  }
  
};

template<typename T, template<typename U> class V>
std::ostream & operator<<(std::ostream &os, const distribution<T,V> &d){
  typedef distributionPrint<distribution<T,V> > printClass;
  assert(printClass::printer() != NULL); printClass::printer()->print(os, d);
  return os;
}




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

  ENABLE_GENERIC_ET(rawDataDistribution, myType);
  
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


template<typename _DataType, template<typename> class _VectorType = basic_vector>
class doubleJackknifeDistribution;

template<typename _DataType, template<typename> class _VectorType = basic_vector>
class jackknifeDistribution: public distribution<_DataType,_VectorType>{
  _DataType variance() const{ assert(0); }
  _DataType standardDeviation() const{ assert(0); };
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<distribution<_DataType,_VectorType> >(*this);
  }
  typedef distribution<_DataType,_VectorType> baseType;
  typedef jackknifeDistribution<_DataType,_VectorType> myType;
public:
  typedef _DataType DataType;
  
  template<typename T>
  using VectorType = _VectorType<T>;
  
  template<typename T>
  using rebase = jackknifeDistribution<T,VectorType>;

  inline DataType best() const{ return this->mean(); }
  
  template<typename DistributionType> //doesn't have to be a distribution, just has to have a .sample and .size method
  void resample(const DistributionType &in){
    struct Op{
      const DistributionType &in;
      Op(const DistributionType &_in): in(_in){}
      inline int size() const{ return in.size(); }
      inline DataType operator()(const int i) const{ return in.sample(i); }
    };
    Op op(in);
    
    const int N = in.size();
    this->resize(N);

    DataType sum = threadedSum(op);
    
    const DataType den(N-1);
    const DataType num = DataType(1.)/den;
    
#pragma omp parallel for
    for(int j=0;j<N;j++)
      this->sample(j) = (sum - in.sample(j))*num;
  }

  DataType standardError() const{
    const int N = this->size();

    struct Op{
      DataType avg;
      const VectorType<DataType> &data;
      Op(const DataType &_avg, const VectorType<DataType> &_data): data(_data), avg(_avg){}
      inline int size() const{ return data.size(); }
      inline DataType operator()(const int i) const{ DataType tmp = data[i] - avg; return tmp * tmp;  }
    };
    Op op(this->mean(),this->_data);
    return sqrt( threadedSum(op) * (double(N-1)/N) );
  }
  
  jackknifeDistribution(): distribution<DataType,VectorType>(){}

  jackknifeDistribution(const jackknifeDistribution &r): baseType(r){}
  
  template<template<typename> class U>
  jackknifeDistribution(const jackknifeDistribution<DataType,U> &r): baseType(r){}
  
  explicit jackknifeDistribution(const int nsample): baseType(nsample){}
  jackknifeDistribution(const int nsample, const DataType &init): baseType(nsample,init){}
  template<typename Initializer>
  jackknifeDistribution(const int nsample, const Initializer &init): baseType(nsample,init){}
  jackknifeDistribution(jackknifeDistribution&& o) noexcept : baseType(std::forward<baseType>(o)){}

  template< template<typename> class U >
  jackknifeDistribution(const rawDataDistribution<DataType,U> &raw): jackknifeDistribution(raw.size()){ this->resample(raw); }
  
  ENABLE_GENERIC_ET(jackknifeDistribution, myType);
  
  jackknifeDistribution & operator=(const jackknifeDistribution &r){ static_cast<baseType*>(this)->operator=(r); return *this; }

  template<typename U=DataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  jackknifeDistribution<typename U::value_type> real() const{
    jackknifeDistribution<typename U::value_type> out(this->size());
#pragma omp parallel for
    for(int i=0;i<this->size();i++) out.sample(i) = this->sample(i).real();
    return out;
  }

  template<template<typename> class U = VectorType>
  doubleJackknifeDistribution<DataType,U> toDoubleJackknife() const;

  static DataType covariance(const jackknifeDistribution<DataType,VectorType> &a, const jackknifeDistribution<DataType,VectorType> &b){
    assert(a.size() == b.size());
    return distribution<DataType,VectorType>::covariance(a,b) * double(a.size()-1);  //like the standard error, the covariance of the jackknife samples is related to the covariance of the underlying distribution by a constant factor
  }

  inline bool operator==(const jackknifeDistribution<DataType,VectorType> &r) const{ return this->baseType::operator==(r); }
  inline bool operator!=(const jackknifeDistribution<DataType,VectorType> &r) const{ return !( *this == r ); }
};

template<typename T, template<typename> class V>
std::ostream & operator<<(std::ostream &os, const jackknifeDistribution<T,V> &d){
  typedef distributionPrint<jackknifeDistribution<T,V> > printClass;
  assert(printClass::printer() != NULL); printClass::printer()->print(os, d);
  return os;
}


//A jackknife distribution that independently propagates it's central value
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




template<typename BaseDataType, template<typename> class BaseVectorType>
class doubleJackknifeDistribution: public distribution<jackknifeDistribution<BaseDataType,BaseVectorType>, basic_vector >{
  jackknifeDistribution<BaseDataType> standardDeviation() const{ assert(0); };
  jackknifeDistribution<BaseDataType> standardError() const{ assert(0); };
  
  typedef distribution<jackknifeDistribution<BaseDataType,BaseVectorType>, basic_vector > baseType;
  typedef doubleJackknifeDistribution<BaseDataType, BaseVectorType> myType;
  
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<baseType>(*this);
  }

public:
  typedef jackknifeDistribution<BaseDataType,BaseVectorType> DataType;
  
  template<typename T>
  using rebase = doubleJackknifeDistribution<T,BaseVectorType>;
  
  template<typename DistributionType> //Assumed to be a raw data distribution
  void resample(const DistributionType &in){
    const int N = in.size();
    this->_data.resize(N);
    for(int i=0;i<N;i++) this->_data[i].resize(N-1);

    struct Op{
      const DistributionType &in;
      Op(const DistributionType &_in): in(_in){}
      inline int size() const{ return in.size(); }
      inline BaseDataType operator()(const int i) const{ return in.sample(i); }
    };
    Op op(in);
    BaseDataType sum = threadedSum(op);
    
    const double num = 1./double(N-2);
#pragma omp parallel for
    for(int i=0;i<N;i++){
      int jj=0;
      for(int j=0;j<N;j++){
	if(i==j) continue;

	this->sample(i).sample(jj++) = (sum - in.sample(i) - in.sample(j))*num;
      }
    }
  }

  doubleJackknifeDistribution(): baseType(){}

  doubleJackknifeDistribution(const doubleJackknifeDistribution &r): baseType(r){}
  
  template<template<typename> class U>
  doubleJackknifeDistribution(const doubleJackknifeDistribution<BaseDataType,U> &r): doubleJackknifeDistribution(r.size()){
    for(int i=0;i<this->size();i++)
      for(int j=0;j<this->size()-1;j++)
	this->sample(i).sample(j) = r.sample(i).sample(j);
  }

  
  explicit doubleJackknifeDistribution(const int nsample): baseType(nsample, DataType(nsample-1)){}
  doubleJackknifeDistribution(const int nsample, const DataType &init): baseType(nsample,init){}
  doubleJackknifeDistribution(const int nsample, const BaseDataType &init): baseType(nsample,DataType(nsample-1,init)){}
  template<typename Initializer>
  doubleJackknifeDistribution(const int nsample, const Initializer &init): baseType(nsample,init){}
  doubleJackknifeDistribution(doubleJackknifeDistribution&& o) noexcept : baseType(std::forward<baseType>(o)){}

  doubleJackknifeDistribution(const rawDataDistribution<BaseDataType> &raw): doubleJackknifeDistribution(raw.size()){
    this->resample(raw);
  }
  
  ENABLE_GENERIC_ET(doubleJackknifeDistribution, doubleJackknifeDistribution<BaseDataType>);

  doubleJackknifeDistribution & operator=(const doubleJackknifeDistribution &r){ static_cast<baseType*>(this)->operator=(r); return *this; }

  template<template<typename> class U = basic_vector>
  static jackknifeDistribution<BaseDataType,U> covariance(const myType &a, const myType &b){
    assert(a.size() == b.size());
    const int nouter = a.size();
    jackknifeDistribution<BaseDataType,U> out(nouter);
    for(int i=0;i<nouter;i++)
      out.sample(i) = DataType::covariance(a.sample(i),b.sample(i));
    return out;
  }
  template<typename U=BaseDataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  doubleJackknifeDistribution<typename U::value_type,BaseVectorType> real() const{
    doubleJackknifeDistribution<typename U::value_type,BaseVectorType> out(this->size());
#pragma omp parallel for
    for(int i=0;i<this->size();i++)
      for(int j=0;j<this->size()-1;j++)
	out.sample(i).sample(j) = this->sample(i).sample(j).real();
    return out;
  }

  jackknifeDistribution<BaseDataType> toJackknife() const{
    jackknifeDistribution<BaseDataType> out(this->size());
    for(int j=0;j<this->size();j++)
      out.sample(j) = this->sample(j).mean();
    return out;
  }

  inline bool operator==(const doubleJackknifeDistribution<BaseDataType,BaseVectorType> &r) const{ return this->baseType::operator==(r); }
  inline bool operator!=(const doubleJackknifeDistribution<BaseDataType,BaseVectorType> &r) const{ return !( *this == r ); }
};

template<typename T, template<typename> class V>
std::ostream & operator<<(std::ostream &os, const doubleJackknifeDistribution<T,V> &d){
  typedef distributionPrint<doubleJackknifeDistribution<T,V> > printClass;
  assert(printClass::printer() != NULL); printClass::printer()->print(os, d);
  return os;
}

template<typename T, template<typename> class V>
struct printStats< doubleJackknifeDistribution<T,V> >{
  inline static std::string centralValue(const doubleJackknifeDistribution<T,V> &d){
    std::ostringstream os; os << "[";
    for(int s=0;s<d.size()-1;s++) os << d.sample(s).best() << ", ";
    os << d.sample(d.size()-1).best() << "]";
    return os.str();
  }
  inline static std::string error(const doubleJackknifeDistribution<T,V> &d){ 
    std::ostringstream os; os << "[";
    for(int s=0;s<d.size()-1;s++) os << d.sample(s).standardError() << ", ";
    os << d.sample(d.size()-1).standardError() << "]";
    return os.str();
  }

};

template<typename DataType, template<typename> class V>
template<template<typename> class U>
doubleJackknifeDistribution<DataType,U> jackknifeDistribution<DataType,V>::toDoubleJackknife() const{
  DataType Jbar = this->mean();
  int N = this->size();
  doubleJackknifeDistribution<DataType,U> out(N);
  for(int j=0;j<N;j++)
    for(int k=0;k<N-1;k++){
      int kp = k < j ? k : k+1;      
      out.sample(j).sample(k) = double(N-1)/(N-2) * ( this->sample(j) + this->sample(kp) - double(N)/(N-1)*Jbar );
    }
  return out;  
}

#include<distribution_IO.h>
#include<distribution_ET.h>

#endif
