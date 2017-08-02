#ifndef _DATA_CONTAINER_H
#define _DATA_CONTAINER_H

#include<vector>
#include<omp.h>
#include<cassert>
#include<cmath>
#include<type_traits>
#include<template_wizardry.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include<utils.h>
#include<generic_ET.h>
#include<distribution_print.h>


template<typename _DataType>
class distribution{
public:
  typedef _DataType DataType;

  template<typename T>
  using rebase = distribution<T>;
protected:  
  std::vector<DataType> _data;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & _data;
  }
public:
  distribution(){}
  distribution(const distribution &r): _data(r._data){}
  explicit distribution(const int nsample): _data(nsample){}
  distribution(const int nsample, const DataType &init): _data(nsample,init){}
  distribution(distribution&& o) noexcept : _data(std::move(o._data)){}

  ENABLE_GENERIC_ET(distribution, distribution<_DataType>);
  
  distribution & operator=(const distribution &r){ _data = r._data; return *this; }
  
  int size() const{ return _data.size(); }

  void resize(const int sz){ _data.resize(sz); } //need default constructor for DataType
  void resize(const int sz, const DataType &init_val){ _data.resize(sz,init_val); }
  
  const std::vector<DataType> &sampleVector() const{ return _data; }
  
  inline const DataType & sample(const int idx) const{ return _data[idx]; }
  inline DataType & sample(const int idx){ return _data[idx]; }

  DataType mean() const{ return threadedSum(_data)/double(_data.size()); }
  
  DataType variance() const{
    const int N = _data.size();
    const DataType avg = mean();

    const int nthread = omp_get_max_threads();

    struct Op{
      DataType avg;
      const std::vector<DataType> &data;
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
  
  static DataType covariance(const distribution<DataType> &a, const distribution<DataType> &b){
    assert(a.size() == b.size());
    DataType avg_a = a.mean();
    DataType avg_b = b.mean();

    struct Op{
      DataType avg_a, avg_b;
      const std::vector<DataType> &data_a, &data_b;
      Op(const DataType &_avg_a, const DataType &_avg_b, const std::vector<DataType> &_data_a, const std::vector<DataType> &_data_b): data_a(_data_a), avg_a(_avg_a), data_b(_data_b), avg_b(_avg_b){}
      inline int size() const{ return data_a.size(); }
      inline DataType operator()(const int i) const{ return (data_a[i]-avg_a)*(data_b[i]-avg_b); }
    };
    Op op(avg_a,avg_b,a.sampleVector(),b.sampleVector());
    return threadedSum(op)/double(a.size());
  }

  template<typename U=DataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  distribution<typename U::value_type> real() const{
    distribution<typename U::value_type> out(this->size());
#pragma omp parallel for
    for(int i=0;i<this->size();i++) out.sample(i) = this->sample(i).real();
    return out;
  }

  inline void zero(){
    for(int i=0;i<this->size();i++) zeroit(this->sample(i));
  }
};

template<typename T>
std::ostream & operator<<(std::ostream &os, const distribution<T> &d){
  assert(distributionPrint<distribution<T> >::printer() != NULL); distributionPrint<distribution<T> >::printer()->print(os, d);
  return os;
}


template<typename _DataType>
class rawDataDistribution: public distribution<_DataType>{
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<distribution<_DataType> >(*this);
  }
  typedef distribution<_DataType> baseType;
public:
  typedef _DataType DataType;
  
  template<typename T>
  using rebase = rawDataDistribution<T>;

  inline DataType best() const{ return this->mean(); } //"central value" of distribution used for printing/plotting
  DataType standardError() const{ return this->standardDeviation()/sqrt(double(this->size()-1)); }
  
  rawDataDistribution(): distribution<DataType>(){}
  rawDataDistribution(const rawDataDistribution &r): baseType(r){}
  explicit rawDataDistribution(const int nsample): baseType(nsample){}
  rawDataDistribution(const int nsample, const DataType &init): baseType(nsample,init){}
  rawDataDistribution(rawDataDistribution&& o) noexcept : baseType(std::forward<baseType>(o)){}

  ENABLE_GENERIC_ET(rawDataDistribution, rawDataDistribution<_DataType>);
  
  rawDataDistribution & operator=(const rawDataDistribution &r){ static_cast<distribution<DataType>*>(this)->operator=(r); return *this; }

  template<typename U=DataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  rawDataDistribution<typename U::value_type> real() const{
    rawDataDistribution<typename U::value_type> out(this->size());
#pragma omp parallel for
    for(int i=0;i<this->size();i++) out.sample(i) = this->sample(i).real();
    return out;
  }
};

template<typename T>
std::ostream & operator<<(std::ostream &os, const rawDataDistribution<T> &d){
  assert(distributionPrint<rawDataDistribution<T> >::printer() != NULL); distributionPrint<rawDataDistribution<T> >::printer()->print(os, d);
  return os;
}


template<typename _DataType>
class doubleJackknifeDistribution;

template<typename _DataType>
class jackknifeDistribution: public distribution<_DataType>{
  _DataType variance() const{ assert(0); }
  _DataType standardDeviation() const{ assert(0); };
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<distribution<_DataType> >(*this);
  }
  typedef distribution<_DataType> baseType;
public:
  typedef _DataType DataType;
  
  template<typename T>
  using rebase = jackknifeDistribution<T>;

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
      const std::vector<DataType> &data;
      Op(const DataType &_avg, const std::vector<DataType> &_data): data(_data), avg(_avg){}
      inline int size() const{ return data.size(); }
      inline DataType operator()(const int i) const{ DataType tmp = data[i] - avg; return tmp * tmp;  }
    };
    Op op(this->mean(),this->_data);
    return sqrt( threadedSum(op) * (double(N-1)/N) );
  }
  
  jackknifeDistribution(): distribution<DataType>(){}
  jackknifeDistribution(const jackknifeDistribution &r): baseType(r){}
  explicit jackknifeDistribution(const int nsample): baseType(nsample){}
  jackknifeDistribution(const int nsample, const DataType &init): baseType(nsample,init){}
  jackknifeDistribution(jackknifeDistribution&& o) noexcept : baseType(std::forward<baseType>(o)){}

  ENABLE_GENERIC_ET(jackknifeDistribution, jackknifeDistribution<_DataType>);
  
  jackknifeDistribution & operator=(const jackknifeDistribution &r){ static_cast<distribution<DataType>*>(this)->operator=(r); return *this; }

  template<typename U=DataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  jackknifeDistribution<typename U::value_type> real() const{
    jackknifeDistribution<typename U::value_type> out(this->size());
#pragma omp parallel for
    for(int i=0;i<this->size();i++) out.sample(i) = this->sample(i).real();
    return out;
  }

  doubleJackknifeDistribution<DataType> toDoubleJackknife() const;

  static DataType covariance(const jackknifeDistribution<DataType> &a, const jackknifeDistribution<DataType> &b){
    assert(a.size() == b.size());
    return distribution<DataType>::covariance(a,b) * double(a.size()-1);  //like the standard error, the covariance of the jackknife samples is related to the covariance of the underlying distribution by a constant factor
  }


};

template<typename T>
std::ostream & operator<<(std::ostream &os, const jackknifeDistribution<T> &d){
  assert(distributionPrint<jackknifeDistribution<T> >::printer() != NULL); distributionPrint<jackknifeDistribution<T> >::printer()->print(os, d);
  return os;
}


//A jackknife distribution that independently propagates it's central value
template<typename _DataType>
class jackknifeCdistribution: public jackknifeDistribution<_DataType>{
  _DataType cen;
  typedef jackknifeDistribution<_DataType> baseType;
public:
  typedef _DataType DataType;
  
  template<typename T>
  using rebase = jackknifeCdistribution<T>;
  
  jackknifeCdistribution() = default;
  jackknifeCdistribution(const jackknifeCdistribution &r) = default;
  explicit jackknifeCdistribution(const int nsample): baseType(nsample){}
  jackknifeCdistribution(const int nsample, const DataType &init): baseType(nsample,init), cen(init){}
  jackknifeCdistribution(jackknifeCdistribution&& o) noexcept : baseType(std::forward<baseType>(o)){}

  typedef jackknifeCdistribution<DataType> ET_tag;
  
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,jackknifeCdistribution<DataType> >::value, int>::type = 0>
  jackknifeCdistribution(U&& expr): jackknifeCdistribution(expr.common_properties()){
#ifdef PARALLELIZE_DISTRIBUTION_ET
#pragma omp parallel for
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
#else
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
#endif
    cen = expr[-1];
  }
  
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,jackknifeCdistribution<DataType> >::value, int>::type = 0>
  jackknifeCdistribution<DataType> & operator=(U&& expr){
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
    this->jackknifeDistribution<_DataType>::resample(in);
    cen = in.mean();
  }
  
#ifdef UKVALENCE_COMPAT
  //UKvalence uses the propagated central value when computing the standard error
  DataType standardError() const{
    const int N = this->size();

    struct Op{
      const DataType &avg;
      const std::vector<DataType> &data;
      Op(const DataType &_avg, const std::vector<DataType> &_data): data(_data), avg(_avg){}
      inline int size() const{ return data.size(); }
      inline DataType operator()(const int i) const{ DataType tmp = data[i] - avg; return tmp * tmp;  }
    };
    Op op(this->cen,this->_data);
    return sqrt( threadedSum(op) * (double(N-1)/N) );
  }
#endif

  void import(const jackknifeDistribution<DataType> &jack){
    this->resize(jack.size());
    for(int i=0;i<this->size();i++) this->sample(i) = jack.sample(i);
    cen = jack.mean();
  }
  template<typename U=DataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  jackknifeCdistribution<typename U::value_type> real() const{
    jackknifeCdistribution<typename U::value_type> out(this->size());
#pragma omp parallel for
    for(int i=0;i<this->size();i++) out.sample(i) = this->sample(i).real();
    return out;
  }
};

template<typename T>
std::ostream & operator<<(std::ostream &os, const jackknifeCdistribution<T> &d){
  assert(distributionPrint<jackknifeCdistribution<T> >::printer() != NULL); distributionPrint<jackknifeCdistribution<T> >::printer()->print(os, d);
  return os;
}


template<typename BaseDataType>
class doubleJackknifeDistribution: public distribution<jackknifeDistribution<BaseDataType> >{
  jackknifeDistribution<BaseDataType> standardDeviation() const{ assert(0); };
  jackknifeDistribution<BaseDataType> standardError() const{ assert(0); };
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<distribution<jackknifeDistribution<BaseDataType> > >(*this);
  }
  typedef distribution<jackknifeDistribution<BaseDataType> > baseType;
public:
  typedef jackknifeDistribution<BaseDataType> DataType;
  
  template<typename T>
  using rebase = doubleJackknifeDistribution<T>;
  
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
  
  explicit doubleJackknifeDistribution(const int nsample): baseType(nsample, DataType(nsample-1)){}
  doubleJackknifeDistribution(const int nsample, const DataType &init): baseType(nsample,init){}
  doubleJackknifeDistribution(const int nsample, const BaseDataType &init): baseType(nsample,DataType(nsample-1,init)){}
  doubleJackknifeDistribution(doubleJackknifeDistribution&& o) noexcept : baseType(std::forward<baseType>(o)){}

  ENABLE_GENERIC_ET(doubleJackknifeDistribution, doubleJackknifeDistribution<BaseDataType>);

  doubleJackknifeDistribution & operator=(const doubleJackknifeDistribution &r){ static_cast<baseType*>(this)->operator=(r); return *this; }

  static jackknifeDistribution<BaseDataType> covariance(const doubleJackknifeDistribution<BaseDataType> &a, const doubleJackknifeDistribution<BaseDataType> &b){
    assert(a.size() == b.size());
    const int nouter = a.size();
    jackknifeDistribution<BaseDataType> out(nouter);
    for(int i=0;i<nouter;i++)
      out.sample(i) = jackknifeDistribution<BaseDataType>::covariance(a.sample(i),b.sample(i));
    return out;
  }
  template<typename U=BaseDataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  doubleJackknifeDistribution<typename U::value_type> real() const{
    doubleJackknifeDistribution<typename U::value_type> out(this->size());
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
};

template<typename T>
std::ostream & operator<<(std::ostream &os, const doubleJackknifeDistribution<T> &d){
  assert(distributionPrint<doubleJackknifeDistribution<T> >::printer() != NULL); distributionPrint<doubleJackknifeDistribution<T> >::printer()->print(os, d);
  return os;
}

template<typename T>
struct printStats< doubleJackknifeDistribution<T> >{
  inline static std::string centralValue(const doubleJackknifeDistribution<T> &d){
    std::ostringstream os; os << "[";
    for(int s=0;s<d.size()-1;s++) os << d.sample(s).best() << ", ";
    os << d.sample(d.size()-1).best() << "]";
    return os.str();
  }
  inline static std::string error(const doubleJackknifeDistribution<T> &d){ 
    std::ostringstream os; os << "[";
    for(int s=0;s<d.size()-1;s++) os << d.sample(s).standardError() << ", ";
    os << d.sample(d.size()-1).standardError() << "]";
    return os.str();
  }

};

template<typename DataType>
doubleJackknifeDistribution<DataType> jackknifeDistribution<DataType>::toDoubleJackknife() const{
  DataType Jbar = this->mean();
  int N = this->size();
  doubleJackknifeDistribution<DataType> out(N);
  for(int j=0;j<N;j++)
    for(int k=0;k<N-1;k++){
      int kp = k < j ? k : k+1;      
      out.sample(j).sample(k) = double(N-1)/(N-2) * ( this->sample(j) + this->sample(kp) - double(N)/(N-1)*Jbar );
    }
  return out;  
}


#include<distribution_ET.h>

#endif
