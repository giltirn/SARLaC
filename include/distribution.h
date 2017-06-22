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

#include<generic_ET.h>

template<typename Operation, typename T>
struct _threadedSumHelper{
  //For getting an initial zero optimally we consider several scenarios:
  //1) It is constructible using a single double or float
  //2) It is default constructible and has a "zero" method
  //3) It is default constructible and does not have a "zero" method but has an operator=(double/float) method
  //4) It is not default constructible and has a "zero" method
  //5) It is not default constructible and does not have a "zero" method but has an operator=(double/float) method

  //Having a copy constructor is assumed
  
  enum { float_constructible = std::is_constructible<T,double>::value || std::is_constructible<T,float>::value,
	 default_constructible = std::is_default_constructible<T>::value,
	 has_zero_method = hasZeroMethod<T>::value,
	 has_equals_method = hasEqualsMethod<T,double>::value || hasEqualsMethod<T,float>::value 
  };

  enum { class1 = float_constructible };
  enum { class2 = !class1 && default_constructible && has_zero_method };
  enum { class3 = !class1 && !class2 && default_constructible && has_equals_method };
  enum { class4 = !class1 && !class2 && !class3 && !default_constructible && has_zero_method };
  enum { class5 = !class1 && !class2 && !class3 && !class4 && !default_constructible && has_equals_method };

  template<typename U>
  struct _truewrap{ enum { value = 1 }; }; //SFINAE only works for deduced types
  
  template<typename U=T, typename std::enable_if<_truewrap<U>::value && class1,int>::type = 0>
    static inline T getZero(const Operation &op){ return T(0.); }

  template<typename U=T, typename std::enable_if<_truewrap<U>::value && class2,int>::type = 0>
  static inline T getZero(const Operation &op){ T ret; ret.zero(); return ret; }

  template<typename U=T, typename std::enable_if<_truewrap<U>::value && class3,int>::type = 0>
  static inline T getZero(const Operation &op){ T ret; ret = 0.; return ret; }
  
  template<typename U=T, typename std::enable_if<_truewrap<U>::value && class4,int>::type = 0>
  static inline T getZero(const Operation &op){ T ret(op(0)); ret.zero(); return ret; }

  template<typename U=T, typename std::enable_if<_truewrap<U>::value && class5,int>::type = 0>
  static inline T getZero(const Operation &op){ T ret(op(0)); ret = 0.; return ret; }
};


template<typename Operation>
auto threadedSum(const Operation &op)->typename std::decay<decltype(op(0))>::type{
  typedef typename std::decay<decltype(op(0))>::type T;
  const int N = op.size();
  const int nthread = omp_get_max_threads();
 
  T init_zero = _threadedSumHelper<Operation,T>::getZero(op);
  std::vector<T> sum(nthread, init_zero);

#pragma omp parallel for
  for(int i=0;i<N;i++)
    sum[omp_get_thread_num()] = sum[omp_get_thread_num()] + op(i);

  for(int t=1;t<nthread;t++)
    sum[0] = sum[0] + sum[t];

  return sum[0];
}

template<typename T>
T threadedSum(const std::vector<T> &v){
  struct Op{
    const std::vector<T> &vv;
    inline const T & operator()(const int idx) const{ return vv[idx]; }
    inline size_t size() const{ return vv.size(); }
    Op(const std::vector<T> &_vv): vv(_vv){}
  };
  Op op(v);
  return threadedSum<Op>(op);
}


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

  inline DataType best() const{ return mean(); } //"central value" of distribution used for printing/plotting
  
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
  
  DataType standardError() const{ return standardDeviation()/sqrt(double(_data.size()-1)); }

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
    for(int i=0;i<this->size();i++) out.sample(i) = this->sample(i).real();
    return out;
  }

  inline void zero(){
    for(int i=0;i<this->size();i++) zeroit(this->sample(i));
  }
  
};

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
    for(int i=0;i<this->size();i++) out.sample(i) = this->sample(i).real();
    return out;
  }
};

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
  jackknifeCdistribution(U&& expr){
    this->resize(expr.common_properties());
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
    cen = expr[-1];
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
    for(int i=0;i<this->size();i++) out.sample(i) = this->sample(i).real();
    return out;
  }
};



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
    for(int i=0;i<this->size();i++)
      for(int j=0;j<this->size()-1;j++)
	out.sample(i).sample(j) = this->sample(i).sample(j).real();
    return out;
  }
};

#include<distribution_ET.h>
#include<distribution_print.h>
#endif
