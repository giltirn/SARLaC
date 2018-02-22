#ifndef _DISTRIBUTION_CLASS_H__
#define _DISTRIBUTION_CLASS_H__

#include<type_traits>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include<config.h>

#include<utils/utils.h>
#include<utils/template_wizardry.h>
#include<ET/generic_ET.h>
#include<serialize/hdf5_serialize.h>
#include<distribution/distribution_print.h>

CPSFIT_START_NAMESPACE

template<typename T>
using basic_vector = std::vector<T>;

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

  typedef distribution<_DataType,_VectorType> myType;
public:
  distribution(){}

  template<template<typename> class U>
  distribution(const distribution<DataType,U> &r): _data(r.sampleVector()){}

  distribution(const distribution &r): _data(r._data){}
  
  explicit distribution(const int nsample): _data(nsample){}
  distribution(const int nsample, const DataType &init): _data(nsample,init){}
  template<typename Initializer>
  distribution(const int nsample, const Initializer &init): _data(nsample){
    for(int i=0;i<nsample;i++) _data[i] = init(i);
  }	       
  distribution(distribution&& o) noexcept : _data(std::move(o._data)){}

  ENABLE_GENERIC_ET(distribution, myType, myType);
  
  distribution & operator=(const distribution &r){ _data = r._data; return *this; }
  
  int size() const{ return _data.size(); }

  void resize(const int sz){ _data.resize(sz); } //need default constructor for DataType
  void resize(const int sz, const DataType &init_val){ _data.resize(sz,init_val); }
  
  const VectorType<DataType> &sampleVector() const{ return _data; }
  VectorType<DataType> &sampleVector(){ return _data; }
  
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

CPSFIT_END_NAMESPACE

#endif
