#ifndef _DATA_CONTAINER_H
#define _DATA_CONTAINER_H

#include<vector>
#include<omp.h>
#include<cassert>
#include<cmath>
#include<template_wizardry.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

template<typename T>
T threadedSum(const std::vector<T> &v){
  const int N = v.size();
  const int nthread = omp_get_max_threads();
    
  T sum[nthread]; for(int i=0;i<nthread;i++) sum[i] = 0.;
#pragma omp parallel for
  for(int i=0;i<N;i++)
    sum[omp_get_thread_num()] += v[i];

  for(int t=1;t<nthread;t++)
    sum[0] += sum[t];

  return sum[0];
}

template<typename Operation>
auto threadedSum(const Operation &op)->decltype( op(0) ){
  typedef decltype( op(0) ) T;
  const int N = op.size();
  const int nthread = omp_get_max_threads();
    
  T sum[nthread]; for(int i=0;i<nthread;i++) sum[i] = 0.;
#pragma omp parallel for
  for(int i=0;i<N;i++)
    sum[omp_get_thread_num()] += op(i);

  for(int t=1;t<nthread;t++)
    sum[0] += sum[t];

  return sum[0];
}


template<typename _DataType>
class distribution{
public:
  typedef _DataType DataType;
protected:  
  std::vector<DataType> _data;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & _data;
  }
public:
  distribution(){}
  distribution(const int nsample): _data(nsample){}
  distribution(const int nsample, const DataType &init): _data(nsample,init){}
  
  int size() const{ return _data.size(); }

  void resize(const int sz){ _data.resize(sz); } //need default constructor for DataType

  const std::vector<DataType> &sampleVector() const{ return _data; }
  
  const DataType & sample(const int idx) const{ return _data[idx]; }
  DataType & sample(const int idx){ return _data[idx]; }

  DataType mean() const{ return threadedSum(_data)/_data.size(); }

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
public:
  typedef _DataType DataType;

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
    return sqrt( threadedSum(op) * double(N-1)/N );
  }
  
  jackknifeDistribution(): distribution<DataType>(){}
  jackknifeDistribution(const int nsample): distribution<DataType>(nsample){}
  jackknifeDistribution(const int nsample, const DataType &init): distribution<DataType>(nsample,init){}
};




#endif
