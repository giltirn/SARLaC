#ifndef _JACKKNIFE_CLASS_H_
#define _JACKKNIFE_CLASS_H_

#include<distribution/raw_data_distribution.h>

//A distribution for single-elimination jackknife data

SARLAC_START_NAMESPACE

template<typename _DataType, template<typename> class _VectorType = basic_vector>
class doubleJackknifeDistribution;

template<typename _DataType, template<typename> class _VectorType = basic_vector>
class jackknifeDistribution: public distribution<_DataType,_VectorType>{
  _DataType variance() const{ assert(0); }

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
  using rebase = jackknifeDistribution<T,_VectorType>;

  typedef int initType;
  
  inline int getInitializer() const{ return this->size(); }

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

  //Note: this is the standard deviation of the jackknife samples, not a property of the population mean
  inline DataType standardDeviation() const{ return this->baseType::standardDeviation(); };
  
  jackknifeDistribution(): baseType(){}

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
  
  ENABLE_GENERIC_ET(jackknifeDistribution, myType, myType);
  
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

  //Convert back to raw data. **NOTE** This is only valid if only linear operations have been performed to the jackknife samples
  template<template<typename> class U = VectorType>
  rawDataDistribution<DataType, U> toRaw() const{
    int N = this->size();
    rawDataDistribution<DataType, U> out(N);
    DataType mean = this->mean();
    for(int j=0;j<N;j++)
      out.sample(j) = N*mean - (N-1)*this->sample(j);
    return out;
  }


  //Covariance of the *means*
  static DataType covariance(const jackknifeDistribution<DataType,VectorType> &a, const jackknifeDistribution<DataType,VectorType> &b){
    assert(a.size() == b.size());
    return distribution<DataType,VectorType>::covariance(a,b) * double(a.size()-1);  //like the standard error, the covariance of the jackknife samples is related to the covariance of the underlying distribution by a constant factor
  }

  inline bool operator==(const jackknifeDistribution<DataType,VectorType> &r) const{ return this->baseType::operator==(r); }
  inline bool operator!=(const jackknifeDistribution<DataType,VectorType> &r) const{ return !( *this == r ); }

  inline void push_back(const DataType &v){ this->sampleVector().push_back(v); }
};

template<typename T, template<typename> class V>
std::ostream & operator<<(std::ostream &os, const jackknifeDistribution<T,V> &d){
  typedef distributionPrint<jackknifeDistribution<T,V> > printClass;
  assert(printClass::printer() != NULL); printClass::printer()->print(os, d);
  return os;
}

template<typename T>
struct is_jackknife{
  enum {value = 0};
};
template<typename T, template<typename> class V>
struct is_jackknife<jackknifeDistribution<T,V> >{
  enum {value=1};
};

//It is straightforward to show that a superjackknife can be treated in the same way as a regular jackknife up to 1/N effects
//This function implements the "boost" of a regular jackknife on subensemble 'ens' to a superjackknife, filling in the gaps with the mean
template<typename T, template<typename> class V>
jackknifeDistribution<T,V> superjackknifeBoost(const jackknifeDistribution<T,V> &v, const int ens, const std::vector<int> &subens_sizes){
  T mu = v.mean();
  int N = 0;
  int ens_off = 0;
  for(int i=0;i<subens_sizes.size();i++){
    if(i == ens) ens_off = N;
    N += subens_sizes[i];
  }
  int Nin = subens_sizes[ens];
  assert(Nin == v.size());
  
  jackknifeDistribution<T,V> out(N,mu);
  for(int s=ens_off;s<ens_off + Nin;s++)
    out.sample(s) = v.sample(s-ens_off);

  return out;
}




SARLAC_END_NAMESPACE

#endif
