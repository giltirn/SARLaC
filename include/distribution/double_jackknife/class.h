#ifndef _DOUBLE_JACKKNIFE_CLASS_H_
#define _DOUBLE_JACKKNIFE_CLASS_H_

#include<distribution/jackknife.h>

//A distribution for double-elimination jackknife data

CPSFIT_START_NAMESPACE


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
    BaseDataType sum; zeroit(sum);
    for(int i=0;i<N;i++){
      sum = sum + in.sample(i);
      this->_data[i].resize(N-1);
    }
    
    const double num = 1./double(N-2);
#pragma omp parallel for
    for(int i=0;i<N;i++){
      for(int j=0;j<i;j++)
	this->sample(i).sample(j) = (sum - in.sample(i) - in.sample(j))*num;
      int jj=i;
      for(int j=i+1;j<N;j++)
	this->sample(i).sample(jj++) = (sum - in.sample(i) - in.sample(j))*num;
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
  
  ENABLE_PARALLEL_GENERIC_ET(doubleJackknifeDistribution, myType, doubleJackknifeDistribution<BaseDataType>);

  doubleJackknifeDistribution & operator=(const doubleJackknifeDistribution &r){ static_cast<baseType*>(this)->operator=(r); return *this; }
  doubleJackknifeDistribution & operator=(doubleJackknifeDistribution &&r){ static_cast<baseType*>(this)->operator=(std::move(r)); return *this; }

  template<template<typename> class U = basic_vector>
  static jackknifeDistribution<BaseDataType,U> covariance(const myType &a, const myType &b){
    assert(a.size() == b.size());
    const int nouter = a.size();
    jackknifeDistribution<BaseDataType,U> out(nouter);
#pragma omp parallel for
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

  //Note this relationship is incorrect (by, I think 1/N^2 effects) if this distribution is related non-linearly to the original resampled distribution
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

//Override basic printStats to print the central values and errors of the sub-distributions
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

//Convert jackknife to double-jackknife (be warned, this is only correct providing you have not performed any non-linear operations on the raw data or jackknife data)
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



CPSFIT_END_NAMESPACE

#endif
