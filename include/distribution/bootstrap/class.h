#ifndef _BOOTSTRAP_CLASS_H_
#define _BOOTSTRAP_CLASS_H_

#include<distribution/raw_data_distribution.h>
#include<random/rng.h>
#include<utils/utils/text_output.h>
//A distribution for bootstrap data

CPSFIT_START_NAMESPACE

template<typename _DataType, template<typename> class _VectorType = basic_vector>
class bootstrapDistribution: public distribution<_DataType,_VectorType>{
  _DataType variance() const{ assert(0); }
  _DataType standardDeviation() const{ assert(0); };

  typedef distribution<_DataType,_VectorType> baseType;
  typedef bootstrapDistribution<_DataType,_VectorType> myType;

#ifdef HAVE_HDF5
  template<typename T, template<typename> class V>
  friend void write(HDF5writer &writer, const bootstrapDistribution<T,V> &value, const std::string &tag);
  template<typename T, template<typename> class V>
  friend void read(HDF5reader &reader, bootstrapDistribution<T,V> &value, const std::string &tag);
#endif

  _DataType avg; //average value, propagated separately
  int _confidence;
public:
  typedef _DataType DataType;
  
  template<typename T>
  using VectorType = _VectorType<T>;
  
  template<typename T>
  using rebase = bootstrapDistribution<T,_VectorType>;

  inline const DataType & propagatedCentral() const{ return avg; }
  inline DataType & propagatedCentral(){ return avg; }
  
  inline const DataType & best() const{ return avg; } //"central value" of distribution used for printing/plotting
  inline DataType & best(){ return avg; }

  template<typename DistributionType> //doesn't have to be a distribution, just has to have a .sample and .size method
  void resample(const DistributionType &in){
    if(!RNG.isInitialized()) error_exit(std::cout << "bootstrapDistribution::resample RNG is not initialized" << std::endl);
    
    this->avg = in.best();
    
    int nraw = in.size();
    std::uniform_int_distribution<> dis(0,nraw-1);

    int boots = this->size();

    assert(!omp_in_parallel());

    for(int b=0;b<boots;b++){
      zeroit(this->sample(b));

      //Draw N samples with replacement
      for(int i=0;i<nraw;i++)
	this->sample(b) = this->sample(b) + in.sample(dis(RNG()));
      
      this->sample(b) = this->sample(b)/nraw;
    }
  }

  int & confidence(){ return _confidence; }
  const int & confidence() const{ return _confidence; }

  //Return the lower and upper bounds of the confidence region (lo, hi)
  std::pair<DataType,DataType> confidenceRegion() const{
    const int N = this->size();

    std::vector<_DataType> sorted(N);
    for(int i=0;i<N;i++) sorted[i] = this->sample(i);
    std::sort(sorted.begin(),sorted.end());

    int omit=(100-_confidence)/2;
    int lo_idx = omit*N/100;
    int hi_idx = N-1-lo_idx;

    return std::pair<DataType,DataType>(sorted[lo_idx], sorted[hi_idx]);
  }

  //return the lower and upper error bars (lower,upper)
  inline std::pair<DataType,DataType> errorBounds() const{
    std::pair<DataType,DataType> bounds = this->confidenceRegion();  
    bounds.first = avg - bounds.first;
    bounds.second = bounds.second - avg;
    return bounds;
  }

  //standardError assuming symmetric errors
  inline DataType standardError() const{
    std::pair<DataType,DataType> bounds = this->confidenceRegion();      
    //( (hi-best) + (best-lo) )/2   
    return (bounds.second-bounds.first)/2.;
  }
 
  bootstrapDistribution(): baseType(), avg(0.0), _confidence(68){}
  bootstrapDistribution(const bootstrapDistribution &r) = default;
  
  template<template<typename> class U>
  bootstrapDistribution(const bootstrapDistribution<DataType,U> &r): baseType(r), _confidence(r._confidence), avg(r.avg){}
  
  struct initType: public OstreamHook{
    int boots;
    int confidence;
    initType(const int b, const int c): boots(b), confidence(c){}
    inline bool operator==(const initType &r) const{ return boots==r.boots && confidence==r.confidence; }
    inline bool operator!=(const initType &r) const{ return !(*this == r); }
    void write(std::ostream &os) const{ os << "(boots=" << boots<< ", confidence=" << confidence <<")"; }
  };
  bootstrapDistribution(const initType &init): baseType(init.boots), _confidence(init.confidence){}

  bootstrapDistribution(const int boots, const int confidence = 68): baseType(boots), _confidence(confidence){}

  bootstrapDistribution(const int boots, const DataType &init, const int confidence = 68): baseType(boots,init), _confidence(confidence){ avg = this->mean(); }

  template<typename Initializer>
  bootstrapDistribution(const int boots, const Initializer &init, const int confidence = 68): baseType(boots,init), _confidence(confidence){ avg = this->mean(); }

  bootstrapDistribution(bootstrapDistribution&& o) noexcept : baseType(std::forward<baseType>(o)), _confidence(o._confidence), avg(o.avg){}

  template< template<typename> class U >
  bootstrapDistribution(const rawDataDistribution<DataType,U> &raw, const int boots, const int confidence = 68): bootstrapDistribution(boots,confidence){ this->resample(raw); }
  

  typedef myType ET_tag;
  
  //common_properties() should return an initType
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,myType >::value, int>::type = 0>
    bootstrapDistribution(U&& expr): bootstrapDistribution(expr.common_properties()){
#ifdef PARALLELIZE_DISTRIBUTION_ET
#pragma omp parallel for
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
#else
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
#endif
    avg = expr[-1];
  }
  
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,myType >::value, int>::type = 0>
  myType & operator=(U&& expr){
    initType init = expr.common_properties();
    this->resize(init.boots);
    this->_confidence = init.confidence;
#ifdef PARALLELIZE_DISTRIBUTION_ET
#pragma omp parallel for
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
#else
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
#endif
    avg = expr[-1];
    return *this;
  }
  
  bootstrapDistribution & operator=(const bootstrapDistribution &r){ _confidence = r._confidence; avg=r.avg; static_cast<baseType*>(this)->operator=(r); return *this; }

  template<typename U=DataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  bootstrapDistribution<typename U::value_type> real() const{
    bootstrapDistribution<typename U::value_type> out(this->size(),this->_confidence);
#pragma omp parallel for
    for(int i=0;i<this->size();i++) out.sample(i) = this->sample(i).real();
    out.avg = this->avg;
    return out;
  }

  inline bool operator==(const bootstrapDistribution<DataType,VectorType> &r) const{ return r._confidence == _confidence && r.avg == avg && this->baseType::operator==(r); }
  inline bool operator!=(const bootstrapDistribution<DataType,VectorType> &r) const{ return !( *this == r ); }
};

template<typename T, template<typename> class V>
std::ostream & operator<<(std::ostream &os, const bootstrapDistribution<T,V> &d){
  typedef distributionPrint<bootstrapDistribution<T,V> > printClass;
  assert(printClass::printer() != NULL); printClass::printer()->print(os, d);
  return os;
}

CPSFIT_END_NAMESPACE

#endif
