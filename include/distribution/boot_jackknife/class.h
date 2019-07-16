#ifndef _BOOT_JACKKNIFE_CLASS_H_
#define _BOOT_JACKKNIFE_CLASS_H_

#include<distribution/jackknife.h>
#include<distribution/bootstrap.h>

//A bootstrap distribution of jackknife distributions! It fills a similar role as the double jackknife, allowing one to obtain an estimate of the covariance matrix in this case
//for each bootstrap sample in a bootstrap fit

CPSFIT_START_NAMESPACE

struct bootJackknifeInitType: public OstreamHook{
  int nsample;
  int boots;
  int confidence;
  bootJackknifeInitType(const int nsample, const int b = bootstrapDistributionOptions::defaultBoots(), const int c = bootstrapDistributionOptions::defaultConfidence()): nsample(nsample), boots(b), confidence(c){}
  inline bool operator==(const bootJackknifeInitType &r) const{ return nsample==r.nsample && boots==r.boots && confidence==r.confidence; }
  inline bool operator!=(const bootJackknifeInitType &r) const{ return !(*this == r); }
  void write(std::ostream &os) const{ os << "(nsample=" << nsample << ", boots=" << boots<< ", confidence=" << confidence <<")"; }
};



template<typename BaseDataType, template<typename> class BaseVectorType = basic_vector>
class bootJackknifeDistribution: public distribution<jackknifeDistribution<BaseDataType,BaseVectorType>, basic_vector >{
  bootstrapDistribution<BaseDataType> standardDeviation() const{ assert(0); };
  bootstrapDistribution<BaseDataType> standardError() const{ assert(0); };
  
  void resize(const int sz){ assert(0); }

  typedef distribution<jackknifeDistribution<BaseDataType,BaseVectorType>, basic_vector > baseType;
  typedef bootJackknifeDistribution<BaseDataType, BaseVectorType> myType;
  
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<baseType>(*this);
  }

  jackknifeDistribution<BaseDataType,BaseVectorType> base_jack; //propagate jackknife of original data through as "central value" 
  int _confidence;
public:
  typedef jackknifeDistribution<BaseDataType,BaseVectorType> DataType;
  
  template<typename T>
  using rebase = bootJackknifeDistribution<T,BaseVectorType>;

  typedef bootJackknifeInitType initType; //struct containing information used in constructor  

  void resize(const int nboots, const int nsample){ 
    this->_data.resize(nboots);
    for(int b=0;b<nboots;b++) this->sample(b).resize(nsample);
    base_jack.resize(nsample);
  }
  void resize(const initType &init){ 
    this->_confidence = init.confidence;
    this->baseType::resize(init.boots);
    for(int b=0;b<init.boots;b++) this->sample(b).resize(init.nsample);
    base_jack.resize(init.nsample);
  }

  inline bootJackknifeInitType getInitializer() const{ return bootJackknifeInitType(this->sample(0).size(), this->size(), this->_confidence); }

  int & confidence(){ return _confidence; }
  const int & confidence() const{ return _confidence; }
  
  inline const DataType & origEnsJackknife() const{ return base_jack; }
  inline DataType & origEnsJackknife(){ return base_jack; }
  
  inline const DataType & best() const{ return origEnsJackknife(); } 
  inline DataType & best(){ return origEnsJackknife(); }


  //Assumed to be a raw data distribution. Note confidence will be default value unless manually set. Boots and samples inferred from table
  template<typename DistributionType> 
  void resample(const DistributionType &in, const std::vector<std::vector<int> > &table){
    int boots = table.size();
    this->baseType::resize(boots);
    
    assert(table[0].size() <= in.size()); //table generation can discard some data
    int nraw = table[0].size();

    DistributionType in_r(nraw);

    for(int b=0;b<boots;b++){
      for(int s=0;s<nraw;s++) in_r.sample(s) = in.sample(table[b][s]); //generate a resampling of the raw data
      this->sample(b).resample(in_r);
    }

    base_jack.resample(in); //propagate jackknife of original data
  }
  
  bootJackknifeDistribution(): baseType(){}

  bootJackknifeDistribution(const bootJackknifeDistribution &r): _confidence(r._confidence), base_jack(r.base_jack), baseType(r){}
  
  template<template<typename> class U>
  bootJackknifeDistribution(const bootJackknifeDistribution<BaseDataType,U> &r): baseType(r.size()), _confidence(r._confidence){
    int nboot = r.size();
    int nsample = r.sample(0).size();

    for(int i=0;i<nboot;i++){
      this->sample(i).resize(nsample);
      for(int j=0;j<nsample;j++)
	this->sample(i).sample(j) = r.sample(i).sample(j);
    }
    base_jack.resize(nsample);
    for(int j=0;j<nsample;j++)
      base_jack.sample(j) = r.origEnsJackknife().sample(j);
  }
  
  explicit bootJackknifeDistribution(const initType &initv): _confidence(initv.confidence), base_jack(initv.nsample), baseType(initv.boots, DataType(initv.nsample)){}

  bootJackknifeDistribution(const initType &initv, const DataType &init):
    _confidence(initv.confidence), baseType(initv.boots,init), base_jack(init){ 
    assert(init.size() == initv.nsample);  
  }

  bootJackknifeDistribution(const initType &initv, const BaseDataType &init):
    _confidence(initv.confidence), baseType(initv.boots,DataType(initv.nsample,init)), base_jack(initv.nsample, init){
  }

  template<typename Initializer>
  bootJackknifeDistribution(const initType &initv, const Initializer &init): _confidence(initv.confidence), baseType(initv.boots,init){ 
    assert(this->sample(0).size() == initv.nsample);  
    base_jack = this->mean();
  }

  bootJackknifeDistribution(bootJackknifeDistribution&& o) noexcept : _confidence(o._confidence), base_jack(std::move(o.base_jack)), baseType(std::forward<baseType>(o)){}

  //Boots and samples inferred from table
  bootJackknifeDistribution(const rawDataDistribution<BaseDataType> &raw, const std::vector<std::vector<int> > &table, 
			    const int confidence = bootstrapDistributionOptions::defaultConfidence()): _confidence(confidence), baseType(){
    this->resample(raw, table);
  }
  
  //ET stuff
  typedef myType ET_tag;
  
  //common_properties() should return an initType
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,myType >::value, int>::type = 0>
  bootJackknifeDistribution(U&& expr): bootJackknifeDistribution(expr.common_properties()){
#pragma omp parallel for
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
    base_jack = expr[-1];
  }
  
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,myType >::value, int>::type = 0>
  myType & operator=(U&& expr){
    initType init = expr.common_properties();
    this->resize(init.boots, init.nsample);
    this->_confidence = init.confidence;
#pragma omp parallel for
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
    base_jack = expr[-1];
    return *this;
  }

  bootJackknifeDistribution & operator=(const bootJackknifeDistribution &r){ 
    _confidence = r._confidence;
    base_jack = r.base_jack ; 
    static_cast<baseType*>(this)->operator=(r); return *this; 
  }
  bootJackknifeDistribution & operator=(bootJackknifeDistribution &&r){ 
    _confidence = r._confidence;
    base_jack = std::move(r.base_jack); 
    static_cast<baseType*>(this)->operator=(std::move(r)); return *this; 
  }

  template<template<typename> class U = basic_vector>
  static bootstrapDistribution<BaseDataType,U> covariance(const myType &a, const myType &b){
    assert(a.size() == b.size());
    assert(a.confidence() == b.confidence());
    const int nboot = a.size();
    bootstrapDistribution<BaseDataType,U> out(bootstrapInitType(nboot, a.confidence()));
#pragma omp parallel for
    for(int i=0;i<nboot;i++)
      out.sample(i) = DataType::covariance(a.sample(i),b.sample(i));

    out.propagatedCentral() = DataType::covariance(a.origEnsJackknife(),b.origEnsJackknife()); //covariance for propagated base mean from jackknife of original unresampled data

    return out;
  }

  template<typename U=BaseDataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  bootJackknifeDistribution<typename U::value_type,BaseVectorType> real() const{
    bootJackknifeDistribution<typename U::value_type,BaseVectorType> out(getInitializer());
#pragma omp parallel for
    for(int i=0;i<this->size();i++)
      for(int j=0;j<this->sample(i).size();j++)
	out.sample(i).sample(j) = this->sample(i).sample(j).real(); 

    for(int i=0;i<base_jack.size();i++)
      out.origEnsJackknife().sample(i) = base_jack.sample(i).real();
   
    return out;
  }

  inline bool operator==(const bootJackknifeDistribution<BaseDataType,BaseVectorType> &r) const{ return this->base_jack == r.base_jack && this->baseType::operator==(r); }
  inline bool operator!=(const bootJackknifeDistribution<BaseDataType,BaseVectorType> &r) const{ return !( *this == r ); }
};

template<typename T, template<typename> class V>
std::ostream & operator<<(std::ostream &os, const bootJackknifeDistribution<T,V> &d){
  typedef distributionPrint<bootJackknifeDistribution<T,V> > printClass;
  assert(printClass::printer() != NULL); printClass::printer()->print(os, d);
  return os;
}

//Override basic printStats to print the central values and errors of the sub-distributions
template<typename T, template<typename> class V>
struct printStats< bootJackknifeDistribution<T,V> >{
  inline static std::string centralValue(const bootJackknifeDistribution<T,V> &d){
    std::ostringstream os; os << "[";
    for(int s=0;s<d.size()-1;s++) os << d.sample(s).best() << ", ";
    os << d.sample(d.size()-1).best() << "]";
    return os.str();
  }
  inline static std::string error(const bootJackknifeDistribution<T,V> &d){ 
    std::ostringstream os; os << "[";
    for(int s=0;s<d.size()-1;s++) os << d.sample(s).standardError() << ", ";
    os << d.sample(d.size()-1).standardError() << "]";
    return os.str();
  }

};


template<typename T, template<typename> class V>
bootJackknifeDistribution<T,V> weightedAvg(const std::vector<bootJackknifeDistribution<T,V> const*> &v){
  bootJackknifeInitType initv(v[0]->sample(0)->size(), v[0]->size(), v[0]->confidence());

  bootJackknifeDistribution<T,V> wavg_dj(initv);
  std::vector<jackknifeDistribution<T,V> const*> towavg_j(v.size());
  for(int s=0;s<initv.boots;s++){ //weighted avg each jackknife distribution, looping over outer index
    for(int i=0;i<v.size();i++) towavg_j[i] = &v[i]->sample(s);
    wavg_dj.sample(s) = CPSfit::weightedAvg(towavg_j);
  }
  
  //Jackknife of original ensemble
  for(int i=0;i<v.size();i++) towavg_j[i] = &v[i]->origEnsJackknife();
  wavg_dj.origEnsJackknife() = CPSfit::weightedAvg(towavg_j);

  return wavg_dj;
}

CPSFIT_END_NAMESPACE

#endif
