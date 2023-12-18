#ifndef _BLOCK_DOUBLE_JACKKNIFE_CLASS_H_
#define _BLOCK_DOUBLE_JACKKNIFE_CLASS_H_

#include<distribution/jackknife.h>

//This is a variant of the double-jackknife intended for estimating covariance matrices from binned data
//If the bin size is b, we drop the b samples associated with the row index, and then consider the single-elimination jackknife of the remainder

//D_ij = ( N<d> - [ \sum_{k = b*i}^{b*(i+1)} d_k ] - d_j ) / (N-b-1)

//******
//NOTE: If original nsample is not a multiple of bin_size, the number of samples used will be cropped to the nearest multiple
//******


SARLAC_START_NAMESPACE


template<typename BaseDataType, template<typename> class BaseVectorType = basic_vector>
class blockDoubleJackknifeDistribution: public distribution<jackknifeDistribution<BaseDataType,BaseVectorType>, basic_vector >{
public:
  typedef jackknifeDistribution<BaseDataType,BaseVectorType> DataType;
  
  template<typename T>
  using rebase = blockDoubleJackknifeDistribution<T,BaseVectorType>;

  typedef std::pair<int,int> initType; //(nsample, bin_size)
private:
  jackknifeDistribution<BaseDataType> standardDeviation() const{ assert(0); };
  jackknifeDistribution<BaseDataType> standardError() const{ assert(0); };
  void resize(const int sz){ assert(0); }
  void resize(const int sz, const DataType &init_val){ assert(0); }
  
  typedef distribution<jackknifeDistribution<BaseDataType,BaseVectorType>, basic_vector > baseType;
  typedef blockDoubleJackknifeDistribution<BaseDataType, BaseVectorType> myType;
  
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<baseType>(*this);
  }

  int bin_size;
  int nsample; 

  inline static int binCrop(const int nsample, const int bin_size){
    return  (nsample / bin_size) * bin_size;
  }

public:
  inline initType getInitializer() const{ return initType(nsample,bin_size); }

  inline int nSamplesUnbinned() const{
    return nsample;
  }
  inline int binSize() const{
    return bin_size;
  }
  
  void resize(const int _nsample, const int _bin_size){
    nsample = binCrop(_nsample, _bin_size);
    bin_size = _bin_size;

    int nbinned = nsample / bin_size;
    int nsub = nsample - bin_size;

    this->_data.resize(nbinned);
    for(int i=0;i<nbinned;i++) this->_data[i].resize(nsub);
  }
  inline void resize(const initType &nb){ //used by ET
    this->resize(nb.first, nb.second);
  }

  template<typename DistributionType> //Assumed to be a raw data distribution
  void resample(const DistributionType &in){
    int in_nsample_crop = binCrop(in.size(), this->bin_size);
    assert(this->nsample == in_nsample_crop);

    //Number of samples may be cropped if in.size() is not a multiple of bin_size
    BaseDataType Nmean = in.sample(0);
    for(int i=1;i<nsample;i++) Nmean = Nmean + in.sample(i); //only include up to cropping

    double nrm = 1./(nsample - bin_size - 1); 

#pragma omp parallel for
    for(int i=0;i<this->size();i++){

      int bin_start = i*bin_size;
      int bin_lessthan = bin_start + bin_size;

      BaseDataType bsum = in.sample(bin_start);
      for(int k=bin_start+1;k<bin_lessthan;k++)
	bsum = bsum + in.sample(k);

      BaseDataType vbase_i = Nmean - bsum;

      //D_ij = ( N<d> - [ \sum_{k = b*i}^{b*(i+1)} d_k ] - d_j ) / (N-b-1)
      int j = 0;
      for(int j_true = 0; j_true < bin_start; j_true++)
	this->sample(i).sample(j++) = (vbase_i - in.sample(j_true)) * nrm;
      for(int j_true = bin_lessthan; j_true < nsample; j_true++)
	this->sample(i).sample(j++) = (vbase_i - in.sample(j_true)) * nrm;
    }
  }

  //This version allows the user to provide a new bin size or set it if the default constructor was used
  template<typename DistributionType> //Assumed to be a raw data distribution
  void resample(const DistributionType &in, const int bin_size){
    this->resize(in.size(), bin_size);
    this->resample(in);
  }

  blockDoubleJackknifeDistribution(): baseType(), bin_size(1), nsample(0){}

  blockDoubleJackknifeDistribution(const blockDoubleJackknifeDistribution &r): baseType(r), bin_size(r.bin_size), nsample(r.nsample){}
  
  template<template<typename> class U>
  blockDoubleJackknifeDistribution(const blockDoubleJackknifeDistribution<BaseDataType,U> &r): blockDoubleJackknifeDistribution(r.nsample, r.bin_size){
    for(int i=0;i<this->size();i++)
      for(int j=0;j<this->sample(i).size();j++)
	this->sample(i).sample(j) = r.sample(i).sample(j);
  }

  explicit blockDoubleJackknifeDistribution(const initType &nb): blockDoubleJackknifeDistribution(nb.first, nb.second){} //used by ET
  
  blockDoubleJackknifeDistribution(const int nsample, const int bin_size): 
    baseType(nsample/bin_size, DataType(binCrop(nsample,bin_size)-bin_size)),
    nsample(binCrop(nsample,bin_size)), bin_size(bin_size){}

  blockDoubleJackknifeDistribution(const int nsample, const int bin_size, const DataType &init): 
    baseType(nsample/bin_size,init), nsample(binCrop(nsample,bin_size)), bin_size(bin_size){
    assert(init.size() == nsample - bin_size);
  }

  blockDoubleJackknifeDistribution(const int nsample, const int bin_size, const BaseDataType &init): 
    blockDoubleJackknifeDistribution(nsample, bin_size, DataType(binCrop(nsample,bin_size)-bin_size, init) ){}

  //Initializer functional should recognize that the number of samples will be cropped to a multiple of the bin size
  template<typename Initializer>
  blockDoubleJackknifeDistribution(const int nsample, const int bin_size, const Initializer &init): 
    baseType(nsample/bin_size,init), nsample(binCrop(nsample,bin_size)), bin_size(bin_size){}

  blockDoubleJackknifeDistribution(blockDoubleJackknifeDistribution&& o) noexcept : baseType(std::forward<baseType>(o)), bin_size(o.bin_size), nsample(o.nsample){}

  blockDoubleJackknifeDistribution(const rawDataDistribution<BaseDataType> &raw, const int bin_size): blockDoubleJackknifeDistribution(raw.size(), bin_size){
    this->resample(raw, bin_size);
  }
  
  ENABLE_PARALLEL_GENERIC_ET(blockDoubleJackknifeDistribution, myType, blockDoubleJackknifeDistribution<BaseDataType> );

  blockDoubleJackknifeDistribution & operator=(const blockDoubleJackknifeDistribution &r){ 
    nsample = r.nsample; bin_size=r.bin_size; static_cast<baseType*>(this)->operator=(r); return *this; 
  }
  blockDoubleJackknifeDistribution & operator=(blockDoubleJackknifeDistribution &&r){ 
    nsample = r.nsample; bin_size=r.bin_size; static_cast<baseType*>(this)->operator=(std::move(r)); return *this; 
  }

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
  blockDoubleJackknifeDistribution<typename U::value_type,BaseVectorType> real() const{
    blockDoubleJackknifeDistribution<typename U::value_type,BaseVectorType> out(this->size());
#pragma omp parallel for
    for(int i=0;i<this->size();i++)
      for(int j=0;j<this->sample(i).size();j++)
	out.sample(i).sample(j) = this->sample(i).sample(j).real();
    return out;
  }

  inline bool operator==(const blockDoubleJackknifeDistribution<BaseDataType,BaseVectorType> &r) const{ 
    return this->nsample == r.nsample && this->bin_size == r.bin_size & this->baseType::operator==(r); 
  }
  inline bool operator!=(const blockDoubleJackknifeDistribution<BaseDataType,BaseVectorType> &r) const{ return !( *this == r ); }
};

template<typename T, template<typename> class V>
std::ostream & operator<<(std::ostream &os, const blockDoubleJackknifeDistribution<T,V> &d){
  typedef distributionPrint<blockDoubleJackknifeDistribution<T,V> > printClass;
  assert(printClass::printer() != NULL); printClass::printer()->print(os, d);
  return os;
}

//Override basic printStats to print the central values and errors of the sub-distributions
template<typename T, template<typename> class V>
struct printStats< blockDoubleJackknifeDistribution<T,V> >{
  inline static std::string centralValue(const blockDoubleJackknifeDistribution<T,V> &d){
    std::ostringstream os; os << "[";
    for(int s=0;s<d.size()-1;s++) os << d.sample(s).best() << ", ";
    os << d.sample(d.size()-1).best() << "]";
    return os.str();
  }
  inline static std::string error(const blockDoubleJackknifeDistribution<T,V> &d){ 
    std::ostringstream os; os << "[";
    for(int s=0;s<d.size()-1;s++) os << d.sample(s).standardError() << ", ";
    os << d.sample(d.size()-1).standardError() << "]";
    return os.str();
  }

};

template<typename T, template<typename> class V>
blockDoubleJackknifeDistribution<T,V> weightedAvg(const std::vector<blockDoubleJackknifeDistribution<T,V> const*> &v){
  int nsample = v[0]->size();
  blockDoubleJackknifeDistribution<T,V> wavg_dj(*v[0]);
  std::vector<jackknifeDistribution<T,V> const*> towavg_j(v.size());
  for(int s=0;s<nsample;s++){ //weighted avg each jackknife distribution, looping over outer index
    for(int i=0;i<v.size();i++) towavg_j[i] = &v[i]->sample(s);
    wavg_dj.sample(s) = SARLaC::weightedAvg(towavg_j);
  }
  return wavg_dj;
}



SARLAC_END_NAMESPACE

#endif
