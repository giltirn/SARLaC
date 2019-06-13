#ifndef _BOOTSTRAP_CLASS_H_
#define _BOOTSTRAP_CLASS_H_

#include<distribution/raw_data_distribution.h>
#include<random/rng.h>
#include<utils/utils/text_output.h>
//A distribution for bootstrap data

CPSFIT_START_NAMESPACE

//To ensure we retain correlations when bootstrapping multiple raw distributions measured on the same ensemble we must use the same mapping of bootstrap sample
//to the set of raw sample indices each time. This is accomplished by using a fixed seed for an RNG that is re-initialized each time. 
//User should change the seed between ensembles to remove fictional correlations

struct bootstrapDistributionOptions{
  static int & defaultBoots(){ static int b = 1000; return b; }
  static int & defaultSeed(){ static int s = 1234; return s; }
  static int & defaultConfidence(){ static int c=68; return c; }
};

template<typename _DataType, template<typename> class _VectorType = basic_vector>
class bootstrapDistribution: public distribution<_DataType,_VectorType>{
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

  inline void zero(){
    zeroit(avg); this->baseType::zero();
  }

  //Compute the m'th standardized moment of the distribution
  DataType standardizedMoment(const int m) const{
    DataType sigma = this->standardDeviation();
    DataType mean = this->mean();
  
    const int N = this->size() == 0 ? 1 : this->size();
    const int nthread = omp_get_max_threads();

    struct Op{
      DataType avg;
      const _VectorType<DataType> &data;
      int m;
      Op(const DataType&_avg, const _VectorType<DataType> &_data, const int m): data(_data), avg(_avg), m(m){}
      inline int size() const{ return data.size(); }
      inline DataType operator()(const int i) const{ DataType tmp = data[i] - avg; return pow(tmp,m);  }
    };
    Op op(mean,this->sampleVector(),m);
    
    return threadedSum(op)/double(N)/pow(sigma,m);
  }

  static std::vector<std::vector<int> > nonoverlappingBlockResampleTable(RNGstore &brng, const int nsample, const int block_size,
									 const int nboots = bootstrapDistributionOptions::defaultBoots()){
    int nblock = nsample / block_size;
    int nsample_b = nblock * block_size; //just in case block size doesn't divide the number of samples equally
    
    std::uniform_int_distribution<> dis(0,nblock-1);
    std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample_b));
    for(int b=0;b<nboots;b++){
      for(int i=0;i<nblock;i++){
  	int block_idx = dis(brng());

  	int to_off = block_size*i;
  	int from_off = block_size*block_idx;

  	for(int ss=0;ss<block_size;ss++)
  	  out[b][to_off + ss] = from_off + ss;
      }
    }
    return out;
  }

  inline static std::vector<std::vector<int> > nonoverlappingBlockResampleTable(const int nsample, const int block_size,
									 const int nboots = bootstrapDistributionOptions::defaultBoots(),
									 const int seed = bootstrapDistributionOptions::defaultSeed()){
    RNGstore brng(seed);
    return nonoverlappingBlockResampleTable(brng, nsample, block_size, nboots);
  }

  static std::vector<std::vector<int> > overlappingBlockResampleTable(RNGstore &brng, const int nsample, const int block_size,
								      const int nboots = bootstrapDistributionOptions::defaultBoots()){
    //Here we use overlapping blocks of size block_size, of which there are nsample-block_size+1. These are the blocks that are resampled
    //We draw floor(n/block_size) of these with replacement and set them down in order, which gives n samples back
    int nblock_ov = nsample-block_size+1;
    int nresample = nsample/block_size;
    int nsample_b = nresample * block_size; //just in case nsample is not evenly divisible by block size

    std::uniform_int_distribution<> dis(0,nblock_ov-1);

    std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample_b));    
    for(int b=0;b<nboots;b++){
      for(int i=0;i<nresample;i++){
	int block_idx = dis(brng());

	int to_off = block_size*i;
	int from_off = block_idx;

	for(int ss=0;ss<block_size;ss++)
	  out[b][to_off + ss] = from_off + ss;
      }
    }
    return out;
  }

  inline static std::vector<std::vector<int> > overlappingBlockResampleTable(const int nsample, const int block_size,
									     const int nboots = bootstrapDistributionOptions::defaultBoots(), 
									     const int seed = bootstrapDistributionOptions::defaultSeed()){
    RNGstore brng(seed);
    return overlappingBlockResampleTable(brng, nsample, block_size, nboots);
  }


  //One problem with the overlapping block bootstrap is that the underlying samples do not appear in the blocks with equal frequency. For example, sample 0 only appears in block 0, whereas sample 1 appears in 0 and 1, sample 2 in 0,1,2 and so on up to the block size. This can apparently skew the mean away from that of the underlying ensemble. In the circular version the configurations are laid out on a circle
  static std::vector<std::vector<int> > circularOverlappingBlockResampleTable(RNGstore &brng, const int nsample, const int block_size,
								      const int nboots = bootstrapDistributionOptions::defaultBoots()){
    //Here we use circularly overlapping blocks of size block_size, of which there are nsample. These are the blocks that are resampled
    //We draw floor(n/block_size) of these with replacement and set them down in order, which gives n samples back
    int nblock_ov = nsample;
    int nresample = nsample/block_size;
    int nsample_b = nresample * block_size; //just in case nsample is not evenly divisible by block size

    std::uniform_int_distribution<> dis(0,nblock_ov-1);

    std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample_b));    
    for(int b=0;b<nboots;b++){
      for(int i=0;i<nresample;i++){
	int block_idx = dis(brng());

	int to_off = block_size*i;
	int from_off = block_idx;

	for(int ss=0;ss<block_size;ss++)
	  out[b][to_off + ss] = (from_off + ss) % nsample;
      }
    }
    return out;
  }

  inline static std::vector<std::vector<int> > circularOverlappingBlockResampleTable(const int nsample, const int block_size,
									     const int nboots = bootstrapDistributionOptions::defaultBoots(), 
									     const int seed = bootstrapDistributionOptions::defaultSeed()){
    RNGstore brng(seed);
    return circularOverlappingBlockResampleTable(brng, nsample, block_size, nboots);
  }



  
  //Generate the mapping between the resampled ensembles and the original. Output indices are [boot][sample]
  static std::vector<std::vector<int> > resampleTable(RNGstore &brng, const int nsample, 
						      const int nboots = bootstrapDistributionOptions::defaultBoots()){
    std::uniform_int_distribution<> dis(0,nsample-1);
    std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample));
    for(int b=0;b<nboots;b++)
      for(int i=0;i<nsample;i++)
	out[b][i] = dis(brng());
    return out;
  }

  inline static std::vector<std::vector<int> > resampleTable(const int nsample, 
							     const int nboots = bootstrapDistributionOptions::defaultBoots(), 
							     const int seed = bootstrapDistributionOptions::defaultSeed()){
    RNGstore brng(seed);
    return resampleTable(brng, nsample, nboots);
  }

  template<typename DistributionType>
  void resample(const DistributionType &in, const std::vector<std::vector<int> > &table){
    assert(table.size() == this->size());
    int boots = this->size();
    assert(table[0].size() == in.size());
    int nraw = in.size();

    for(int b=0;b<boots;b++){
      zeroit(this->sample(b));

      //Draw N samples with replacement
      for(int i=0;i<nraw;i++)
	this->sample(b) = this->sample(b) + in.sample(table[b][i]);
      
      this->sample(b) = this->sample(b)/nraw;
    }
  }

  //This version generates the mapping on-the-fly. The same seed should be used for all data
  template<typename DistributionType> //doesn't have to be a distribution, just has to have a .sample and .size method
  void resample(const DistributionType &in, const int seed = bootstrapDistributionOptions::defaultSeed()){
    RNGstore brng(seed);
    
    this->avg = in.best();
    
    int nraw = in.size();
    std::uniform_int_distribution<> dis(0,nraw-1);

    int boots = this->size();

    for(int b=0;b<boots;b++){
      zeroit(this->sample(b));

      //Draw N samples with replacement
      for(int i=0;i<nraw;i++)
	this->sample(b) = this->sample(b) + in.sample(dis(brng()));
      
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
 
  bootstrapDistribution(const bootstrapDistribution &r) = default;
  
  template<template<typename> class U>
  bootstrapDistribution(const bootstrapDistribution<DataType,U> &r): baseType(r), _confidence(r._confidence), avg(r.avg){}
  
  struct initType: public OstreamHook{
    int boots;
    int confidence;
    initType(const int b = bootstrapDistributionOptions::defaultBoots(), const int c = bootstrapDistributionOptions::defaultConfidence()): boots(b), confidence(c){}
    inline bool operator==(const initType &r) const{ return boots==r.boots && confidence==r.confidence; }
    inline bool operator!=(const initType &r) const{ return !(*this == r); }
    void write(std::ostream &os) const{ os << "(boots=" << boots<< ", confidence=" << confidence <<")"; }
  };
  bootstrapDistribution(const initType &init = initType()): baseType(init.boots), _confidence(init.confidence){}

  bootstrapDistribution(const DataType &initv, const initType &init = initType()): baseType(init.boots,initv), _confidence(init.confidence){ avg = this->mean(); }

  template<typename Initializer>
  bootstrapDistribution(const Initializer &initf, const initType &init = initType()): baseType(init.boots,initf), _confidence(init.confidence){ avg = this->mean(); }

  bootstrapDistribution(bootstrapDistribution&& o) noexcept : baseType(std::forward<baseType>(o)), _confidence(o._confidence), avg(o.avg){}

  template< template<typename> class U >
  bootstrapDistribution(const rawDataDistribution<DataType,U> &raw, const initType &init = initType()): bootstrapDistribution(init){ this->resample(raw); }
  

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
