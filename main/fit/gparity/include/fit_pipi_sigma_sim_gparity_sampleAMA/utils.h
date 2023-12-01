#ifndef FIT_SAMPLEAMA_UTILS_H
#define FIT_SAMPLEAMA_UTILS_H

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

jackknifeDistributionD errorWeightedAverage(const jackknifeDistributionD &a, const jackknifeDistributionD &b){
  double erra = a.standardError();
  double wa = 1/erra/erra;

  double errb = b.standardError();
  double wb = 1/errb/errb;

  return (wa * a + wb * b) * (1./(wa + wb));
} 

//Double-jackknife is essentially a regular jackknife performed on a reduced ensemble where we have thrown away a sample i. Thus the error-weighted average
//should use the standard errors for each outer sample individually 
inline doubleJackknifeDistributionD errorWeightedAverage(const doubleJackknifeDistributionD &a, const doubleJackknifeDistributionD &b){
  doubleJackknifeDistributionD out(a.size());
#pragma omp parallel for
  for(int i=0;i<a.size();i++)
    out.sample(i) = errorWeightedAverage(a.sample(i),b.sample(i));
  return out;
}

template<typename DistributionType>
correlationFunction<double,DistributionType> errorWeightedAverage(const correlationFunction<double,DistributionType> &a, const correlationFunction<double,DistributionType> &b){
  assert(a.size() == b.size());
  correlationFunction<double,DistributionType> out(a.size());
  for(int i=0;i<a.size();i++){
    assert(a.coord(i) == b.coord(i));
    out.coord(i) = a.coord(i);
    out.value(i) = errorWeightedAverage(a.value(i), b.value(i));
  }
  return out;
}

inline std::set<int> truncate(const std::set<int> &iset, const int size){
  return std::set<int>( iset.begin(), std::next(iset.begin(),size) );
}

/*Accessor should have
  static DistributionType & at(CorrelationFunctionType &corr, const int elem)
  static const DistributionType & at(const CorrelationFunctionType &corr, const int elem)
  static CorrelationFunctionType errorWeightedAverage(const CorrelationFunctionType &a, const CorrelationFunctionType &b)
  static size_t size(const CorrelationFunctionType &corr)  //number of elements used in conjunction with accessor
  static std::string printCoord(const CorrelationFunctionType &corr, const int elem)
  static std::string printValue(const DistributionType &dist)
*/


template<typename CorrelationFunctionType, typename Accessor>
CorrelationFunctionType combineResampledDataSetsGeneric(const std::map<DataTag, CorrelationFunctionType> &sjack,
							const std::map<DataTag, std::map<int, DataLocationInfo const*> > &data_info_map){
  typedef typename std::decay<decltype(Accessor::at(  *( (CorrelationFunctionType*)(NULL) ), 0 ))>::type DistributionType;

  size_t size  = Accessor::size(sjack.begin()->second); 
  for(auto it=sjack.begin(); it != sjack.end(); it++) assert(Accessor::size(it->second) == size);

  //3 possible allowed scenarios:
  //1) Have asymm + correction
  //2) Have symm
  //3) Have (asymm + correction) + symm, if so do err-weighted avg of the two

  bool have_asymm = sjack.find(DataTag::AsymmOnly) != sjack.end();
  bool have_symm = sjack.find(DataTag::SymmOnly) != sjack.end();
  bool have_correction = (sjack.find(DataTag::AsymmCorr) != sjack.end()) && (sjack.find(DataTag::SymmCorr) != sjack.end());

  have_asymm = have_asymm && getNsamplesWithTag(DataTag::AsymmOnly, data_info_map) > 0;
  have_symm = have_symm && getNsamplesWithTag(DataTag::SymmOnly, data_info_map) > 0;
  have_correction = have_correction && getNsamplesWithTag(DataTag::AsymmCorr, data_info_map) > 0  &&   getNsamplesWithTag(DataTag::SymmCorr, data_info_map) > 0;
  
  bool have_asymm_and_corr = have_asymm && have_correction;

#define E(c,t) Accessor::at(c,t)
#define C(c,t) Accessor::printCoord(c,t)
#define P(c,t) Accessor::printValue(Accessor::at(c,t))
#define PD(d) Accessor::printValue(d)

  CorrelationFunctionType asymm_corrected;
  if(have_asymm_and_corr){
    const CorrelationFunctionType &asymm_only = sjack.find(DataTag::AsymmOnly)->second;
    const CorrelationFunctionType &symm_corr = sjack.find(DataTag::SymmCorr)->second;
    const CorrelationFunctionType &asymm_corr = sjack.find(DataTag::AsymmCorr)->second;

    asymm_corrected = asymm_only;  
    std::cout << "combineResampledDataSetsGeneric have asymm and correction configs - performing sampleAMA correction\n";
    std::cout << "sampleAMA correction:\n";
    for(int i=0;i<size;i++){
      E(asymm_corrected, i) = E(asymm_corrected, i) + E(symm_corr, i) - E(asymm_corr, i);
    
      DistributionType reldiff = 2. * (E(asymm_corrected, i) - E(asymm_only, i)) / (E(asymm_corrected, i) + E(asymm_only, i));
      
      std::cout << C(asymm_only, i) << " " << P(asymm_only, i) << " " << P(asymm_corrected, i) << " " << PD(reldiff) << std::endl;
    }
  }

  if(have_asymm_and_corr && have_symm){
    std::cout << "combineResampledDataSetsGeneric have sampleAMA corrected data *and* symmetric data - computing error weighted average of data sets\n";
    const CorrelationFunctionType &symm_only = sjack.find(DataTag::SymmOnly)->second;
    CorrelationFunctionType errw = Accessor::errorWeightedAverage(asymm_corrected, symm_only);
    std::cout << "Error weighted avg of corrected and symm-only:\n";
    for(int i=0;i<size;i++){
      std::cout << C(errw,i) << " " << P(asymm_corrected, i) << " " << P(symm_only, i) << " " << P(errw, i) << std::endl;
    }
    return errw;    
  }else if(have_asymm_and_corr){
    return asymm_corrected;
  }else if(have_symm){
    std::cout << "combineResampledDataSetsGeneric have symmetric data only:\n";
    const CorrelationFunctionType &symm_only = sjack.find(DataTag::SymmOnly)->second;
    for(int i=0;i<size;i++){
      std::cout << C(symm_only,i) << " " << P(symm_only, i) << std::endl;
    }
    return symm_only;
  }else error_exit(std::cout << "combineResampledDataSets have neither symmetric or corrected-asymmetric data!\n"); 

#undef E
#undef C
#undef P
#undef PD
}

template<typename DistributionType>
struct combDistPrint{};

template<>
struct combDistPrint<jackknifeDistributionD>{
  static inline std::string print(const jackknifeDistributionD &d){ std::ostringstream os; os << d; return os.str(); }
};
template<>
struct combDistPrint<doubleJackknifeDistributionD>{ //use outer sample 0 as representative. Could also take mean to convert to single jackknife but at cost of flops
  static inline std::string print(const doubleJackknifeDistributionD &d){ std::ostringstream os; os << "[DJ i=0] " << d.sample(0); return os.str(); }
};

template<typename DistributionType>
struct correlationFunctionAccessor{
  typedef correlationFunction<double,DistributionType> CorrFunc;
  inline static DistributionType & at(CorrFunc &corr, const int t){ return corr.value(t); }
  inline static const DistributionType & at(const CorrFunc &corr, const int t){ return corr.value(t); }
  inline static CorrFunc errorWeightedAverage(const CorrFunc &a, const CorrFunc &b){ return SARLaC::errorWeightedAverage(a,b); }
  inline static size_t size(const CorrFunc &corr){ return corr.size(); }
  inline static std::string printCoord(const CorrFunc &corr, const int t){ return anyToStr(t); }
  inline static std::string printValue(const DistributionType &dist){ return combDistPrint<DistributionType>::print(dist); }
};

template<typename DistributionType>
correlationFunction<double,DistributionType> combineResampledDataSets(const std::map<DataTag, correlationFunction<double,DistributionType> > &sjack,
								      const std::map<DataTag, std::map<int, DataLocationInfo const*> > &data_info_map){
  return combineResampledDataSetsGeneric<correlationFunction<double,DistributionType>, correlationFunctionAccessor<DistributionType> >(sjack, data_info_map);
}

SARLAC_END_NAMESPACE

#endif


