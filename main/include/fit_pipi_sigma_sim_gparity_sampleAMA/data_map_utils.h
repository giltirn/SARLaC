#ifndef _GLOBAL_DATA_MAP_UTILS_H
#define _GLOBAL_DATA_MAP_UTILS_H

#include<config.h>
#include<utils/macros.h>

#include "data_map.h"
#include "enums.h"

CPSFIT_START_NAMESPACE

//Get the mapping of an outer sample index to the sample index of a distribution measured on a subset
//Use -1 for entries that do not exist in the sub ensemble
std::vector<int> getSubEnsIdxMap(const std::set<int> &subsens_outer_cfgs, const int full_size){
  std::vector<int> out(full_size, -1);
  int s=0;
  for(auto it=subsens_outer_cfgs.begin(); it != subsens_outer_cfgs.end(); ++it) out[*it] = s++;
  return out;
}

std::vector<int> getDataOuterConfigMap(const std::map<int, DataLocationInfo const*> &data_loc_map, const int full_size){
  std::vector<int> out(full_size, -1);
  int s=0;
  for(auto it=data_loc_map.begin(); it != data_loc_map.end(); ++it) out[it->first] = s++;
  return out;
}

inline int getNsamplesWithTag(const DataTag tag, const std::map<DataTag, std::map<int, DataLocationInfo const*> > &data_info_map){
  auto it = data_info_map.find(tag);
  assert( it != data_info_map.end() );
  return it->second.size();
}


template<typename T>
inline void superJackknifeResample(jackknifeDistribution<T> &out, const rawDataDistribution<T> &data, const std::vector<int> &idxmap, const int subens_size, const int full_size){
  out = superJackknifeResampleGen(data, idxmap, subens_size, full_size);
}
template<typename T>
inline void superJackknifeResample(doubleJackknifeDistribution<T> &out, const rawDataDistribution<T> &data, const std::vector<int> &idxmap, const int subens_size, const int full_size){
  out = superDoubleJackknifeResampleGen(data, idxmap, subens_size, full_size);
}


template<typename resampledCorrFuncType>
void superJackknifeResample(resampledCorrFuncType &out,
			    const rawCorrelationFunction &in, 
			    const std::map<int, DataLocationInfo const*> &data_info_map, const int full_ens_size){
  typedef typename resampledCorrFuncType::DataType resampledDistType;
  out.resize(in.size());
  std::vector<int> idxmap = getDataOuterConfigMap(data_info_map, full_ens_size); //mapping of outer to sample index, -1 for entries not present
  for(int t=0;t<in.size();t++){
    out.coord(t) = in.coord(t);
    superJackknifeResample(out.value(t), in.value(t), idxmap, data_info_map.size(), full_ens_size);
  }
}


CPSFIT_END_NAMESPACE

#endif



// template<typename DistributionType>
// struct getElemDataType{
//   typedef std::decay<decltype(  iterate<DistributionType>::at(0, *( (DistributionType*)NULL ) )   )>::type type;
// };

// template<typename resampledDistributionType>
// struct superJackknifeResampler{};

// template<typename T>
// struct superJackknifeResampler<jackknifeDistribution<T> >{
//   inline static jackknifeDistribution<T> resample(const rawDataDistribution<T> &data, const std::vector<int> &idxmap, const int subens_size, const int full_size){
//     return superJackknifeResampleGen(data, idxmap, subens_size, full_size);
//   }
// };

// template<typename T>
// struct superJackknifeResampler<doubleJackknifeDistribution<T> >{
//   inline static doubleJackknifeDistribution<T> resample(const rawDataDistribution<T> &data, const std::vector<int> &idxmap, const int subens_size, const int full_size){
//     return superDoubleJackknifeResampleGen(data, idxmap, subens_size, full_size);
//   }
// };
