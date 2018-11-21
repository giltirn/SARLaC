#ifndef READ_SIGMA_BUBBLE_H_
#define READ_SIGMA_BUBBLE_H_

#include<config.h>
#include<utils/macros.h>

#include<pipi_common/read_data_sigma.h>

CPSFIT_START_NAMESPACE

struct SigmaSelfMapReadPolicy{
  std::vector<DataLocationInfo const*> dinfo_vec;
  subStringReplace fmt;
  
  SigmaSelfMapReadPolicy(const std::map<int, DataLocationInfo const*> &dinfo_map): dinfo_vec(dinfo_map.size()){
    int sample = 0;
    for(auto it=dinfo_map.begin(); it!= dinfo_map.end(); it++) dinfo_vec[sample++] = it->second;

    std::string ffmt = "traj_<CONF>_sigmaself_mom<PQUARK>_v2";
    fmt.chunkString(ffmt, { subStringSpecify("<CONF>"), subStringSpecify("<PQUARK>") });
  }

  int nsample() const{ return dinfo_vec.size(); }
  
  std::string filename(const int sample, const threeMomentum &pquark) const{
    std::ostringstream os;
    os << dinfo_vec[sample]->directory << '/';
    fmt.replace(os, { anyToStr(dinfo_vec[sample]->inner_config), momStr(pquark) } );
    return os.str();
  }
};

//expect sigmaSelfContraction or sigmaSelfContractionZ
template<typename sigmaSelfContractionType>
void readSigmaSelf(sigmaSelfContractionType &into,
		   const int Lt,
		   const std::map<int, DataLocationInfo const*> &data_info_map){
  SigmaSelfMapReadPolicy rp(data_info_map);
  readSigmaSelf(into, Lt, rp);
}

template<typename resampledBubbleDataType>
void superJackknifeResample(resampledBubbleDataType &out,
			    const sigmaSelfContraction &in,
			    const std::map<int, DataLocationInfo const*> &data_info_map, const int full_ens_size){
  
  out.setup(in.getLt());
  std::vector<int> idxmap = getDataOuterConfigMap(data_info_map, full_ens_size); //mapping of outer to sample index, -1 for entries not present
  for(int t=0;t<in.getLt();t++)
    superJackknifeResample(out(t), in(t), idxmap, data_info_map.size(), full_ens_size);
}


//expect sigmaSelfContractionJack or sigmaSelfContractionDoubleJack
template<typename resampledBubbleDataType>
void getSigmaBubble(resampledBubbleDataType &out,
		    const int Lt,
		    const std::map<DataTag, std::map<int, DataLocationInfo const*> > &data_info_map, const int full_ens_size){

  std::vector<DataTag> loop_tags = { DataTag::AsymmOnly, DataTag::AsymmCorr, DataTag::SymmCorr, DataTag::SymmOnly };

  std::map<DataTag, resampledBubbleDataType> sjack;

  for(auto dtag = loop_tags.begin(); dtag != loop_tags.end(); dtag++){
    const std::map<int, DataLocationInfo const*> &subens = data_info_map.find(*dtag)->second;

    sigmaSelfContraction raw;
    readSigmaSelf(raw, Lt, subens);

    superJackknifeResample(sjack[*dtag], raw, subens, full_ens_size);
  }

  out = combineResampledBubble(sjack, data_info_map);
}

CPSFIT_END_NAMESPACE

#endif
