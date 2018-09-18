#ifndef READ_SIGMA2PT_H
#define READ_SIGMA2PT_H

#include<config.h>
#include<utils/macros.h>

#include <pipi_common/read_data_sigma.h>
#include <fit_sigmasigma_gparity/resampled_correlator.h>

CPSFIT_START_NAMESPACE

struct Sigma2ptMapReadPolicy{
  std::vector<DataLocationInfo const*> dinfo_vec;
  subStringReplace fmt;
  
  Sigma2ptMapReadPolicy(const std::map<int, DataLocationInfo const*> &dinfo_map): dinfo_vec(dinfo_map.size()){
    int sample = 0;
    for(auto it=dinfo_map.begin(); it!= dinfo_map.end(); it++) dinfo_vec[sample++] = it->second;

    std::string ffmt = "traj_<CONF>_sigmacorr_mompsrc<PSRC_QUARK>psnk<PSNK_QUARK>_v2";
    fmt.chunkString(ffmt, { subStringSpecify("<CONF>"), subStringSpecify("<PSRC_QUARK>"), subStringSpecify("<PSNK_QUARK>") });
  }

  int nsample() const{ return dinfo_vec.size(); }
  
  std::string filename(const int sample, const threeMomentum &psrc_quark, const threeMomentum &psnk_quark) const{
    std::ostringstream os;
    os << dinfo_vec[sample]->directory << '/';
    fmt.replace(os, { anyToStr(dinfo_vec[sample]->inner_config), momStr(psrc_quark), momStr(psnk_quark) } );
    return os.str();
  }
};

template<typename resampledCorrFuncType>
void getSigma2pt(resampledCorrFuncType &out, const int Lt, const std::map<DataTag, std::map<int, DataLocationInfo const*> > &data_info_map, const int full_ens_size){
  std::vector<DataTag> loop_tags = { DataTag::SymmOnly, DataTag::AsymmOnly, DataTag::AsymmCorr, DataTag::SymmCorr };

  std::map<DataTag, resampledCorrFuncType> sjack;

  for(auto dtag = loop_tags.begin(); dtag != loop_tags.end(); dtag++){
    const std::map<int, DataLocationInfo const*> &subens = data_info_map.find(*dtag)->second;

    figureData raw;
    Sigma2ptMapReadPolicy rd(subens);
    readSigmaSigma(raw, Lt, rd);
  
    rawCorrelationFunction raw_srcavg = sourceAverage(raw);
    
    superJackknifeResample(sjack[*dtag], raw_srcavg, subens, full_ens_size);
  }

  out = combineResampledDataSets(sjack, data_info_map);
}

template<typename resampledCorrFuncType, typename resampledSigmaBubbleType>
void performSigma2ptVacuumSubtraction(resampledCorrFuncType &out,
				      const resampledCorrFuncType &in, 
				      const resampledSigmaBubbleType &sigma_bubble){
  auto vacsub = computeSigmaVacSub(sigma_bubble);
  out = resampledCorrFuncType(in.size(), [&](const int t){ return typename resampledCorrFuncType::ElementType(in.coord(t), in.value(t) - vacsub.value(t)); }); 
}
      
CPSFIT_END_NAMESPACE

#endif
