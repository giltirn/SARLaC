#ifndef _GLOBAL_DATA_MAP_READ_PIPI_TO_SIGMA_H
#define _GLOBAL_DATA_MAP_READ_PIPI_TO_SIGMA_H

#include<config.h>
#include<utils/macros.h>

#include <fit_pipitosigma_gparity/read_data.h>
#include <fit_pipitosigma_gparity/resampled_correlator.h>

CPSFIT_START_NAMESPACE

//For pipi->sigma, the "sigma" data set was measured with sources on every timeslice whereas the "extended" data set was measured with tstep=8. Thus to ensure correctness we should
//1) Ensure the 4 sets AsymmOnly, AsymmCorr, SymmCorr and SymmOnly each contain only data with a consistent tstep
//       < Fortunately as it currently stands, this is true without intervention. However it should be kept in mind! >
//2) Perform the source average prior to performing the sample-AMA correction

//---> If we do end up with mixed sets, they should be broken up and the error-weighted average performed under the superjackknife

//Note that the data for tstep=8 includes measuring the disconnected part with this tstep. We can reduce the statistical error on these data by
//1) Subtracting from the raw data the disconnected part, leaving just the connected part. This is possible because we have the pipi and sigma self contractions also stored.
//2) Performing the source average of the connected part with tstep=8
//3) Reconstructing the disconnected part for all timeslices and then performing the source average of the disconnected part with tstep=1
//4) Adding the two parts together


struct PiPiToSigmaMapReadPolicy{
  std::vector<DataLocationInfo const*> dinfo_vec;
  subStringReplace fmt;
  
  PiPiToSigmaMapReadPolicy(const std::map<int, DataLocationInfo const*> &dinfo_map): dinfo_vec(dinfo_map.size()){
    int sample = 0;
    for(auto it=dinfo_map.begin(); it!= dinfo_map.end(); it++) dinfo_vec[sample++] = it->second;

    std::string ffmt = "traj_<CONF>_pipitosigma_sigmawdagmom<MOM_QUARK_SIGMA>_pionmom<MOM_PI>_v2";
    fmt.chunkString(ffmt, { subStringSpecify("<CONF>"), subStringSpecify("<MOM_QUARK_SIGMA>"), subStringSpecify("<MOM_PI>") });
  }

  int nsample() const{ return dinfo_vec.size(); }
  
  std::string filename(const int sample, const threeMomentum &mom_quark_sigma, const threeMomentum &mom_pi) const{
    std::ostringstream os;
    os << dinfo_vec[sample]->directory << '/';
    fmt.replace(os, { anyToStr(dinfo_vec[sample]->inner_config), momStr(mom_quark_sigma), momStr(mom_pi) } );
    return os.str();
  }
};
  
template<typename resampledCorrFuncType>
void getPiPiToSigma(resampledCorrFuncType &out, const int Lt, const int tsep_pipi, const int tstep_src, 
		    const std::map<DataTag, std::map<int, DataLocationInfo const*> > &data_info_map, const int full_ens_size){
  std::vector<DataTag> loop_tags = { DataTag::AsymmOnly, DataTag::AsymmCorr, DataTag::SymmCorr, DataTag::SymmOnly };

  std::map<DataTag, resampledCorrFuncType> sjack;

  for(auto dtag = loop_tags.begin(); dtag != loop_tags.end(); dtag++){
    const std::map<int, DataLocationInfo const*> &subens = data_info_map.find(*dtag)->second;
    int nsample_ens = subens.size();
    
    bubbleDataZ raw_pipi_bubble;
    getA1projectedSourcePiPiBubble(raw_pipi_bubble, Lt, tsep_pipi, subens);
    
    sigmaSelfContractionZ raw_sigma_bubble;
    readSigmaSelf(raw_sigma_bubble, Lt, subens);

    PiPiToSigmaMapReadPolicy rd(subens);
    rawCorrelationFunction raw_srcavg = readReconstructPiPiToSigmaWithDisconnAllTsrc(rd, Lt, tstep_src,raw_pipi_bubble, raw_sigma_bubble);

    superJackknifeResample(sjack[*dtag], raw_srcavg, subens, full_ens_size);
  }

  out = combineResampledDataSets(sjack, data_info_map);
}

template<typename resampledCorrFuncType, typename resampledSigmaBubbleType, typename resampledPiPiBubbleDataAllMomentaType>
void performPiPiToSigmaVacuumSubtraction(resampledCorrFuncType &out,
					 const resampledCorrFuncType &in, 
					 const resampledSigmaBubbleType &sigma_bubble, const resampledPiPiBubbleDataAllMomentaType &pipi_bubble,
					 const int Lt, const std::vector<threeMomentum> &pion_mom){

  auto pipi_srcbub_proj = A1projectSourcePiPiBubble(pipi_bubble, pion_mom);
  auto vacsub = computePiPiToSigmaVacSub(sigma_bubble, pipi_srcbub_proj);
  out = resampledCorrFuncType(Lt, [&](const int t){ return typename resampledCorrFuncType::ElementType(in.coord(t), in.value(t) - vacsub.value(t)); }); 
}


CPSFIT_END_NAMESPACE

#endif
