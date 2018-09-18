#ifndef READ_PIPI_BUBBLE_H_
#define READ_PIPI_BUBBLE_H_

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//PiPi bubble format has a run-tag dependent format
struct PiPiBubbleMapReadPolicy{
  std::vector<DataLocationInfo const*> dinfo_vec;
  int tsep_pipi;
  threeMomentum p_tot;
  SourceOrSink src_snk;
  
  PiPiBubbleMapReadPolicy(const int tsep_pipi, const threeMomentum &p_tot, const SourceOrSink src_snk,
			  const std::map<int, DataLocationInfo const*> &dinfo_map): tsep_pipi(tsep_pipi), p_tot(p_tot), src_snk(src_snk), dinfo_vec(dinfo_map.size()){    
    int sample = 0;
    for(auto it=dinfo_map.begin(); it!= dinfo_map.end(); it++) dinfo_vec[sample++] = it->second;
  }

  int nsample() const{ return dinfo_vec.size(); }
  
  std::string filename(const int sample, const threeMomentum &p_pi1) const{
    const static std::map<std::string, std::string> runtag_fmt_map = {  {"orig", "traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_mom<PB>" },
									{"correction", "traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_mom<PB>_symm" },
									{"extended", "traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_pi1mom<PB>_pi2mom<PA>_symm" } };
    bubbleFilenamePolicyGeneric fp(runtag_fmt_map.find(dinfo_vec[sample]->tag)->second,
				   p_tot, src_snk);

    return fp(dinfo_vec[sample]->directory, dinfo_vec[sample]->inner_config, p_pi1, tsep_pipi);
  }
};

template<typename bubbleDataAllMomentaType>
void readPiPiBubble(bubbleDataAllMomentaType &raw_data, const int tsep_pipi, const int Lt,
		const std::map<int, DataLocationInfo const*> &dinfo_map,
		const std::vector<threeMomentum> &pion_momenta,
		const PiPiCorrelatorSelector &corr_select){
  PiPiBubbleMapReadPolicy rp_src(tsep_pipi, {0,0,0}, Source, dinfo_map);
  PiPiBubbleMapReadPolicy rp_snk(tsep_pipi, {0,0,0}, Sink, dinfo_map);
  readPiPiBubble(raw_data, Lt, tsep_pipi, rp_src, rp_snk, pion_momenta, corr_select);
}

//For use in pipi->sigma
//expected ContainerType = bubbleData or bubbleDataZ
template<typename ContainerType>
inline void getA1projectedSourcePiPiBubble(ContainerType &out, const int Lt, const int tsep_pipi, const std::vector<threeMomentum> &pion_mom, const std::map<int, DataLocationInfo const*> &data_info_map){
  PiPiBubbleMapReadPolicy rp(tsep_pipi, {0,0,0}, Source, data_info_map);
  getA1projectedSourcePiPiBubble(out,Lt,tsep_pipi,pion_mom, rp);
}

template<typename resampledBubbleDataType>
void superJackknifeResample(resampledBubbleDataType &out,
			    const bubbleData &in,
			    const std::map<int, DataLocationInfo const*> &data_info_map, const int full_ens_size){
  
  out.setup(Source, in.getLt(), in.getTsepPiPi(), full_ens_size);
  std::vector<int> idxmap = getDataOuterConfigMap(data_info_map, full_ens_size); //mapping of outer to sample index, -1 for entries not present
  for(int t=0;t<in.getLt();t++)
    superJackknifeResample(out(t), in(t), idxmap, data_info_map.size(), full_ens_size);
}

// ,   
// 	 typename std::enable_if< 
// 	   std::is_same<resampledBubbleDataType, bubbleDataJack>::value || std::is_same<resampledBubbleDataType, bubbleDataDoubleJack>::value,
// 	   void
// 	 >::type


template<typename bubbleDataType>
struct bubbleDataAccessor{
  typedef typename std::decay<decltype( ( (bubbleDataType*)NULL )->at(0) )>::type DistributionType;
  inline static DistributionType & at(bubbleDataType &corr, const int t){ return corr(t); }
  inline static const DistributionType & at(const bubbleDataType &corr, const int t){ return corr(t); }
  inline static bubbleDataType errorWeightedAverage(const bubbleDataType &a, const bubbleDataType &b){
    bubbleDataType out(a);
    int Lt = a.getLt();
    for(int t=0;t<Lt;t++)
      out(t) = CPSfit::errorWeightedAverage(a(t),b(t));
    return out;
  }
  inline static size_t size(const bubbleDataType &corr){ return corr.getLt(); }
  inline static std::string printCoord(const bubbleDataType &corr, const int t){ return anyToStr(t); }
  inline static std::string printValue(const DistributionType &dist){ return combDistPrint<DistributionType>::print(dist); }
};

template<typename bubbleDataType>
bubbleDataType combineResampledBubble(const std::map<DataTag, bubbleDataType> &sjack,
				      const std::map<DataTag, std::map<int, DataLocationInfo const*> > &data_info_map){
  return combineResampledDataSetsGeneric<bubbleDataType, bubbleDataAccessor<bubbleDataType> >(sjack, data_info_map);
}


template<typename bubbleDataAllMomentaType>
struct bubbleDataAllMomentaAccessor{
  typedef typename std::decay<decltype( ( (bubbleDataAllMomentaType*)NULL )->begin()->second(0) )>::type DistributionType;
  inline static DistributionType & at(bubbleDataAllMomentaType &corr, const int i){ 
    //mapping t + Lt*pidx
    int t = i % corr.getLt();
    int pidx = i / corr.getLt();
    return std::next(corr.begin(), pidx)->second(t);
  }
  inline static const DistributionType & at(const bubbleDataAllMomentaType &corr, const int i){
    int t = i % corr.getLt();
    int pidx = i / corr.getLt();
    return std::next(corr.begin(), pidx)->second(t);
  }   
  inline static bubbleDataAllMomentaType errorWeightedAverage(const bubbleDataAllMomentaType &a, const bubbleDataAllMomentaType &b){
    assert(a.getNmomenta() == b.getNmomenta());
    
    bubbleDataAllMomentaType out(a);
    int Lt = a.getLt(); assert(b.getLt() == Lt); assert(out.getLt() == Lt);
    auto ait = a.begin();
    auto bit = b.begin();
    auto oit = out.begin();

    while(oit != out.end()){
      for(int t=0;t<Lt;t++)
	oit->second(t) = CPSfit::errorWeightedAverage(ait->second(t),bit->second(t));
      ++ait; ++bit; ++oit;
    }
    return out;
  }
  inline static size_t size(const bubbleDataAllMomentaType &corr){ return corr.getNmomenta()*corr.getLt(); }
  inline static std::string printCoord(const bubbleDataAllMomentaType &corr, const int i){ 
    int t = i % corr.getLt();
    int pidx = i / corr.getLt();
    std::ostringstream os;  os << "[p=" << std::next(corr.begin(),pidx)->first  << ", t=" << t << "]";
    return os.str();
  }
  inline static std::string printValue(const DistributionType &dist){ return combDistPrint<DistributionType>::print(dist); }
};

template<typename bubbleDataAllMomentaType>
bubbleDataAllMomentaType combineResampledBubbleAllMomenta(const std::map<DataTag, bubbleDataAllMomentaType> &sjack,
							  const std::map<DataTag, std::map<int, DataLocationInfo const*> > &data_info_map){
  return combineResampledDataSetsGeneric<bubbleDataAllMomentaType, bubbleDataAllMomentaAccessor<bubbleDataAllMomentaType> >(sjack, data_info_map);
}



template<typename resampledBubbleDataAllMomentaType>
void getPiPiBubble(resampledBubbleDataAllMomentaType &out,
		   const int Lt, const int tsep_pipi, const std::vector<threeMomentum> &pion_mom,
		   const std::map<DataTag, std::map<int, DataLocationInfo const*> > &data_info_map, const int full_ens_size,
		   const PiPiProjector proj_src = PiPiProjector::A1, const PiPiProjector proj_snk = PiPiProjector::A1){
  PiPiCorrelatorBasicSelector corr_select(proj_src, proj_snk,PiPiMomAllowed::All,{0,0,0});

  std::vector<DataTag> loop_tags = { DataTag::AsymmOnly, DataTag::AsymmCorr, DataTag::SymmCorr, DataTag::SymmOnly };

  std::map<DataTag, resampledBubbleDataAllMomentaType> sjack;
  for(auto dtag = loop_tags.begin(); dtag != loop_tags.end(); dtag++){
    const auto &subens = data_info_map.find(*dtag)->second;
    bubbleDataAllMomenta raw;
    readPiPiBubble(raw, tsep_pipi, Lt, subens, pion_mom, corr_select);

    sjack[*dtag].setup(Lt,tsep_pipi,full_ens_size);
    for(auto it = raw.begin(); it != raw.end(); it++)
      superJackknifeResample(sjack[*dtag](it->first), it->second, subens, full_ens_size);
  }
  
  out = combineResampledBubbleAllMomenta(sjack, data_info_map);
}


CPSFIT_END_NAMESPACE

#endif

