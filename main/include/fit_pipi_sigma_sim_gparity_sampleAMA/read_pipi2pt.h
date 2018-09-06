#ifndef READ_PIPI_2PT_H_
#define READ_PIPI_2PT_H_

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//PiPi figure format has a run-tag dependent format, plus the extended requires a momentum combination mapping
struct PiPiFigureMapReadPolicy{
  std::vector<DataLocationInfo const*> dinfo_vec;
  int tsep_pipi;
  threeMomentum p_tot;
  std::vector<threeMomentum> pion_mom;

  PiPiFigureMapReadPolicy(const int tsep_pipi, const std::vector<threeMomentum> &pion_mom, const threeMomentum &p_tot,
			  const std::map<int, DataLocationInfo const*> &dinfo_map): tsep_pipi(tsep_pipi), pion_mom(pion_mom), p_tot(p_tot), dinfo_vec(dinfo_map.size()){    
    int sample = 0;
    for(auto it=dinfo_map.begin(); it!= dinfo_map.end(); it++) dinfo_vec[sample++] = it->second;
  }

  int nsample() const{ return dinfo_vec.size(); }
  
  std::string filename(const int sample, const char fig, const threeMomentum &psnk, const threeMomentum &psrc) const{
    const static std::map<std::string, std::string> runtag_fmt_map = {  {"orig", "traj_<TRAJ>_Figure<FIG>_sep<TSEP_PIPI>_mom<P1SRC>_mom<P1SNK>" },
									{"correction", "traj_<TRAJ>_Figure<FIG>_sep<TSEP_PIPI>_mom<P1SRC>_mom<P1SNK>_symm" },
									{"extended", "traj_<TRAJ>_Figure<FIG>_sep<TSEP_PIPI>_p1src<P1SRC>_p2src<P2SRC>_p1snk<P1SNK>_p2snk<P2SNK>_symm" } };

    const std::string &tag = dinfo_vec[sample]->tag;
    if(tag == "extended"){
      //Momentum pair mapping for extended data set
      static std::unique_ptr<PiPiSymmetrySubsetFigureFileMapping> ffn;
      static bool ffn_setup = false;
      if(!ffn_setup){ //mapping stays the same, but we need a sample config and directory to find available data
	ffn.reset(new PiPiSymmetrySubsetFigureFileMapping(dinfo_vec[sample]->directory, runtag_fmt_map.find(tag)->second,
							  dinfo_vec[sample]->inner_config, tsep_pipi, pion_mom, p_tot));
	ffn_setup = true;
      }
      return (*ffn)(dinfo_vec[sample]->directory, fig, dinfo_vec[sample]->inner_config, psnk, psrc, tsep_pipi);
    }else{
      figureFilenamePolicyGeneric fp(runtag_fmt_map.find(tag)->second, p_tot);
      return fp(dinfo_vec[sample]->directory, fig, dinfo_vec[sample]->inner_config, psnk, psrc, tsep_pipi);
    }
  }
};


void readPiPi2pt(rawCorrelationFunction &pipi_raw,
		 const int tsep_pipi, const int tstep_pipi, const int Lt,
		 const std::map<int, DataLocationInfo const*> &dinfo_map,
		 const std::vector<threeMomentum> &pion_mom,
		 const PiPiProjector proj_src = PiPiProjector::A1, const PiPiProjector proj_snk = PiPiProjector::A1){
  PiPiCorrelatorBasicSelector corr_select(proj_src, proj_snk,PiPiMomAllowed::All,{0,0,0});
  
  PiPiFigureMapReadPolicy rp(tsep_pipi, pion_mom, {0,0,0}, dinfo_map);

  //Read C, D, R diagrams
  figureDataAllMomenta raw_data;
  char figs[3] = {'C','D','R'};
  
  for(int f=0;f<3;f++){
    readFigure(raw_data, figs[f], Lt, pion_mom, corr_select, rp);
    zeroUnmeasuredSourceTimeslices(raw_data, figs[f], tstep_pipi);
  }

  bubbleDataAllMomenta raw_bubble_data;
  readBubble(raw_bubble_data, tsep_pipi, Lt, dinfo_map, pion_mom, corr_select);
  
  computeV(raw_data, raw_bubble_data, tsep_pipi, pion_mom, corr_select);

  //Combine diagrams to construct raw correlator
  const int bin_size = 1;
  const int isospin = 0;
  getRawPiPiCorrFunc(pipi_raw, raw_data, corr_select, isospin, pion_mom, bin_size); //binned, source-averaged pipi 2pt data
}

template<typename resampledCorrFuncType>
void getPiPi2pt(resampledCorrFuncType &out,
		const int Lt, const int tsep_pipi, const int tstep_pipi, 
		const std::map<DataTag, std::map<int, DataLocationInfo const*> > &data_info_map, const int full_ens_size,
		const std::vector<threeMomentum> &pion_mom,
		const PiPiProjector proj_src = PiPiProjector::A1, const PiPiProjector proj_snk = PiPiProjector::A1){
  std::vector<DataTag> loop_tags = { DataTag::AsymmOnly, DataTag::AsymmCorr, DataTag::SymmCorr, DataTag::SymmOnly };

  std::map<DataTag, resampledCorrFuncType> sjack;

  for(auto dtag = loop_tags.begin(); dtag != loop_tags.end(); dtag++){
    const std::map<int, DataLocationInfo const*> &subens = data_info_map.find(*dtag)->second;

    rawCorrelationFunction pipi_raw;
    readPiPi2pt(pipi_raw, tsep_pipi, tstep_pipi, Lt, subens, pion_mom, proj_src, proj_snk);

    superJackknifeResample(sjack[*dtag], pipi_raw, subens, full_ens_size);
  }

  out = combineResampledDataSets(sjack, data_info_map);
}

template<typename resampledCorrFuncType, typename resampledBubbleDataAllMomentaType>
void performPiPi2ptVacuumSubtraction(resampledCorrFuncType &out,
				     const resampledCorrFuncType &in, const resampledBubbleDataAllMomentaType &bubble,
				     const int Lt, const int tsep_pipi, const std::vector<threeMomentum> &pion_mom,
				     const PiPiProjector proj_src = PiPiProjector::A1, const PiPiProjector proj_snk = PiPiProjector::A1){
  PiPiCorrelatorBasicSelector corr_select(proj_src, proj_snk,PiPiMomAllowed::All,{0,0,0});
  auto vacsub = computeVprojectSourceAvg(bubble, tsep_pipi, corr_select, pion_mom);
  out = resampledCorrFuncType(Lt, [&](const int t){ return typename resampledCorrFuncType::ElementType(in.coord(t), in.value(t) - 3.*vacsub.value(t)); });  //3V
}

CPSFIT_END_NAMESPACE

#endif
