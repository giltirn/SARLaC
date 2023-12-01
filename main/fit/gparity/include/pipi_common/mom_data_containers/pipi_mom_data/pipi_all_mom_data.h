#pragma once

#include<config.h>
#include<utils/macros.h>

#include "pipi_bubble_mom_data.h"
#include "pipi_figure_mom_data.h"
#include "pipi_figureV_construct.h"

CPSFIT_START_NAMESPACE

//Read the raw contraction data. No rotational-state projection is done, but we do avoid reading data that won't be needed, hence the projectors input
//bubbleDataAllMomentaType = bubbleDataAllMomenta or bubbleDataAllMomentaZ
template<typename FigureReadPolicy, typename BubbleReadPolicy, typename bubbleDataAllMomentaType>
void readRawPiPi2ptData(figureDataAllMomenta &raw_data, bubbleDataAllMomentaType &raw_bubble_data,
			const FigureReadPolicy &frp, const BubbleReadPolicy &brp_src, const BubbleReadPolicy &brp_snk,
			const int Lt, const int tstep_pipi, const int tsep_pipi, 
			const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk){
  const char figs[3] = {'C','D','R'};
  for(int f=0;f<3;f++){
    readPiPi2ptFigure(raw_data, figs[f],  Lt,  proj_src, proj_snk, frp);

    //Some of Daiqian's old data was measured on every source timeslice while the majority was measured every 8. To fix this discrepancy we explicitly zero the abnormal data  
    zeroUnmeasuredSourceTimeslices(raw_data, figs[f], tstep_pipi);
  }

  readPiPiBubble(raw_bubble_data, Lt, tsep_pipi, brp_src, brp_snk, proj_src, proj_snk);

  //Populate the V diagrams from the bubble data
  computePiPi2ptFigureV(raw_data, raw_bubble_data, tsep_pipi, proj_src, proj_snk);
}
//Call the above using the standard read policies
template<typename FigureFilenamePolicy, typename BubbleFilenamePolicy, typename bubbleDataAllMomentaType>
void readRawPiPi2ptData(figureDataAllMomenta &raw_data, bubbleDataAllMomentaType &raw_bubble_data,
			const FigureFilenamePolicy &ffn, const BubbleFilenamePolicy &bfn_src, const BubbleFilenamePolicy &bfn_snk,
			const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, 
			const int Lt, const int tstep_pipi,
			const int tsep_pipi, const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk){
  PiPiBubbleBasicReadPolicy<BubbleFilenamePolicy> brp_src(bfn_src, data_dir, traj_start, traj_inc, traj_lessthan, tsep_pipi);
  PiPiBubbleBasicReadPolicy<BubbleFilenamePolicy> brp_snk(bfn_snk, data_dir, traj_start, traj_inc, traj_lessthan, tsep_pipi);
  PiPiFigureBasicReadPolicy<FigureFilenamePolicy> frp(ffn, data_dir, traj_start, traj_inc, traj_lessthan, tsep_pipi);
  readRawPiPi2ptData(raw_data, raw_bubble_data, frp, brp_src, brp_snk, Lt, tstep_pipi, tsep_pipi, proj_src, proj_snk);
}


CPSFIT_END_NAMESPACE
