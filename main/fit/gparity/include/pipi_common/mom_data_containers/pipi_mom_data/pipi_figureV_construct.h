#pragma once

#include "pipi_bubble_mom_data.h"

CPSFIT_START_NAMESPACE

//The V diagrams is computed offline frome the bubble data. We only compute for pion momenta that are going to be used in the rotational-state projection
template<typename DataAllMomentumType, typename BubbleDataType>
void computePiPi2ptFigureV(DataAllMomentumType &raw_data, const BubbleDataType &raw_bubble_data, const int tsep_pipi, const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk){
  (std::cout << "Computing V diagrams with BubbleDataType = " << printType<BubbleDataType>() << " and " << omp_get_max_threads() << " threads\n").flush(); 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed V diagrams with BubbleType = ") + printType<BubbleDataType>() + " in %w s\n");
				 
  const int Lt = raw_bubble_data.getLt();
  raw_data.setup(Lt);

  int npsrc = proj_src.nMomenta();
  int npsnk = proj_snk.nMomenta();

  for(int i=0;i<npsrc;i++)
    for(int j=0;j<npsnk;j++)
      raw_data('V',momComb(proj_snk.momentum(j),proj_src.momentum(i)));

#pragma omp parallel for
  for(int pp=0;pp<npsrc*npsnk;pp++){
    //psnkidx + npsnk * psrcidx
    
    threeMomentum psrc = proj_src.momentum(pp / npsnk);
    threeMomentum psnk = proj_snk.momentum(pp % npsnk);
    
    //Bsnk(tsrc + tsep + tsep_pipi, p1_snk) Bsrc(tsrc, p1_src)
    const auto &Bp1_snk = raw_bubble_data(Sink, psnk);
    const auto &Bp1_src  = raw_bubble_data(Source, psrc);
    
    auto &into = raw_data('V',momComb(psnk,psrc));
    for(int tsrc=0;tsrc<Lt;tsrc++)
      for(int tsep=0;tsep<Lt;tsep++)
	into(tsrc,tsep) = Bp1_snk( (tsrc + tsep) % Lt ) * Bp1_src( tsrc );
  }
}

inline void computePiPi2ptFigureV(figureDataAllMomenta &raw_data, const bubbleDataAllMomentaZ &raw_bubble_data, const int tsep_pipi, const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk){
  auto re = reIm(raw_bubble_data,0);
  computePiPi2ptFigureV(raw_data, re, tsep_pipi, proj_src, proj_snk);
}

CPSFIT_END_NAMESPACE
