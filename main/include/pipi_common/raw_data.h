#ifndef _PIPI_RAW_DATA_H_
#define _PIPI_RAW_DATA_H_

#include<config.h>
#include<utils/macros.h>

#include "mom_data_containers.h"
#include "mom_project.h"
#include "read_data_pipi.h"
CPSFIT_START_NAMESPACE

void zeroUnmeasuredSourceTimeslices(figureDataAllMomenta &data, const char fig, const int tstep_pipi){
  for(figureDataAllMomenta::iterator it = data.begin(fig); it != data.end(fig); it++){
    figureData &f = it->second;
    int Lt = f.getLt();
    for(int tsrc=0;tsrc<Lt;tsrc++)
      if(tsrc % tstep_pipi != 0)	  
	for(int tsep=0;tsep<Lt;tsep++) zeroit(f(tsrc,tsep));

  }    
}

//The V diagrams is computed offline frome the bubble data. We only compute for pion momenta that are going to be used in the rotational-state projection
template<typename DataAllMomentumType, typename BubbleDataType>
void computePiPi2ptFigureV(DataAllMomentumType &raw_data, const BubbleDataType &raw_bubble_data, const int tsep_pipi, const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk){
  (std::cout << "Computing V diagrams with BubbleDataType = " << printType<BubbleDataType>() << " and " << omp_get_max_threads() << " threads\n").flush(); 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed V diagrams with BubbleType = ") + printType<BubbleDataType>() + " in %w s\n");
				 
  const int Lt = raw_bubble_data.getLt();
  const int Nsample = raw_bubble_data.getNsample();
  raw_data.setup(Lt,Nsample);

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



inline std::string checkpointFilename(const std::string &stub, const std::string &extra_descr){
  std::ostringstream filename;
  filename << stub;
  if(extra_descr.size() > 0) filename << "_" << extra_descr;
  filename << ".hdf5";
  return filename.str();
}

void saveRawDataCheckpoint(const figureDataAllMomenta &raw_data, const bubbleDataAllMomenta &raw_bubble_data, const std::string &filename_stub, const std::string &extra_descr){
  saveHDF5checkpoint(raw_data, raw_bubble_data, checkpointFilename(filename_stub, extra_descr) );
}
void loadRawDataCheckpoint(figureDataAllMomenta &raw_data, bubbleDataAllMomenta &raw_bubble_data, const std::string &filename_stub, const std::string &extra_descr){
  loadHDF5checkpoint(raw_data, raw_bubble_data, checkpointFilename(filename_stub, extra_descr) );
}



//Read the raw contraction data. No rotational-state projection is done, but we do avoid reading data that won't be needed, hence the projectors input
//bubbleDataAllMomentaType = bubbleDataAllMomenta or bubbleDataAllMomentaZ
template<typename FigureReadPolicy, typename BubbleReadPolicy, typename bubbleDataAllMomentaType>
void readRawPiPi2ptData(figureDataAllMomenta &raw_data, bubbleDataAllMomentaType &raw_bubble_data,
			const FigureReadPolicy &frp, const BubbleReadPolicy &brp_src, const BubbleReadPolicy &brp_snk,
			const int Lt, const int tstep_pipi, const int tsep_pipi, 
			const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk){
  figureDataPolicies::useFileCache() = true;
  const char figs[3] = {'C','D','R'};
  for(int f=0;f<3;f++){
    readPiPi2ptFigure(raw_data, figs[f],  Lt,  proj_src, proj_snk, frp);

    //Some of Daiqian's old data was measured on every source timeslice while the majority was measured every 8. To fix this discrepancy we explicitly zero the abnormal data  
    zeroUnmeasuredSourceTimeslices(raw_data, figs[f], tstep_pipi);
    figureDataPolicies::getFileCache().clear();
  }
  figureDataPolicies::useFileCache() = false;

  readPiPiBubble(raw_bubble_data, Lt, tsep_pipi, brp_src, brp_snk, proj_src, proj_snk);

  //Populate the V diagrams from the bubble data
  computePiPi2ptFigureV(raw_data, raw_bubble_data, tsep_pipi, proj_src, proj_snk);
}

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




//Perform the rotational state projection
template<typename DataAllMomentumType>
typename DataAllMomentumType::ContainerType project(const char fig, const DataAllMomentumType &raw_data, 
						    const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk){
  std::cout << "Computing projection of figure " << fig << " with DataAllMomentumType = " << printType<DataAllMomentumType>() << "\n"; 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed projection of figure ") + fig + " with DataAllMomentumType = " + printType<DataAllMomentumType>() + " in %w s\n");

  typename DataAllMomentumType::ContainerType out(raw_data.getLt(), raw_data.getNsample()); out.zero();

  for(int psnki=0; psnki<proj_snk.nMomenta();psnki++){
    threeMomentum psnk = proj_snk.momentum(psnki);

    for(int psrci=0; psrci<proj_src.nMomenta();psrci++){
      threeMomentum psrc = proj_src.momentum(psrci);

      out = out + std::real(proj_snk.coefficient(psnki)*proj_src.coefficient(psrci)) * raw_data(fig, momComb(psnk,psrc));
    }
  }

  return out;
}



CPSFIT_END_NAMESPACE

#endif
