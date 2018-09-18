#ifndef _PIPI_RAW_DATA_H_
#define _PIPI_RAW_DATA_H_

#include<config.h>
#include<utils/macros.h>

#include "mom_data_containers.h"
#include "mom_project.h"

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
void computePiPi2ptFigureV(DataAllMomentumType &raw_data, const BubbleDataType &raw_bubble_data, const int tsep_pipi, const std::vector<threeMomentum> &pion_momenta, 
	      const PiPiCorrelatorSelector &corr_select){
  (std::cout << "Computing V diagrams with BubbleDataType = " << printType<BubbleDataType>() << " and " << omp_get_max_threads() << " threads\n").flush(); 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed V diagrams with BubbleType = ") + printType<BubbleDataType>() + " in %w s\n");
				 
  const int Lt = raw_bubble_data.getLt();
  const int Nsample = raw_bubble_data.getNsample();
  raw_data.setup(Lt,Nsample);

  const int nmom = pion_momenta.size();
  
  //Populate output
  double dummy;

  std::vector<std::pair<int,int> > todo;

  for(int psnk=0;psnk<nmom;psnk++)
    for(int psrc=0;psrc<nmom;psrc++)
      if(corr_select(dummy,pion_momenta[psrc],pion_momenta[psnk])){
	auto &into = raw_data('V',momComb(pion_momenta[psnk],pion_momenta[psrc]));
	todo.push_back(std::make_pair(psrc,psnk));
      }

#pragma omp parallel for
  for(int pp=0;pp<todo.size();pp++){
    int psrc = todo[pp].first;
    int psnk = todo[pp].second;
    
    //Bsnk(tsrc + tsep + tsep_pipi, p1_snk) Bsrc(tsrc, p1_src)
    const auto &Bp1_snk = raw_bubble_data(Sink, pion_momenta[psnk] );
    const auto &Bp1_src  = raw_bubble_data(Source,  pion_momenta[psrc] );
    
    auto &into = raw_data('V',momComb(pion_momenta[psnk],pion_momenta[psrc]));
    for(int tsrc=0;tsrc<Lt;tsrc++)
      for(int tsep=0;tsep<Lt;tsep++)
	into(tsrc,tsep) = Bp1_snk( (tsrc + tsep) % Lt ) * Bp1_src( tsrc );
  }
}

inline std::string checkpointFilename(const std::string &stub, const std::string &extra_descr){
  std::ostringstream filename;
  filename << stub;
  if(extra_descr.size() > 0) filename << "_" << extra_descr;
  filename << ".hdf5";
  return filename.str();
}

void checkpointRawData(const figureDataAllMomenta &raw_data, const bubbleDataAllMomenta &raw_bubble_data, const std::string &filename_stub, const std::string &extra_descr){
  saveHDF5checkpoint(raw_data, raw_bubble_data, checkpointFilename(filename_stub, extra_descr) );
}

//Read the raw contraction data. No rotational-state projection is done, but we do avoid reading data that won't be needed, hence the discriminators input
template<typename FigureFilenamePolicy, typename BubbleFilenamePolicy>
void readRawData(figureDataAllMomenta &raw_data, bubbleDataAllMomenta &raw_bubble_data,
		 const FigureFilenamePolicy &ffn, const BubbleFilenamePolicy &bfn_src, const BubbleFilenamePolicy &bfn_snk,
		 const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, 
		 const int Lt, const int tstep_pipi,
		 const int tsep_pipi, const std::vector<threeMomentum> &pion_momenta, const PiPiCorrelatorSelector &corr_select,
		 const std::string &descr, bool load_hdf5_data_checkpoint = false, const std::string &load_hdf5_data_checkpoint_stub = ""){
  std::cout << "Reading raw data " << descr << std::endl;

  if(load_hdf5_data_checkpoint){
    loadHDF5checkpoint(raw_data, raw_bubble_data, checkpointFilename(load_hdf5_data_checkpoint_stub, descr));
  }else{
    readPiPi2ptFigure(raw_data, 'C', data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, ffn, pion_momenta, corr_select);
    readPiPi2ptFigure(raw_data, 'D', data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, ffn, pion_momenta, corr_select);
    readPiPi2ptFigure(raw_data, 'R', data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, ffn, pion_momenta, corr_select);
    readPiPiBubble(raw_bubble_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, bfn_src, bfn_snk, pion_momenta, corr_select);
  }
  //Do the stuff below even if reading from checkpoint because some older checkpoints were saved prior to these operations being performed

  //Populate the V diagrams from the bubble data
  computePiPi2ptFigureV(raw_data, raw_bubble_data, tsep_pipi, pion_momenta, corr_select);

  //Some of Daiqian's old data was measured on every source timeslice while the majority was measured every 8. To fix this discrepancy we explicitly zero the abnormal data  
  zeroUnmeasuredSourceTimeslices(raw_data, 'C', tstep_pipi);
  zeroUnmeasuredSourceTimeslices(raw_data, 'D', tstep_pipi);
  zeroUnmeasuredSourceTimeslices(raw_data, 'R', tstep_pipi);
}

CPSFIT_END_NAMESPACE

#endif
