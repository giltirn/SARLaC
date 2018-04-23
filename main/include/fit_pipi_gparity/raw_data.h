#ifndef _PIPI_RAW_DATA_H_
#define _PIPI_RAW_DATA_H_

void zeroUnmeasuredSourceTimeslices(figureDataAllMomenta &data, const char fig, const int tstep_pipi){
  for(figureDataAllMomenta::iterator it = data.begin(fig); it != data.end(fig); it++){
    figureData &f = it->second;
    int Lt = f.getLt();
    for(int tsrc=0;tsrc<Lt;tsrc++)
      if(tsrc % tstep_pipi != 0)	  
	for(int tsep=0;tsep<Lt;tsep++) zeroit(f(tsrc,tsep));

  }    
}

template<typename DataAllMomentumType, typename BubbleDataType>
void computeV(DataAllMomentumType &raw_data, const BubbleDataType &raw_bubble_data, const int tsep_pipi, const std::vector<threeMomentum> &pion_momenta, 
	      const PiPiProject &proj_src, const PiPiProject &proj_snk, const PiPiMomAllow &allow){
  (std::cout << "Computing V diagrams with BubbleDataType = " << printType<BubbleDataType>() << " and " << omp_get_max_threads() << " threads\n").flush(); 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed V diagrams with BubbleType = ") + printType<BubbleDataType>() + " in %w s\n");
				 
  const int Lt = raw_bubble_data.getLt();
  const int Nsample = raw_bubble_data.getNsample();
  raw_data.setup(Lt,Nsample);

  const int nmom = pion_momenta.size();
  
  //Populate output
  double dummy; std::complex<double> zdummy;

  std::vector<std::pair<int,int> > todo;

  for(int psnk=0;psnk<nmom;psnk++)
    for(int psrc=0;psrc<nmom;psrc++)
      if(proj_src(zdummy,pion_momenta[psrc]) && proj_snk(zdummy,pion_momenta[psnk]) && allow(dummy,pion_momenta[psrc],pion_momenta[psnk])){
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

void checkpointRawData(const figureDataAllMomenta &raw_data, const bubbleDataAllMomenta &raw_bubble_data, const Args &args, const CMDline &cmdline, const std::string &extra_descr){
  if(cmdline.save_hdf5_data_checkpoint)
    saveHDF5checkpoint(raw_data, raw_bubble_data, checkpointFilename(cmdline.save_hdf5_data_checkpoint_stub, extra_descr) );
}

//Read the contraction data
template<typename FigureFilenamePolicy, typename BubbleFilenamePolicy>
void readRawData(figureDataAllMomenta &raw_data, bubbleDataAllMomenta &raw_bubble_data, const Args &args, const CMDline &cmdline, 
		 const FigureFilenamePolicy &ffn, const BubbleFilenamePolicy &bfn_src, const BubbleFilenamePolicy &bfn_snk,
		 const PiPiProject &proj_src, const PiPiProject &proj_snk, const PiPiMomAllow &allow,
		 const std::string &extra_descr = ""){
  std::cout << "Reading raw data " << extra_descr << std::endl;

  if(cmdline.load_hdf5_data_checkpoint){
    loadHDF5checkpoint(raw_data, raw_bubble_data, checkpointFilename(cmdline.load_hdf5_data_checkpoint_stub, extra_descr));
  }else{
    readFigure(raw_data, 'C', args.data_dir, args.tsep_pipi, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, ffn, args.pion_momenta, proj_src, proj_snk, allow);
    readFigure(raw_data, 'D', args.data_dir, args.tsep_pipi, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, ffn, args.pion_momenta, proj_src, proj_snk, allow);
    readFigure(raw_data, 'R', args.data_dir, args.tsep_pipi, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, ffn, args.pion_momenta, proj_src, proj_snk, allow);
    readBubble(raw_bubble_data, args.data_dir, args.tsep_pipi, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, bfn_src, bfn_snk, args.pion_momenta, proj_src, proj_snk, allow);
  }
  //Do the stuff below even if reading from checkpoint because some older checkpoints were saved prior to these operations being performed

  //Populate the V diagrams from the bubble data
  computeV(raw_data, raw_bubble_data, args.tsep_pipi, args.pion_momenta, proj_src, proj_snk, allow);

  //Some of Daiqian's old data was measured on every source timeslice while the majority was measured every 8. To fix this discrepancy we explicitly zero the abnormal data  
  zeroUnmeasuredSourceTimeslices(raw_data, 'C', args.tstep_pipi);
  zeroUnmeasuredSourceTimeslices(raw_data, 'D', args.tstep_pipi);
  zeroUnmeasuredSourceTimeslices(raw_data, 'R', args.tstep_pipi);
}



#endif
