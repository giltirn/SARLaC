#ifndef _PIPI_MAIN_H__
#define _PIPI_MAIN_H__

//Combine the computation of the V diagram with A2 projection and source average to avoid large intermediate data storage
template<typename BubbleDataType>
auto computeVprojectSourceAvg(const BubbleDataType &raw_bubble_data, const int tsep_pipi, const PiPiProject &proj_src, const PiPiProject &proj_snk, const PiPiMomAllow &allow, const std::vector<threeMomentum> &pion_momenta)
  ->correlationFunction<double,typename std::decay<decltype(raw_bubble_data(Source,*((threeMomentum*)NULL))(0))>::type>{

  (std::cout << "Computing projected, src-averaged V diagrams with BubbleDataType = " << printType<BubbleDataType>() << " and " << omp_get_max_threads() << " threads\n").flush(); 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed projected, src-averaged V diagrams with BubbleType = ") + printType<BubbleDataType>() + " in %w s\n");

  const int nmom = pion_momenta.size();
  
  const int Lt = raw_bubble_data.getLt();
  const int Nsample = raw_bubble_data.getNsample();

  typedef typename std::decay<decltype(raw_bubble_data(Source,*((threeMomentum*)NULL))(0))>::type  DistributionType;
  correlationFunction<double,DistributionType> out(Lt,
					    [&](const int t)
					    {
					      return typename correlationFunction<double,DistributionType>::ElementType(double(t), DistributionType(Nsample,0.));
					    }
					    );
  int nthr = omp_get_max_threads();
  std::vector<correlationFunction<double,DistributionType> > thr_sum(nthr, out);
  
  std::vector<std::pair<int,int> > todo;
  double dummy; std::complex<double> zdummy;
  for(int psnk=0;psnk<nmom;psnk++)
    for(int psrc=0;psrc<nmom;psrc++)
      if(proj_src(zdummy,pion_momenta[psrc]) && proj_snk(zdummy,pion_momenta[psnk]) && allow(dummy,pion_momenta[psrc],pion_momenta[psnk])){
	todo.push_back(std::make_pair(psrc,psnk));
      }

#pragma omp parallel for
  for(int pp=0;pp<todo.size();pp++){
    int me = omp_get_thread_num();
    int psrc = todo[pp].first;
    int psnk = todo[pp].second;

    std::complex<double> csrc, csnk;
    double m;
    proj_src(csrc, pion_momenta[psrc]);
    proj_snk(csnk, pion_momenta[psnk]);
    allow(m,pion_momenta[psrc],pion_momenta[psnk]);

    double coeff = m*std::real(csrc * csnk);

    const auto &Bp1_snk = raw_bubble_data(Sink, pion_momenta[psnk] ); 
    const auto &Bp1_src  = raw_bubble_data(Source,  pion_momenta[psrc] );

    for(int tsrc=0;tsrc<Lt;tsrc++)
      for(int tsep=0;tsep<Lt;tsep++)
	thr_sum[me].value(tsep) = thr_sum[me].value(tsep) + coeff * Bp1_snk( (tsrc + tsep) % Lt ) * Bp1_src( tsrc );
  }
  for(int tsep=0;tsep<Lt;tsep++){
    for(int thr=0;thr<nthr;thr++)
      out.value(tsep) = out.value(tsep) + thr_sum[thr].value(tsep);
    out.value(tsep) = out.value(tsep)/double(Lt);
  }
  return out;  
}

bubbleDataDoubleJackAllMomenta binDoubleJackknifeResampleBubble(const bubbleDataAllMomenta &bubbles, const int bin_size){
  int Lt = bubbles.getLt();
  bubbleDataDoubleJackAllMomenta out(Lt, bubbles.getTsepPiPi(), bubbles.getNsample()/bin_size);

  for(auto it = bubbles.begin(); it != bubbles.end(); it++){
    const auto & key = it->first;
    const typename bubbleDataAllMomenta::ContainerType & raw = it->second;
    for(int t=0;t<Lt;t++)
      out(key)(t).resample(raw(t).bin(bin_size));
  }
  return out;
}

template<typename T>
inline correlationFunction<double,T> fold(const correlationFunction<double,T> &f, const int tsep_pipi){
  const int Lt = f.size();
  correlationFunction<double,T> out(Lt);
  const int Tref = Lt-2*tsep_pipi;
  for(int t=0;t<Lt;t++){
    out.coord(t) = f.coord(t);
    out.value(t) = ( f.value(t) + f.value( (Tref-t+Lt) % Lt ) )/2.;
  }
  return out;
}

//Read and combine/double-jack resample data from original files or a checkpoint of the entire data set
template<typename FigureFilenamePolicy, typename BubbleFilenamePolicy>
doubleJackCorrelationFunction generateData(const PiPiProject &proj_src, const PiPiProject &proj_snk, const PiPiMomAllow &allow, 
					   const int isospin, const Args &args, const CMDline &cmdline,
					   const FigureFilenamePolicy &ffn, const BubbleFilenamePolicy &bfn_src, const BubbleFilenamePolicy &bfn_snk, const std::string &extra_descr){

  //Read the data
  bubbleDataAllMomenta raw_bubble_data;
  rawCorrelationFunction pipi_raw;
  {
    figureDataAllMomenta raw_data;
    readRawData(raw_data, raw_bubble_data, args, cmdline, ffn, bfn_src, bfn_snk, proj_src, proj_snk, allow, extra_descr);
    checkpointRawData(raw_data, raw_bubble_data, args, cmdline, extra_descr);
    getRawPiPiCorrFunc(pipi_raw, raw_data, raw_bubble_data, proj_src, proj_snk, allow, isospin, args.pion_momenta, args.bin_size, extra_descr);
  }

  const int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc/args.bin_size;
  
  doubleJackCorrelationFunction pipi_dj(args.Lt,
					[&](const int t)
					{
					  typename doubleJackCorrelationFunction::ElementType out(t, doubleJackknifeDistributionD(nsample));
					  out.second.resample(pipi_raw.value(t).bin(args.bin_size));
					  return out;
					}
					);
  
  if(isospin == 0 && args.do_vacuum_subtraction){
    bubbleDataDoubleJackAllMomenta dj_bubble_data = binDoubleJackknifeResampleBubble(raw_bubble_data, args.bin_size);
    doubleJackCorrelationFunction A2_realavg_V_dj = computeVprojectSourceAvg(dj_bubble_data,args.tsep_pipi,proj_src,proj_snk,allow,args.pion_momenta);
    pipi_dj = pipi_dj - 3*A2_realavg_V_dj;
  }
  
  pipi_dj = fold(pipi_dj, args.tsep_pipi);
  return pipi_dj;
}

doubleJackCorrelationFunction generateData(const PiPiProject &proj_src, const PiPiProject &proj_snk, const PiPiMomAllow &allow, const int isospin, const Args &args, const CMDline &cmdline){
  if(args.total_mom.size() == 1 && args.total_mom[0] == threeMomentum({0,0,0})){
    readFigureStationaryPolicy ffn(cmdline.use_symmetric_quark_momenta);
    readBubbleStationaryPolicy bfn_src(cmdline.use_symmetric_quark_momenta,Source);
    readBubbleStationaryPolicy bfn_snk(cmdline.use_symmetric_quark_momenta,Sink);
    return generateData(proj_src, proj_snk, allow, isospin, args, cmdline, ffn, bfn_src, bfn_snk,"");
  }else{ //average over the different total momenta provided
    if(cmdline.use_symmetric_quark_momenta) error_exit(std::cout << "getData use_symmetric_quark_momenta not supported for p_tot != (0,0,0)\n");

    doubleJackCorrelationFunction out;
    for(int p=0;p<args.total_mom.size();p++){
      assert(args.total_mom[p] != threeMomentum({0,0,0}));
      readFigureTianleComovingPolicy ffn(args.total_mom[p]);
      readBubbleTianleComovingPolicy bfn_src(args.total_mom[p],Source);
      readBubbleTianleComovingPolicy bfn_snk(-args.total_mom[p],Sink);
      PiPiProjectAllowOnlyExistingPionMom proj_src_f(proj_src, args.total_mom[p], args.pion_momenta); //make sure the pion momenta are in the computed set
      PiPiProjectAllowOnlyExistingPionMom proj_snk_f(proj_snk, -args.total_mom[p], args.pion_momenta);

      if(p==0) out = generateData(proj_src_f, proj_snk_f, allow, isospin, args, cmdline, ffn, bfn_src, bfn_snk, std::string("ptot") + momStr(args.total_mom[p]) );
      else out = out + generateData(proj_src_f, proj_snk_f, allow, isospin, args, cmdline, ffn, bfn_src, bfn_snk, std::string("ptot") + momStr(args.total_mom[p]));
    }
    out = out * (1./args.total_mom.size());
    return std::move(out);
  }
}

doubleJackCorrelationFunction generateData(const Args &args, const CMDline &cmdline){
  PiPiProject *proj_src = getProjector(args.proj_src);
  PiPiProject *proj_snk = getProjector(args.proj_snk);
  PiPiMomAllow* allow = getMomPairFilter(args.allowed_mom);

  doubleJackCorrelationFunction out = generateData(*proj_src, *proj_snk, *allow, args.isospin, args,cmdline);
  delete proj_src; delete proj_snk; delete allow;

  return std::move(out);
}

//Load a checkpoint of precomputed resampled correlation function or create it from raw data
doubleJackCorrelationFunction getData(const Args &args, const CMDline &cmdline){
  doubleJackCorrelationFunction data;
  if(cmdline.load_combined_data){
#ifdef HAVE_HDF5
    HDF5reader reader(cmdline.load_combined_data_file);
    read(reader, data, "data");
#else
    error_exit("getData: Reading amplitude data requires HDF5\n");
#endif
  }else{
    data = generateData(args,cmdline);
  }

  if(cmdline.save_combined_data){
#ifdef HAVE_HDF5
    HDF5writer writer(cmdline.save_combined_data_file);
    write(writer, data, "data");
#else
    error_exit("getData: Saving amplitude data requires HDF5\n");
#endif
  }
  return std::move(data);
}



#endif
