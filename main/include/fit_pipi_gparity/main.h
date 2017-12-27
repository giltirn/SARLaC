#ifndef _PIPI_MAIN_H__
#define _PIPI_MAIN_H__

void zeroUnmeasuredSourceTimeslices(figureDataAllMomenta &data, const char fig, const int tstep_pipi){
  for(figureDataAllMomenta::iterator it = data.begin(fig); it != data.end(fig); it++){
    figureData &f = it->second;
    int Lt = f.getLt();
    for(int tsrc=0;tsrc<Lt;tsrc++)
      if(tsrc % tstep_pipi != 0)	  
	for(int tsep=0;tsep<Lt;tsep++) zeroit(f(tsrc,tsep));

  }    
}

template<typename DataAllMomentumType>
typename DataAllMomentumType::ContainerType projectA2(const char fig, const DataAllMomentumType &raw_data){
  std::cout << "Computing A2 projection of figure " << fig << " with DataAllMomentumType = " << printType<DataAllMomentumType>() << "\n"; 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed A2 projection of figure ") + fig + " with DataAllMomentumType = " + printType<DataAllMomentumType>() + " in %w s\n");
  threeMomentum R[8] = { {1,1,1}, {-1,-1,-1},
			 {1,1,-1}, {-1,-1,1},
			 {1,-1,1}, {-1,1,-1},
			 {-1,1,1}, {1,-1,-1} };
  
  typename DataAllMomentumType::ContainerType out(raw_data.getLt(), raw_data.getNsample()); out.zero();

  for(int psnk=0;psnk<8;psnk++)
    for(int psrc=0;psrc<8;psrc++)
      out = out + raw_data(fig, momComb(R[psnk], R[psrc]));

  out = out / double(64);
  return out;
}

template<typename DataAllMomentumType, typename BubbleDataType>
void computeV(DataAllMomentumType &raw_data, const BubbleDataType &raw_bubble_data, const int tsep_pipi){
  (std::cout << "Computing V diagrams with BubbleDataType = " << printType<BubbleDataType>() << " and " << omp_get_max_threads() << " threads\n").flush(); 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed V diagrams with BubbleType = ") + printType<BubbleDataType>() + " in %w s\n");
				 
  threeMomentum R[8] = { {1,1,1}, {-1,-1,-1},
			 {1,1,-1}, {-1,-1,1},
			 {1,-1,1}, {-1,1,-1},
			 {-1,1,1}, {1,-1,-1} };
  const int Lt = raw_bubble_data.getLt();
  const int Nsample = raw_bubble_data.getNsample();
  raw_data.setup(Lt,Nsample);
  
  //Populate output
  for(int psnk=0;psnk<8;psnk++)
    for(int psrc=0;psrc<8;psrc++)
      auto &into = raw_data('V',momComb(R[psnk],R[psrc]));

#pragma omp parallel for
  for(int pp=0;pp<8*8;pp++){
    int psnk = pp / 8;
    int psrc = pp % 8;

    //for(int psnk=0;psnk<8;psnk++)
    //for(int psrc=0;psrc<8;psrc++){  
      //B(tsrc + tsep + tsep_pipi, -p1_snk) B(tsrc, p1_src)
      const auto &Bmp1_snk = raw_bubble_data( -R[psnk] );
      const auto &Bp1_src  = raw_bubble_data(  R[psrc] );

      auto &into = raw_data('V',momComb(R[psnk],R[psrc]));
      for(int tsrc=0;tsrc<Lt;tsrc++)
	for(int tsep=0;tsep<Lt;tsep++)
	  into(tsrc,tsep) = Bmp1_snk( (tsrc + tsep + tsep_pipi) % Lt ) * Bp1_src( tsrc );
  }
}



  
template<typename FigureDataType>
auto sourceAverage(const FigureDataType & data)->correlationFunction<double,typename std::decay<decltype(data(0,0))>::type>{
  typedef typename std::decay<decltype(data(0,0))>::type DistributionType;

  int Lt = data.getLt();
  int nsample = data.getNsample();
  correlationFunction<double,DistributionType> into(data.getLt(),
					     [nsample](int i) {  return typename correlationFunction<double,DistributionType>::ElementType(i, DistributionType(nsample,0.)); }
					     );
  std::vector<int> tsrc_include;
  for(int tsrc=0;tsrc<Lt;tsrc++){
    bool is_nonzero = !data.isZero(tsrc);
    if(is_nonzero)
      tsrc_include.push_back(tsrc);
  }
  double N(tsrc_include.size());

  for(int tsep=0;tsep<Lt;tsep++){
    auto & v = into.value(tsep);
    v.zero();
    for(int i=0;i<tsrc_include.size();i++)
      v = v + data(tsrc_include[i],tsep);
    v = v/N;
  }
  return into;
}

//Combine the computation of the V diagram with A2 projection and source average to avoid large intermediate data storage
template<typename BubbleDataType>
auto computeVprojectA2sourceAvg(const BubbleDataType &raw_bubble_data, const int tsep_pipi)->correlationFunction<double,typename std::decay<decltype(raw_bubble_data(*((threeMomentum*)NULL))(0))>::type>{
  (std::cout << "Computing A2-projected, src-averaged V diagrams with BubbleDataType = " << printType<BubbleDataType>() << " and " << omp_get_max_threads() << " threads\n").flush(); 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed A2-projected, src-averaged V diagrams with BubbleType = ") + printType<BubbleDataType>() + " in %w s\n");

  threeMomentum R[8] = { {1,1,1}, {-1,-1,-1},
			 {1,1,-1}, {-1,-1,1},
			 {1,-1,1}, {-1,1,-1},
			 {-1,1,1}, {1,-1,-1} };
  const int Lt = raw_bubble_data.getLt();
  const int Nsample = raw_bubble_data.getNsample();

  typedef typename std::decay<decltype(raw_bubble_data(*((threeMomentum*)NULL))(0))>::type  DistributionType;
  correlationFunction<double,DistributionType> out(Lt,
					    [&](const int t)
					    {
					      return typename correlationFunction<double,DistributionType>::ElementType(double(t), DistributionType(Nsample,0.));
					    }
					    );
  int nthr = omp_get_max_threads();
  std::vector<correlationFunction<double,DistributionType> > thr_sum(nthr, out);
  
#pragma omp parallel for
  for(int pp=0;pp<8*8;pp++){
    int me = omp_get_thread_num();
    int psnk = pp / 8;
    int psrc = pp % 8;

    const auto &Bmp1_snk = raw_bubble_data( -R[psnk] );
    const auto &Bp1_src  = raw_bubble_data(  R[psrc] );

    for(int tsrc=0;tsrc<Lt;tsrc++)
      for(int tsep=0;tsep<Lt;tsep++)
	thr_sum[me].value(tsep) = thr_sum[me].value(tsep) + Bmp1_snk( (tsrc + tsep + tsep_pipi) % Lt ) * Bp1_src( tsrc );
  }
  for(int tsep=0;tsep<Lt;tsep++){
    for(int thr=0;thr<nthr;thr++)
      out.value(tsep) = out.value(tsep) + thr_sum[thr].value(tsep);
    out.value(tsep) = out.value(tsep)/double(8*8*Lt);
  }
  return out;  
}
  



bubbleDataDoubleJackAllMomenta doubleJackknifeResampleBubble(const bubbleDataAllMomenta &bubbles){
  int Lt = bubbles.getLt();
  bubbleDataDoubleJackAllMomenta out(Lt,bubbles.getNsample());

  for(auto it = bubbles.begin(); it != bubbles.end(); it++){
    const threeMomentum & mom = it->first;
    const typename bubbleDataAllMomenta::ContainerType & raw = it->second;
    for(int t=0;t<Lt;t++)
      out(mom)(t).resample(raw(t));
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

void outputRawData(const std::string &filename, const correlationFunction<double,rawDataDistributionD> &data, const double coeff){
  std::ofstream of(filename.c_str());
  of << std::setprecision(11) << std::scientific;
  int Lt = data.size();
  int nsample = data.value(0).size();

  for(int t=0;t<data.size();t++)
    for(int s=0;s<nsample;s++)
      of << t << " " << s << " " << coeff * data.value(t).sample(s) << " " << 0. << std::endl;
  of.close();
}

doubleJackCorrelationFunction generateData(const Args &args, const CMDline &cmdline){
  const int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc;
  figureDataAllMomenta raw_data;
  bubbleDataAllMomenta raw_bubble_data;

  if(cmdline.load_data_checkpoint){
    loadCheckpoint<boost::archive::binary_iarchive>(raw_data, raw_bubble_data, cmdline.load_data_checkpoint_file);
  }else if(cmdline.load_hdf5_data_checkpoint){
    loadHDF5checkpoint(raw_data, raw_bubble_data, cmdline.load_hdf5_data_checkpoint_file);
  }else if(cmdline.load_text_data_checkpoint){
    loadCheckpoint<boost::archive::text_iarchive>(raw_data, raw_bubble_data, cmdline.load_text_data_checkpoint_file);    
  }else{
    readFigure(raw_data, 'C', args.data_dir, args.tsep_pipi, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, cmdline);
    readFigure(raw_data, 'D', args.data_dir, args.tsep_pipi, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, cmdline);
    readFigure(raw_data, 'R', args.data_dir, args.tsep_pipi, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, cmdline);
    readBubble(raw_bubble_data, args.data_dir, args.tsep_pipi, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, cmdline);
  }

  if(cmdline.save_data_checkpoint){
    saveCheckpoint<boost::archive::binary_oarchive>(raw_data, raw_bubble_data, cmdline.save_data_checkpoint_file);
  }
  if(cmdline.save_text_data_checkpoint){
    saveCheckpoint<boost::archive::text_oarchive>(raw_data, raw_bubble_data, cmdline.save_text_data_checkpoint_file);
  }
  if(cmdline.save_hdf5_data_checkpoint){
    saveHDF5checkpoint(raw_data, raw_bubble_data, cmdline.save_hdf5_data_checkpoint_file);
  }
  
  //Some of Daiqian's old data was measured on every source timeslice while the majority was measured every 8. To fix this discrepancy we explicitly zero the abnormal data  
  zeroUnmeasuredSourceTimeslices(raw_data, 'C', args.tstep_pipi);
  zeroUnmeasuredSourceTimeslices(raw_data, 'D', args.tstep_pipi);
  zeroUnmeasuredSourceTimeslices(raw_data, 'R', args.tstep_pipi);

  computeV(raw_data, raw_bubble_data, args.tsep_pipi);

  figureData A2_C = projectA2('C', raw_data);
  figureData A2_D = projectA2('D', raw_data);
  figureData A2_R = projectA2('R', raw_data);
  figureData A2_V = projectA2('V', raw_data);

  typedef correlationFunction<double,rawDataDistributionD> rawCorrelationFunction;
  
  rawCorrelationFunction A2_realavg_C = sourceAverage(A2_C);
  rawCorrelationFunction A2_realavg_D = sourceAverage(A2_D);
  rawCorrelationFunction A2_realavg_R = sourceAverage(A2_R);
  rawCorrelationFunction A2_realavg_V = sourceAverage(A2_V);

  outputRawData("raw_data_Cpart.dat", A2_realavg_C, 1.);
  outputRawData("raw_data_Dpart.dat", A2_realavg_D, 2.);
  outputRawData("raw_data_Rpart.dat", A2_realavg_R, -6.);
  outputRawData("raw_data_Vpart.dat", A2_realavg_V, 3.);
  
  rawCorrelationFunction pipi_raw = 2*A2_realavg_D + A2_realavg_C - 6*A2_realavg_R + 3*A2_realavg_V;

  std::cout << "Raw data:\n" << pipi_raw << std::endl;

  doubleJackCorrelationFunction pipi_dj(args.Lt,
					[&pipi_raw,nsample](const int t)
					{
					  typename doubleJackCorrelationFunction::ElementType out(t, doubleJackknifeDistributionD(nsample));
					  out.second.resample(pipi_raw.value(t));
					  return out;
					}
					);
  
  if(args.do_vacuum_subtraction){
    bubbleDataDoubleJackAllMomenta dj_bubble_data = doubleJackknifeResampleBubble(raw_bubble_data);
    doubleJackCorrelationFunction A2_realavg_V_dj = computeVprojectA2sourceAvg(dj_bubble_data,args.tsep_pipi);     //sourceAverage(A2_V_dj);
    pipi_dj = pipi_dj - 3*A2_realavg_V_dj;
  }
  
  pipi_dj = fold(pipi_dj, args.tsep_pipi);
  return pipi_dj;
}

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
