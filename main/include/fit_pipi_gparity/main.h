#ifndef _PIPI_MAIN_H__
#define _PIPI_MAIN_H__

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

//User provides a set of total source pipi momenta. The data are read, rotational-state projected, resampled then averaged over the provided set of total momenta. The resulting correlation function is returned
doubleJackCorrelationFunction generateData(const Args &args, const CMDline &cmdline){
  if(args.total_mom.size() == 1 && args.total_mom[0] == threeMomentum({0,0,0})){
    PiPiProject *proj_src = getProjector(args.proj_src, args.total_mom[0]);
    PiPiProject *proj_snk = getProjector(args.proj_snk, args.total_mom[0]);
    PiPiMomAllow* allow = getMomPairFilter(args.allowed_mom, args.total_mom[0]);

    readFigureStationaryPolicy ffn(cmdline.use_symmetric_quark_momenta);
    readBubbleStationaryPolicy bfn_src(cmdline.use_symmetric_quark_momenta,Source);
    readBubbleStationaryPolicy bfn_snk(cmdline.use_symmetric_quark_momenta,Sink);
    doubleJackCorrelationFunction out = generateData(*proj_src, *proj_snk, *allow, args.isospin, args, cmdline, ffn, bfn_src, bfn_snk,"");

    delete proj_src; delete proj_snk; delete allow;
    return std::move(out);
  }else{ //average over the different total momenta provided
    if(cmdline.use_symmetric_quark_momenta) error_exit(std::cout << "getData use_symmetric_quark_momenta not supported for p_tot != (0,0,0)\n");

    doubleJackCorrelationFunction out;
    for(int p=0;p<args.total_mom.size();p++){
      assert(args.total_mom[p] != threeMomentum({0,0,0}));

      PiPiProject *proj_src = getProjector(args.proj_src, args.total_mom[p]);
      PiPiProject *proj_snk = getProjector(args.proj_snk, -args.total_mom[p]);
      PiPiMomAllow* allow = getMomPairFilter(args.allowed_mom, args.total_mom[p]);

      readFigureTianleComovingPolicy ffn(args.total_mom[p]);
      readBubbleTianleComovingPolicy bfn_src(args.total_mom[p],Source);
      readBubbleTianleComovingPolicy bfn_snk(-args.total_mom[p],Sink);
      PiPiProjectAllowOnlyExistingPionMom proj_src_f(*proj_src, args.total_mom[p], args.pion_momenta); //make sure the pion momenta are in the computed set
      PiPiProjectAllowOnlyExistingPionMom proj_snk_f(*proj_snk, -args.total_mom[p], args.pion_momenta);

      if(p==0) out = generateData(proj_src_f, proj_snk_f, *allow, args.isospin, args, cmdline, ffn, bfn_src, bfn_snk, std::string("ptot") + momStr(args.total_mom[p]) );
      else out = out + generateData(proj_src_f, proj_snk_f, *allow, args.isospin, args, cmdline, ffn, bfn_src, bfn_snk, std::string("ptot") + momStr(args.total_mom[p]));

      delete proj_src; delete proj_snk; delete allow;
    }
    out = out * (1./args.total_mom.size());
    return std::move(out);
  }
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
