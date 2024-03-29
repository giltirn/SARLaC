#include <pipi_common/pipi_common.h>

using namespace SARLaC;

#include <fit_pipi_gparity_sampleAMA/args.h>
#include <fit_pipi_gparity_sampleAMA/cmdline.h>

struct rawData{ //raw, unbinned data
  bubbleDataAllMomenta raw_bubble_data_sloppy_S;
  rawDataCorrelationFunctionD pipi_raw_sloppy_S;

  bubbleDataAllMomenta raw_bubble_data_sloppy_C;
  rawDataCorrelationFunctionD pipi_raw_sloppy_C;

  bubbleDataAllMomenta raw_bubble_data_exact_C;
  rawDataCorrelationFunctionD pipi_raw_exact_C;

  int nS; //unbinned
  int nC;

  rawData(const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk, const int isospin, const ArgsSampleAMA &args, const CMDlineSampleAMA &cmdline){
    bubbleDataAllMomenta* raw_bubble_data[3] = { &raw_bubble_data_sloppy_S, &raw_bubble_data_sloppy_C, &raw_bubble_data_exact_C };
    rawDataCorrelationFunctionD* pipi_raw[3] = { &pipi_raw_sloppy_S, &pipi_raw_sloppy_C, &pipi_raw_exact_C };

    const char ens[3] = { 'S', 'C', 'C' };
    const SloppyExact se[3] = { Sloppy, Sloppy, Exact };
    const std::string descr[3] = { "Sloppy_S", "Sloppy_C", "Exact_C" };
    
    for(int i=0;i<3;i++){
      figureDataAllMomenta raw_data;
      readFigureStationaryPolicy ffn(se[i] == Exact);
      readBubbleStationaryPolicy bfn_src(se[i] == Exact,Source);
      readBubbleStationaryPolicy bfn_snk(se[i] == Exact,Sink);

      int traj_start, traj_inc, traj_lessthan;
      args.traj_info(traj_start, traj_inc, traj_lessthan, ens[i]);
      
      std::string checkpoint_filename;
      bool load_checkpoint = cmdline.load_checkpoint(checkpoint_filename,se[i],ens[i]);

      if(load_checkpoint) loadRawDataCheckpoint(raw_data, *raw_bubble_data[i], checkpoint_filename, descr[i]);
      else{
	std::cout << "Reading raw data " << descr[i] << std::endl;
	
	readRawPiPi2ptData(raw_data, *raw_bubble_data[i], 
			   ffn, bfn_src, bfn_snk, 
			   args.data_dir(ens[i]), traj_start, traj_inc, traj_lessthan, 
			   args.Lt, args.tstep_pipi, 
			   args.tsep_pipi, proj_src, proj_snk);
      }	

      bool save_checkpoint = cmdline.save_checkpoint(checkpoint_filename,se[i],ens[i]);

      if(save_checkpoint) saveRawDataCheckpoint(raw_data, *raw_bubble_data[i], checkpoint_filename, descr[i]);

      combineRawPiPiContractionsOpts opts;
      opts.output_raw_data_bin_size = args.bin_size;
      combineRawPiPiContractions(*pipi_raw[i], raw_data, proj_src, proj_snk, isospin, opts);
    }

    nS = pipi_raw_sloppy_S.value(0).size();
    nC = pipi_raw_sloppy_C.value(0).size();
  }
};
    



template<typename DistributionType>
bubbleDataAllMomentaBase<bubbleDataBase<DistributionType> > resampleCorrectBubbleSampleAMA(const rawData &raw, const int nS, const int nC, const int bin_size){
  const int Lt = raw.raw_bubble_data_sloppy_S.getLt();
  const int tsep_pipi = raw.raw_bubble_data_sloppy_S.getTsepPiPi();

  const int nSb = nS/bin_size;
  const int nCb = nC/bin_size;

  bubbleDataAllMomentaBase<bubbleDataBase<DistributionType> > out(Lt,tsep_pipi);

  for(auto it = raw.raw_bubble_data_sloppy_S.begin(); it != raw.raw_bubble_data_sloppy_S.end(); ++it){
    for(int t=0;t<Lt;t++){      
      DistributionType sloppy_S = sampleAMAresample<DistributionType>::resample(raw.raw_bubble_data_sloppy_S(it->first)(t).bin(bin_size),'S',nSb,nCb);
      DistributionType sloppy_C = sampleAMAresample<DistributionType>::resample(raw.raw_bubble_data_sloppy_C(it->first)(t).bin(bin_size),'C',nSb,nCb);
      DistributionType exact_C = sampleAMAresample<DistributionType>::resample(raw.raw_bubble_data_exact_C(it->first)(t).bin(bin_size),'C',nSb,nCb);
      out(it->first)(t) = sloppy_S + exact_C - sloppy_C;
    }
  }
  return out;
}

template<typename DistributionType>
void resampleCombineData(correlationFunction<double,DistributionType> &pipi, 
			 const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk, const int isospin,
			 const rawData &raw, const ArgsSampleAMA &args, const CMDlineSampleAMA &cmdline){
  typedef correlationFunction<double,DistributionType> CorrFunc;
  typedef typename CorrFunc::ElementType Elem;
  pipi.resize(args.Lt);
  const int nS = raw.nS;
  const int nC = raw.nC;

  const int nSb = nS/args.bin_size;
  const int nCb = nC/args.bin_size;

  DistributionType pipi_sloppy_S_r, pipi_sloppy_C_r, pipi_exact_C_r;
  for(int t=0;t<args.Lt;t++){
    pipi.coord(t) = t;

    pipi_sloppy_S_r = sampleAMAresample<DistributionType>::resample(raw.pipi_raw_sloppy_S.value(t).bin(args.bin_size),'S',nSb,nCb);
    pipi_sloppy_C_r = sampleAMAresample<DistributionType>::resample(raw.pipi_raw_sloppy_C.value(t).bin(args.bin_size),'C',nSb,nCb);
    pipi_exact_C_r = sampleAMAresample<DistributionType>::resample(raw.pipi_raw_exact_C.value(t).bin(args.bin_size),'C',nSb,nCb);

    pipi.value(t) = pipi_sloppy_S_r + pipi_exact_C_r - pipi_sloppy_C_r;
  }
  if(isospin == 0 && args.do_vacuum_subtraction){
    auto bubble_data_corrected_r = resampleCorrectBubbleSampleAMA<DistributionType>(raw,nS,nC,args.bin_size);
    CorrFunc A2_realavg_V_r = computePiPi2ptFigureVprojectSourceAvg(bubble_data_corrected_r,args.tsep_pipi, proj_src, proj_snk);
    pipi = pipi - 3*A2_realavg_V_r;
  }
}

void generateData(jackknifeCorrelationFunctionD &pipi_j, doubleJackknifeCorrelationFunctionD &pipi_dj, 
		  const PiPiProjectorBase &proj_src,  const PiPiProjectorBase &proj_snk, const int isospin,
		  const ArgsSampleAMA &args, const CMDlineSampleAMA &cmdline){
  rawData raw(proj_src, proj_snk, isospin, args,cmdline);
  resampleCombineData<jackknifeDistributionD>(pipi_j,proj_src,proj_snk, isospin, raw,args,cmdline);
  resampleCombineData<doubleJackknifeDistributionD>(pipi_dj,proj_src,proj_snk, isospin, raw,args,cmdline);
  
  pipi_j = fold(pipi_j, 2*args.tsep_pipi);
  pipi_dj = fold(pipi_dj, 2*args.tsep_pipi);
}

void getData(const PiPiProjectorBase &proj_src,  const PiPiProjectorBase &proj_snk, const int isospin, 
	     jackknifeCorrelationFunctionD &pipi_j, doubleJackknifeCorrelationFunctionD &pipi_dj, const ArgsSampleAMA &args, const CMDlineSampleAMA &cmdline){
  if(cmdline.load_combined_data){
#ifdef HAVE_HDF5
    HDF5reader reader(cmdline.load_combined_data_file);
    read(reader, pipi_j, "pipi_j");
    read(reader, pipi_dj, "pipi_dj");
#else
    error_exit("getData: Reading amplitude data requires HDF5\n");
#endif
  }else{
    generateData(pipi_j,pipi_dj, proj_src,proj_snk, isospin,args,cmdline);
  }

  if(cmdline.save_combined_data){
#ifdef HAVE_HDF5
    HDF5writer writer(cmdline.save_combined_data_file);
    write(writer, pipi_j, "pipi_j");
    write(writer, pipi_dj, "pipi_dj");
#else
    error_exit("getData: Saving amplitude data requires HDF5\n");
#endif
  }
}

int main(const int argc, const char* argv[]){
  ArgsSampleAMA args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  {
    const std::string arg_file = argv[1];
    parse(args,arg_file);
    std::cout << "Read arguments: \n" << args << std::endl;
  }
  
  CMDlineSampleAMA cmdline(argc,argv,2);

  const threeMomentum ptot = {0,0,0};
  std::unique_ptr<PiPiProjectorBase> proj_src(getProjector(args.proj_src,ptot));
  std::unique_ptr<PiPiProjectorBase> proj_snk(getProjector(args.proj_snk,ptot));

  //Read resampled
  jackknifeCorrelationFunctionD pipi_j;
  doubleJackknifeCorrelationFunctionD pipi_dj;

  getData(*proj_src, *proj_snk, args.isospin, pipi_j,pipi_dj,args,cmdline);
   
  //Filter out the data that is to be fitted
  const int nsample = pipi_j.value(0).size();
  
  filterXrange<double> trange(args.t_min,args.t_max);
  
  doubleJackknifeCorrelationFunctionD pipi_dj_inrange;
  jackknifeCorrelationFunctionD pipi_j_inrange;
  for(int d=0;d<pipi_dj.size();d++)
    if(trange.accept(pipi_dj.coord(d),pipi_dj.value(d) )){
      pipi_dj_inrange.push_back(pipi_dj[d]);
      pipi_j_inrange.push_back(pipi_j[d]);
    }
  
  //Perform the fit
  pipiFitOptions opt; cmdline.Export(opt);
  std::pair<jackknifeDistributionD,jackknifeDistributionD> Epipi_and_const = fit(pipi_j_inrange,pipi_dj_inrange,
										 args.fitfunc, args.correlated, args.Lt, args.tsep_pipi, args.Ascale, args.Cscale, opt);

  plot(pipi_j,Epipi_and_const.first,Epipi_and_const.second,
       args.t_min, args.t_max, args.effective_energy, args.Lt, args.tsep_pipi, args.Ascale, args.Cscale);

  std::cout << "Done\n";
  return 0;
}

