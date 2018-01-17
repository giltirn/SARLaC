#ifndef FIT_KTOPIPI_GPARITY_SAMPLEAMA_MAIN_H__
#define FIT_KTOPIPI_GPARITY_SAMPLEAMA_MAIN_H__

class sampleAMA_resampler{
  char ens;
  int nS;
  int nC;
public:
  sampleAMA_resampler(){}
  sampleAMA_resampler(char _ens, int _nS, int _nC): ens(_ens), nS(_nS), nC(_nC){}

  template<typename DistributionType>
  inline void resample(DistributionType &out, const rawDataDistributionD &in) const{ 
    out = sampleAMAresample<DistributionType>::resample(in,ens,nS,nC);
  }
};


struct allInputs{
  SampleAMAargs args;
  SampleAMAcmdLine cmdline;

  Args args_S;
  Args args_C;
  CMDline cmdline_sloppy_S;
  CMDline cmdline_sloppy_C;
  CMDline cmdline_exact_C;

  int nS;
  int nC;

  sampleAMA_resampler resampler_S;
  sampleAMA_resampler resampler_C;

  allInputs(const SampleAMAargs &_args, const SampleAMAcmdLine &_cmdline): args(_args), cmdline(_cmdline){
    nS = (args.traj_lessthan_S - args.traj_start_S)/args.traj_inc/args.bin_size;
    nC = (args.traj_lessthan_C - args.traj_start_C)/args.traj_inc/args.bin_size;

    resampler_S = sampleAMA_resampler('S',nS,nC);
    resampler_C = sampleAMA_resampler('C',nS,nC);

    args_S = args.toArgs('S');
    args_C = args.toArgs('C');
    
    cmdline_sloppy_S = cmdline.toCMDline('S',Sloppy);
    cmdline_sloppy_C = cmdline.toCMDline('C',Sloppy);
    cmdline_exact_C = cmdline.toCMDline('C',Exact);
  }
};

//#define PRINT_CORRECTION

template<typename resampledDistributionType>
resampledDistributionType resampleCorrect(const rawDataDistributionD &sloppy_S, const rawDataDistributionD &sloppy_C, const rawDataDistributionD &exact_C, 
					  const sampleAMA_resampler &resampler_S, const sampleAMA_resampler &resampler_C, const std::string &descr = ""){
  resampledDistributionType out_r, sloppy_S_r, sloppy_C_r, exact_C_r;
  resampler_S.resample(sloppy_S_r, sloppy_S);
  resampler_C.resample(sloppy_C_r, sloppy_C);
  resampler_C.resample(exact_C_r, exact_C);
  
  out_r = sloppy_S_r + exact_C_r - sloppy_C_r;

#ifdef PRINT_CORRECTION
  if(descr != ""){
    resampledDistributionType diff = out_r - sloppy_S_r;
    resampledDistributionType reldiff = (out_r - sloppy_S_r)/sloppy_S_r;
    std::cout << descr << " corrected:" << out_r << " sloppy:" << sloppy_S_r << " diff:" << diff << " reldiff:" << reldiff << std::endl;
  }
#endif
  return out_r;
}

template<typename T, typename U>
struct printOnlyIfType{
  inline static std::string str(const std::string &str){ return ""; }
};
template<typename T>
struct printOnlyIfType<T,T>{
  inline static std::string str(const std::string &str){ return str; }
};


struct allBubbleData{
  NumericTensor<rawDataDistributionD,1> bubble_sloppy_S;  
  NumericTensor<rawDataDistributionD,1> bubble_sloppy_S_binned;

  NumericTensor<rawDataDistributionD,1> bubble_sloppy_C;
  NumericTensor<rawDataDistributionD,1> bubble_sloppy_C_binned;

  NumericTensor<rawDataDistributionD,1> bubble_exact_C;
  NumericTensor<rawDataDistributionD,1> bubble_exact_C_binned;

  //full superjackknife
  NumericTensor<jackknifeDistributionD,1> bubble_j;
  NumericTensor<doubleJackknifeDistributionD,1> bubble_dj;


  allBubbleData(const allInputs &inputs): bubble_j({inputs.args.Lt}), bubble_dj({inputs.args.Lt}){
    bubble_sloppy_S = getA2projectedBubble(inputs.args_S,inputs.cmdline_sloppy_S);
    bubble_sloppy_C = getA2projectedBubble(inputs.args_C,inputs.cmdline_sloppy_C);
    bubble_exact_C = getA2projectedBubble(inputs.args_C,inputs.cmdline_exact_C);

    bubble_sloppy_S_binned = bin(bubble_sloppy_S,inputs.args.bin_size);
    bubble_sloppy_C_binned = bin(bubble_sloppy_C,inputs.args.bin_size);
    bubble_exact_C_binned = bin(bubble_exact_C,inputs.args.bin_size);

    for(int t=0;t<inputs.args.Lt;t++){
      bubble_j(&t) = resampleCorrect<jackknifeDistributionD>(bubble_sloppy_S_binned(&t), bubble_sloppy_C_binned(&t), bubble_exact_C_binned(&t), 
							     inputs.resampler_S, inputs.resampler_C, stringize("Bubble(%d)",t));
      bubble_dj(&t) = resampleCorrect<doubleJackknifeDistributionD>(bubble_sloppy_S_binned(&t), bubble_sloppy_C_binned(&t), bubble_exact_C_binned(&t), 
							     inputs.resampler_S, inputs.resampler_C);
    }
  }
  
  BubbleData extractUnbinnedBubble(const char ens, const SloppyExact se) const{
    BubbleData out;
    if(ens == 'S'){
      assert(se == Sloppy);
      out.bubble = bubble_sloppy_S;
    }else{
      out.bubble = se == Sloppy ? bubble_sloppy_C : bubble_exact_C;
    }
    return out;
  }

};
template<typename DistributionType> struct getResampledBubbleSampleAMA{};
template<> struct getResampledBubbleSampleAMA<jackknifeDistributionD>{ 
  static inline const NumericTensor<jackknifeDistributionD,1> &get(const allBubbleData &bubble_data){ return bubble_data.bubble_j; }  
};
template<> struct getResampledBubbleSampleAMA<doubleJackknifeDistributionD>{
  static inline const NumericTensor<doubleJackknifeDistributionD,1> &get(const allBubbleData &bubble_data){ return bubble_data.bubble_dj; }  
};


struct allRawData{
  RawKtoPiPiData raw_sloppy_S;
  RawKtoPiPiData raw_sloppy_C;
  RawKtoPiPiData raw_exact_C;

  allRawData(const allBubbleData &bubble_data, const int tsep_k_pi, const allInputs &inputs){
    BubbleData bubble_sloppy_S = bubble_data.extractUnbinnedBubble('S',Sloppy);
    BubbleData bubble_sloppy_C = bubble_data.extractUnbinnedBubble('C',Sloppy);
    BubbleData bubble_exact_C = bubble_data.extractUnbinnedBubble('C',Exact);
    
    std::cout << "allRawData loading sloppy_S" << std::endl;
    raw_sloppy_S = RawKtoPiPiData(tsep_k_pi, bubble_sloppy_S, inputs.args_S, inputs.cmdline_sloppy_S);
    std::cout << "allRawData loading sloppy_C" << std::endl;
    raw_sloppy_C = RawKtoPiPiData(tsep_k_pi, bubble_sloppy_C, inputs.args_C, inputs.cmdline_sloppy_C);
    std::cout << "allRawData loading exact_C" << std::endl;
    raw_exact_C = RawKtoPiPiData(tsep_k_pi, bubble_exact_C, inputs.args_C, inputs.cmdline_exact_C);
  }
};




template<typename resampledDistributionType>
void computeAlphaAndVacuumSubtractionsSampleAMA(NumericTensor<resampledDistributionType,1> &alpha,
						NumericTensor<resampledDistributionType,1> &A0_type4_srcavg_vacsub,
						NumericTensor<resampledDistributionType,1> &mix4_srcavg_vacsub,
						const allRawData &raw,
						const NumericTensor<resampledDistributionType,1> &bubble_rs,
						const int q,
						const int tsep_k_pi,
						const allInputs &inputs){
  typedef printOnlyIfType<resampledDistributionType,jackknifeDistributionD> printer;

  //Compute mix4 double-jackknife and tK average
  resampledDistributionType zro = bubble_rs({0}); zeroit(zro);
  
  const std::vector<int> &type4_nonzerotK = raw.raw_sloppy_S.nonzerotK(4);
  const int Lt = inputs.args.Lt;

#ifndef PRINT_CORRECTION
#pragma omp parallel for
#endif
  for(int t=0;t<Lt;t++){

    mix4_srcavg_vacsub(&t) = zro;
    A0_type4_srcavg_vacsub(&t) = zro;

    resampledDistributionType mix4_nobub_srcavg = zro;
    resampledDistributionType A0_type4_nobub_srcavg = zro;
    
    for(int ii=0;ii<type4_nonzerotK.size();ii++){
      const int tK = type4_nonzerotK[ii];
      const int tB = (tK + tsep_k_pi + inputs.args.tsep_pipi) % Lt;
      
      resampledDistributionType mix4_nobub_rs = 
	resampleCorrect<resampledDistributionType>(raw.raw_sloppy_S.mix4_alltK_nobub({tK,t}), raw.raw_sloppy_C.mix4_alltK_nobub({tK,t}), raw.raw_exact_C.mix4_alltK_nobub({tK,t}),
						   inputs.resampler_S, inputs.resampler_C, printer::str(stringize("mix4(nobub)(tK=%d,t=%d)",tK,t)) );
      
      resampledDistributionType mix4_vacsub_rs = mix4_nobub_rs * bubble_rs(&tB);
      
      resampledDistributionType A0_type4_nobub_rs =
	resampleCorrect<resampledDistributionType>(raw.raw_sloppy_S.A0_type4_alltK_nobub({q,tK,t}), raw.raw_sloppy_C.A0_type4_alltK_nobub({q,tK,t}), raw.raw_exact_C.A0_type4_alltK_nobub({q,tK,t}),
						   inputs.resampler_S, inputs.resampler_C, printer::str(stringize("Q=%d type4(nobub)(tK=%d,t=%d)",q+1,tK,t)) );
      
      resampledDistributionType A0_type4_vacsub_rs = A0_type4_nobub_rs * bubble_rs(&tB);

      mix4_srcavg_vacsub(&t) = mix4_srcavg_vacsub(&t) + mix4_vacsub_rs;
      mix4_nobub_srcavg = mix4_nobub_srcavg + mix4_nobub_rs;

      A0_type4_srcavg_vacsub(&t) = A0_type4_srcavg_vacsub(&t) + A0_type4_vacsub_rs;
      A0_type4_nobub_srcavg = A0_type4_nobub_srcavg + A0_type4_nobub_rs;
    }
    double n(type4_nonzerotK.size());
    
    mix4_srcavg_vacsub(&t) = mix4_srcavg_vacsub(&t)/n;
    A0_type4_srcavg_vacsub(&t) = A0_type4_srcavg_vacsub(&t)/n;
    mix4_nobub_srcavg = mix4_nobub_srcavg/n;
    A0_type4_nobub_srcavg = A0_type4_nobub_srcavg/n;

    alpha(&t) = A0_type4_nobub_srcavg/mix4_nobub_srcavg;
  }
}


template<typename DistributionType>
NumericTensor<DistributionType,1> resampleAverageTypeDataSampleAMA(const int type,
								   const int q,
								   const allRawData &raw,
								   const allInputs &inputs){
  typedef printOnlyIfType<DistributionType,jackknifeDistributionD> printer;
  const std::vector<int> &typedata_nonzerotK = raw.raw_sloppy_S.nonzerotK(type);
  NumericTensor<DistributionType,1> out({inputs.args.Lt}); //[t]
  DistributionType sloppy_C, exact_C;
  for(int t=0;t<inputs.args.Lt;t++){
    rawDataDistributionD raw_sloppy_S_avg, raw_sloppy_C_avg, raw_exact_C_avg;
    average(raw_sloppy_S_avg, [&](const int i){ return raw.raw_sloppy_S.A0_alltK(type)({q,typedata_nonzerotK[i],t}); }, typedata_nonzerotK.size());
    average(raw_sloppy_C_avg, [&](const int i){ return raw.raw_sloppy_C.A0_alltK(type)({q,typedata_nonzerotK[i],t}); }, typedata_nonzerotK.size());
    average(raw_exact_C_avg,  [&](const int i){ return raw.raw_exact_C.A0_alltK(type)({q,typedata_nonzerotK[i],t}); }, typedata_nonzerotK.size());
    out(&t) = resampleCorrect<DistributionType>(raw_sloppy_S_avg,raw_sloppy_C_avg,raw_exact_C_avg,
						inputs.resampler_S,inputs.resampler_C, printer::str(stringize("Q=%d type%d (tK avg)(t=%d)",q+1,type,t)) );
  }
  return out;
}

template<typename DistributionType>
NumericTensor<DistributionType,1> resampleAverageMixDiagramSampleAMA(const int type,
								  const allRawData &raw,
								  const allInputs &inputs){
  typedef printOnlyIfType<DistributionType,jackknifeDistributionD> printer;
  const std::vector<int> &typedata_nonzerotK = raw.raw_sloppy_S.nonzerotK(type);
  NumericTensor<DistributionType,1> out({inputs.args.Lt}); //[t]
  DistributionType sloppy_C, exact_C;
  for(int t=0;t<inputs.args.Lt;t++){
    rawDataDistributionD raw_sloppy_S_avg, raw_sloppy_C_avg, raw_exact_C_avg;
    average(raw_sloppy_S_avg, [&](const int i){ return raw.raw_sloppy_S.mix_alltK(type)({typedata_nonzerotK[i],t}); }, typedata_nonzerotK.size());
    average(raw_sloppy_C_avg, [&](const int i){ return raw.raw_sloppy_C.mix_alltK(type)({typedata_nonzerotK[i],t}); }, typedata_nonzerotK.size());
    average(raw_exact_C_avg,  [&](const int i){ return raw.raw_exact_C.mix_alltK(type)({typedata_nonzerotK[i],t}); }, typedata_nonzerotK.size());
    out(&t) = resampleCorrect<DistributionType>(raw_sloppy_S_avg,raw_sloppy_C_avg,raw_exact_C_avg,
						inputs.resampler_S,inputs.resampler_C, printer::str(stringize("mix%d (tK avg)(t=%d)",type,t)) );
  }
  return out;
}




template<typename DistributionType>
NumericTensor<DistributionType,1> computeQamplitudeSampleAMA(const int q, const int tsep_k_pi, const allRawData &raw, const allBubbleData &bubble_data, const allInputs &inputs, const std::string &descr){
  const int Lt = inputs.args.Lt;
  //Compute alpha and type4/mix4 vacuum subtractions
  std::cout << "Computing " << descr << " alpha and vacuum subtractions\n";
  NumericTensor<DistributionType,1> alpha_r({Lt}), A0_type4_srcavg_vacsub_r({Lt}), mix4_srcavg_vacsub_r({Lt}); //[t]
  computeAlphaAndVacuumSubtractionsSampleAMA(alpha_r, A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r,
					     raw, getResampledBubbleSampleAMA<DistributionType>::get(bubble_data),q,tsep_k_pi,inputs);

  //Compute tK-averages type4 and mix4 diagrams from data including bubble-------------//
  std::cout << "Computing " << descr << " tK averages and mix diagrams\n";
  IndexedContainer<NumericTensor<DistributionType,1>, 4, 1> A0_srcavg_r; //[t]
  IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> mix_srcavg_r; //[t]
  for(int i=1;i<=4;i++) A0_srcavg_r(i) = resampleAverageTypeDataSampleAMA<DistributionType>(i,q,raw,inputs); //[t]
  for(int i=3;i<=4;i++) mix_srcavg_r(i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,raw,inputs);

  //Subtract the pseudoscalar operators and mix4 vacuum term
  std::cout << "Subtracting pseudoscalar operators and mix4 vacuum term under " << descr << "\n";
  A0_srcavg_r(3) = A0_srcavg_r(3).transform([&](int const* t, const DistributionType &from){ return DistributionType(from - alpha_r(t)*mix_srcavg_r(3)(t)); }); 
  A0_srcavg_r(4) = A0_srcavg_r(4).transform([&](int const* t, const DistributionType &from){
      return DistributionType(from - alpha_r(t)*( mix_srcavg_r(4)(t) - mix4_srcavg_vacsub_r(t) ) );
    }); 

  //Perform the type 4 vacuum subtraction
  std::cout << "Performing type-4 vacuum subtraction\n";
  A0_srcavg_r(4) = A0_srcavg_r(4) - A0_type4_srcavg_vacsub_r;

  //Get the full double-jackknife amplitude
  std::cout << "Computing full amplitudes\n";
  NumericTensor<DistributionType,1> A0_full_srcavg_r = A0_srcavg_r(1) + A0_srcavg_r(2) + A0_srcavg_r(3) + A0_srcavg_r(4);

  return A0_full_srcavg_r;
}

  
//Read and prepare the data for a particular tsep_k_pi_idx
void getDataSampleAMA(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
		      std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
		      const allBubbleData &bubble_data, const int tsep_k_pi_idx, const allInputs &inputs){
  int tsep_k_pi = inputs.args.tsep_k_pi[tsep_k_pi_idx];
  std::cout << "Getting data for tsep_k_pi = " <<  tsep_k_pi << std::endl;
  printMem("getData called");
  
  allRawData raw(bubble_data, tsep_k_pi, inputs);
  
  for(int q=0;q<10;q++){
    std::cout << "Starting Q" << q+1 << std::endl;

    printMem("Starting new Q");
 
    NumericTensor<doubleJackknifeDistributionD,1> A0_full_srcavg_dj = computeQamplitudeSampleAMA<doubleJackknifeDistributionD>(q, tsep_k_pi, raw, bubble_data, inputs, "double jackknife");
    NumericTensor<jackknifeDistributionD,1> A0_full_srcavg_j = computeQamplitudeSampleAMA<jackknifeDistributionD>(q, tsep_k_pi, raw, bubble_data, inputs, "single jackknife");

    //Insert data into output containers    
    for(int t=0;t<inputs.args.Lt;t++){
      A0_all_j[q].push_back(amplitudeDataCoord(t,tsep_k_pi), A0_full_srcavg_j(&t));
      A0_all_dj[q].push_back(amplitudeDataCoord(t,tsep_k_pi), A0_full_srcavg_dj(&t));
    }
  }
}



void getDataSampleAMA(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
		      std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
		      const allInputs &inputs){
  if(inputs.cmdline.load_amplitude_data){
#ifdef HAVE_HDF5
    HDF5reader reader(inputs.cmdline.load_amplitude_data_file);
    read(reader, A0_all_j, "A0_all_j");
    read(reader, A0_all_dj, "A0_all_dj");
#else
    error_exit("getData: Reading amplitude data requires HDF5\n");
#endif
  }else{
    //Read the bubble data
    allBubbleData bubble_data(inputs);
    scratch scratch_store(inputs.cmdline.use_scratch, inputs.cmdline.use_scratch_stub, inputs.cmdline.use_existing_scratch_files, inputs.args.tsep_k_pi);

    for(int tsep_k_pi_idx=0;tsep_k_pi_idx<inputs.args.tsep_k_pi.size();tsep_k_pi_idx++){
      if(scratch_store.doSkipLoad(tsep_k_pi_idx)) continue;
      getDataSampleAMA(A0_all_j,A0_all_dj,bubble_data,tsep_k_pi_idx,inputs);
      scratch_store.writeScratch(A0_all_j, A0_all_dj, tsep_k_pi_idx);
    }	

    scratch_store.reloadScratch(A0_all_j, A0_all_dj); 
  }

  if(inputs.cmdline.save_amplitude_data){
#ifdef HAVE_HDF5
    HDF5writer writer(inputs.cmdline.save_amplitude_data_file);
    write(writer, A0_all_j, "A0_all_j");
    write(writer, A0_all_dj, "A0_all_dj");
#else
    error_exit("getData: Saving amplitude data requires HDF5\n");
#endif
  }
}

void checkpointRawOnly(const allInputs &inputs){
  allBubbleData bubble_data(inputs);
  for(int tsep_k_pi_idx=0;tsep_k_pi_idx<inputs.args.tsep_k_pi.size();tsep_k_pi_idx++){
    int tsep_k_pi = inputs.args.tsep_k_pi[tsep_k_pi_idx];
    std::cout << "Getting data for tsep_k_pi = " <<  tsep_k_pi << std::endl;
    allRawData raw(bubble_data, tsep_k_pi, inputs);
  }
}
#endif
