#ifndef _KTOPIPI_SAMPLE_AMA_DATA_STRUCTS_H_
#define _KTOPIPI_SAMPLE_AMA_DATA_STRUCTS_H_

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

  const NumericTensor<rawDataDistributionD,1> & rawBubble(const char ens, const SloppyExact se) const{
    if(ens == 'S'){
      assert(se == Sloppy);
      return bubble_sloppy_S_binned;
    }else{
      return se == Sloppy ? bubble_sloppy_C_binned : bubble_exact_C_binned;
    }
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

  const RawKtoPiPiData &getRaw(const char ens, const SloppyExact se) const{
    if(ens == 'S'){
      assert(se == Sloppy);
      return raw_sloppy_S;
    }else{
      return se == Sloppy ? raw_sloppy_C : raw_exact_C;
    }
  }
};


#endif
