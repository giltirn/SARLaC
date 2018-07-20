#ifndef _KTOPIPI_SAMPLE_AMA_ALPHA_VAC_SUB_H_
#define _KTOPIPI_SAMPLE_AMA_ALPHA_VAC_SUB_H_

//This version computes alpha and the vacuum subtractions for each (sloppy/exact, S/C) combination separately
template<typename resampledDistributionType>
void computeAlphaAndVacuumSubtractionsSampleAMASeparately(NumericTensor<resampledDistributionType,1> &alpha,
						NumericTensor<resampledDistributionType,1> &A0_type4_srcavg_vacsub,
						NumericTensor<resampledDistributionType,1> &mix4_srcavg_vacsub,
						const allRawData &raw,
						const allBubbleData &bubble_data,
						const int q,
						const int tsep_k_pi,
						const char ens, const SloppyExact se,
						const allInputs &inputs){
  typedef printOnlyIfType<resampledDistributionType,jackknifeDistributionD> printer;

  //Compute mix4 double-jackknife and tK average
  resampledDistributionType zro(inputs.nS + inputs.nC); zeroit(zro);
  
  const RawKtoPiPiData &raw_data = raw.getRaw(ens,se);
  const NumericTensor<rawDataDistributionD,1> &raw_bubble = bubble_data.rawBubble(ens,se);

  const std::vector<int> &type4_nonzerotK = raw_data.nonzerotK(4);
  const int Lt = inputs.args.Lt;

  alpha = NumericTensor<resampledDistributionType,1>({Lt});
  A0_type4_srcavg_vacsub = NumericTensor<resampledDistributionType,1>({Lt});
  mix4_srcavg_vacsub = NumericTensor<resampledDistributionType,1>({Lt});

  NumericTensor<resampledDistributionType,1> bubble_rs({Lt}, [&](const int *t){ return sampleAMAresample<resampledDistributionType>::resample(raw_bubble(t),ens,inputs.nS,inputs.nC); });

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
      const int tB = (tK + tsep_k_pi) % Lt;
      
      resampledDistributionType mix4_nobub_rs = sampleAMAresample<resampledDistributionType>::resample(raw_data.mix4_alltK_nobub({tK,t}),ens,inputs.nS,inputs.nC);

      resampledDistributionType mix4_vacsub_rs = mix4_nobub_rs * bubble_rs(&tB);
      
      resampledDistributionType A0_type4_nobub_rs =sampleAMAresample<resampledDistributionType>::resample(raw_data.A0_type4_alltK_nobub({q,tK,t}),ens,inputs.nS,inputs.nC);
      
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

//This version computes alpha and the vacuum subtractions by performing the AMA correction first
template<typename resampledDistributionType>
void computeAlphaAndVacuumSubtractionsSampleAMACorrected(NumericTensor<resampledDistributionType,1> &alpha,
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
      const int tB = (tK + tsep_k_pi) % Lt;
      
      resampledDistributionType mix4_nobub_rs = 
	sampleAMAresampleCorrect<resampledDistributionType>(raw.raw_sloppy_S.mix4_alltK_nobub({tK,t}), raw.raw_sloppy_C.mix4_alltK_nobub({tK,t}), raw.raw_exact_C.mix4_alltK_nobub({tK,t}),
						   inputs.resampler_S, inputs.resampler_C, printer::str(stringize("mix4(nobub)(tK=%d,t=%d)",tK,t)) );
      
      resampledDistributionType mix4_vacsub_rs = mix4_nobub_rs * bubble_rs(&tB);
      
      resampledDistributionType A0_type4_nobub_rs =
	sampleAMAresampleCorrect<resampledDistributionType>(raw.raw_sloppy_S.A0_type4_alltK_nobub({q,tK,t}), raw.raw_sloppy_C.A0_type4_alltK_nobub({q,tK,t}), raw.raw_exact_C.A0_type4_alltK_nobub({q,tK,t}),
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


#endif
