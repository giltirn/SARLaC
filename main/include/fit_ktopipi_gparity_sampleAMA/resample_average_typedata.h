#ifndef _KTOPIPI_SAMPLE_AMA_RESAMPLE_AVERAGE_TYPE_DATA_H_
#define _KTOPIPI_SAMPLE_AMA_RESAMPLE_AVERAGE_TYPE_DATA_H_



//Perform the source average, sample-AMA resampling and sample-AMA correction
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
    out(&t) = sampleAMAresampleCorrect<DistributionType>(raw_sloppy_S_avg,raw_sloppy_C_avg,raw_exact_C_avg,
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
    out(&t) = sampleAMAresampleCorrect<DistributionType>(raw_sloppy_S_avg,raw_sloppy_C_avg,raw_exact_C_avg,
						inputs.resampler_S,inputs.resampler_C, printer::str(stringize("mix%d (tK avg)(t=%d)",type,t)) );
  }
  return out;
}


//For a specific ens and status, perform the source average and sample-AMA resampling
template<typename DistributionType>
NumericTensor<DistributionType,1> resampleAverageTypeDataSampleAMA(const int type,
								   const int q,
								   const char ens, const SloppyExact se,
								   const allRawData &raw,
								   const allInputs &inputs){
  typedef printOnlyIfType<DistributionType,jackknifeDistributionD> printer;
  const RawKtoPiPiData &raw_data = raw.getRaw(ens,se);
  const std::vector<int> &typedata_nonzerotK = raw_data.nonzerotK(type);
  NumericTensor<DistributionType,1> out({inputs.args.Lt}); //[t]

  for(int t=0;t<inputs.args.Lt;t++){
    rawDataDistributionD raw_avg;
    average(raw_avg, [&](const int i){ return raw_data.A0_alltK(type)({q,typedata_nonzerotK[i],t}); }, typedata_nonzerotK.size());
    out(&t) = sampleAMAresample<DistributionType>::resample(raw_avg,ens,inputs.nS,inputs.nC);
  }
  return out;
}

template<typename DistributionType>
NumericTensor<DistributionType,1> resampleAverageMixDiagramSampleAMA(const int type,
								     const char ens, const SloppyExact se,
								     const allRawData &raw,
								     const allInputs &inputs){
  typedef printOnlyIfType<DistributionType,jackknifeDistributionD> printer;
  const RawKtoPiPiData &raw_data = raw.getRaw(ens,se);
  const std::vector<int> &typedata_nonzerotK = raw_data.nonzerotK(type);
  NumericTensor<DistributionType,1> out({inputs.args.Lt}); //[t]

  for(int t=0;t<inputs.args.Lt;t++){
    rawDataDistributionD raw_avg;
    average(raw_avg, [&](const int i){ return raw_data.mix_alltK(type)({typedata_nonzerotK[i],t}); }, typedata_nonzerotK.size());
    out(&t) = sampleAMAresample<DistributionType>::resample(raw_avg,ens,inputs.nS,inputs.nC);     
  }
  return out;
}


//Accessor should be a struct with   
//const rawDataDistributionD & operator()(const int tK, const int t, const RawKtoPiPiData &raw) const
//std::string descr(const int t) const
//const std::vector<int> & nonzerotK(const RawKtoPiPiData &raw) const
template<typename DistributionType, typename Accessor>
NumericTensor<DistributionType,1> resampleAverageSampleAMA(const allRawData &raw,
							   const allInputs &inputs,
							   const Accessor &accessor
							   ){
  typedef printOnlyIfType<DistributionType,jackknifeDistributionD> printer;
  const std::vector<int> &typedata_nonzerotK = accessor.nonzerotK(raw.raw_sloppy_S);
  NumericTensor<DistributionType,1> out({inputs.args.Lt}); //[t]
  DistributionType sloppy_C, exact_C;
  for(int t=0;t<inputs.args.Lt;t++){
    rawDataDistributionD raw_sloppy_S_avg, raw_sloppy_C_avg, raw_exact_C_avg;
    average(raw_sloppy_S_avg, [&](const int i){ return accessor(typedata_nonzerotK[i], t, raw.raw_sloppy_S); }, typedata_nonzerotK.size());
    average(raw_sloppy_C_avg, [&](const int i){ return accessor(typedata_nonzerotK[i], t, raw.raw_sloppy_C); }, typedata_nonzerotK.size());
    average(raw_exact_C_avg,  [&](const int i){ return accessor(typedata_nonzerotK[i], t, raw.raw_exact_C); }, typedata_nonzerotK.size());
    out(&t) = sampleAMAresampleCorrect<DistributionType>(raw_sloppy_S_avg,raw_sloppy_C_avg,raw_exact_C_avg,
						inputs.resampler_S,inputs.resampler_C, printer::str(accessor.descr(t)) );
  }
  return out;
}

#endif
