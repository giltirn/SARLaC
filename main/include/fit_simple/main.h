#ifndef _FIT_SIMPLE_MAIN_H_
#define _FIT_SIMPLE_MAIN_H_

template<typename DistributionType>
correlationFunction<double, DistributionType> resampleAndCombine(const std::vector<rawDataCorrelationFunctionD> &channels_raw, const Args &args, const CMDline &cmdline){
  const int nchannel = channels_raw.size();

  std::vector<correlationFunction<double, DistributionType> > channels_r(nchannel);
  for(int i=0;i<nchannel;i++){
    channels_r[i] = correlationFunction<double, DistributionType>(args.Lt, [&](const int t){
	return typename correlationFunction<double, DistributionType>::ElementType(t, DistributionType(channels_raw[i].value(t).bin(args.bin_size)));								  });
  }
  correlationFunction<double, DistributionType> out(args.Lt);

  applyCombination(out,channels_r,args.combination);
  applyTimeDep(out, args.outer_time_dep, args.Lt);

  return out;
}

template<>
correlationFunction<double, blockDoubleJackknifeDistributionD> resampleAndCombine<blockDoubleJackknifeDistributionD>(const std::vector<rawDataCorrelationFunctionD> &channels_raw, const Args &args, const CMDline &cmdline){
  const int nchannel = channels_raw.size();

  std::vector<correlationFunction<double, blockDoubleJackknifeDistributionD> > channels_r(nchannel);
  for(int i=0;i<nchannel;i++){
    channels_r[i] = correlationFunction<double, blockDoubleJackknifeDistributionD>(args.Lt, [&](const int t){
	return typename correlationFunction<double, blockDoubleJackknifeDistributionD>::ElementType(t, blockDoubleJackknifeDistributionD(channels_raw[i].value(t), args.bin_size));
      });
  }
  correlationFunction<double, blockDoubleJackknifeDistributionD> out(args.Lt);

  applyCombination(out,channels_r,args.combination);
  applyTimeDep(out, args.outer_time_dep, args.Lt);

  return out;
}
  
#endif
