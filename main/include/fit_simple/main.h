#ifndef _FIT_SIMPLE_MAIN_H_
#define _FIT_SIMPLE_MAIN_H_

template<typename DistributionType>
correlationFunction<double, DistributionType> resampleAndCombine(const std::vector<rawDataCorrelationFunctionD> &channels_raw,
								 const int Lt, const int bin_size, Combination combination, TimeDependence outer_time_dep){
  const int nchannel = channels_raw.size();

  std::vector<correlationFunction<double, DistributionType> > channels_r(nchannel);
  for(int i=0;i<nchannel;i++){
    channels_r[i] = correlationFunction<double, DistributionType>(Lt, [&](const int t){
	return typename correlationFunction<double, DistributionType>::ElementType(t, DistributionType(channels_raw[i].value(t).bin(bin_size)));								  });
  }
  correlationFunction<double, DistributionType> out(Lt);

  applyCombination(out,channels_r,combination);
  applyTimeDep(out, outer_time_dep, Lt);

  return out;
}

template<>
correlationFunction<double, blockDoubleJackknifeDistributionD> 
resampleAndCombine<blockDoubleJackknifeDistributionD>(const std::vector<rawDataCorrelationFunctionD> &channels_raw, 
						      const int Lt, const int bin_size, Combination combination, TimeDependence outer_time_dep){ 
  const int nchannel = channels_raw.size();

  std::vector<correlationFunction<double, blockDoubleJackknifeDistributionD> > channels_r(nchannel);
  for(int i=0;i<nchannel;i++){
    channels_r[i] = correlationFunction<double, blockDoubleJackknifeDistributionD>(Lt, [&](const int t){
	return typename correlationFunction<double, blockDoubleJackknifeDistributionD>::ElementType(t, blockDoubleJackknifeDistributionD(channels_raw[i].value(t), bin_size));
      });
  }
  correlationFunction<double, blockDoubleJackknifeDistributionD> out(Lt);

  applyCombination(out,channels_r,combination);
  applyTimeDep(out, outer_time_dep, Lt);

  return out;
}

void transformRaw(std::vector<rawDataCorrelationFunctionD> &channels_raw, const Args &args, const CMDline &cmdline){
  int nchannel = channels_raw.size();

  if(cmdline.remove_samples_in_range){
    std::cout << "Removing samples in range [" << cmdline.remove_samples_in_range_start << ", " <<  cmdline.remove_samples_in_range_lessthan << ")" << std::endl;
    for(int c=0;c<nchannel;c++)
      for(int t=0;t<channels_raw[c].size();t++)
	channels_raw[c].value(t) = removeSamplesInRange(channels_raw[c].value(t), cmdline.remove_samples_in_range_start, cmdline.remove_samples_in_range_lessthan);
  }
  if(cmdline.scramble_raw_data){ //useful as a check to see if binning is actually doing anything more than reducing resolution on the covariance matrix
    int nsample = channels_raw[0].value(0).size();
    if(!RNG.isInitialized()) RNG.initialize(1234);
    std::vector<int> reord(nsample);
    std::list<int> rem; 
    for(int i=0;i<nsample;i++) rem.push_back(i);
    
    for(int i=0;i<nsample;i++){
      int off = (int)uniformRandom<float>(0,rem.size());
      auto it = std::next(rem.begin(), off);
      reord[i] = *it;
      rem.erase(it);
    }
    
    //Check indices are unique
    std::cout << "Reordered samples: ";
    std::set<int> con;
    for(int i=0;i<nsample;i++){
      std::cout << reord[i] << " ";
      con.insert(reord[i]);
      }
    std::cout << std::endl;
    assert(con.size() == nsample); 
    
    for(int c=0;c<nchannel;c++)
      for(int t=0;t<channels_raw[c].size();t++)
	channels_raw[c].value(t) = rawDataDistributionD(nsample, [&](const int s){ return channels_raw[c].value(t).sample(reord[s]); });
  }
}

  
#endif
