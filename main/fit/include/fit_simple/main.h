#ifndef _FIT_SIMPLE_MAIN_H_
#define _FIT_SIMPLE_MAIN_H_

template<typename DistributionType, typename BinResampler>
correlationFunction<double, DistributionType> resampleAndCombine(const std::vector<rawDataCorrelationFunctionD> &channels_raw,
								 const int Lt, Combination combination, TimeDependence outer_time_dep, const BinResampler &bin_resampler){
  const int nchannel = channels_raw.size();

  std::vector<correlationFunction<double, DistributionType> > channels_r(nchannel);
  for(int i=0;i<nchannel;i++){
    channels_r[i].resize(Lt);
    for(int t=0;t<Lt;t++){
      channels_r[i].coord(t) = t;
      bin_resampler.binResample(channels_r[i].value(t), channels_raw[i].value(t));
    }
  }

  correlationFunction<double, DistributionType> out(Lt);

  applyCombination(out,channels_r,combination);
  applyTimeDep(out, outer_time_dep, Lt);

  return out;
}

template<typename DistributionType>
inline correlationFunction<double, DistributionType> resampleAndCombine(const std::vector<rawDataCorrelationFunctionD> &channels_raw,
									const int Lt, const int bin_size, Combination combination, TimeDependence outer_time_dep){
  return ::resampleAndCombine<DistributionType>(channels_raw, Lt, combination, outer_time_dep, basicBinResampler(bin_size));
}


struct transformOptions{
  bool remove_samples_in_range;
  int remove_samples_in_range_start;
  int remove_samples_in_range_lessthan;

  bool scramble_raw_data;
};


void transformRaw(std::vector<rawDataCorrelationFunctionD> &channels_raw, const transformOptions &opt){
  int nchannel = channels_raw.size();

  if(opt.remove_samples_in_range){
    std::cout << "Removing samples in range [" << opt.remove_samples_in_range_start << ", " <<  opt.remove_samples_in_range_lessthan << ")" << std::endl;
    for(int c=0;c<nchannel;c++)
      for(int t=0;t<channels_raw[c].size();t++)
	channels_raw[c].value(t) = removeSamplesInRange(channels_raw[c].value(t), opt.remove_samples_in_range_start, opt.remove_samples_in_range_lessthan);
  }
  if(opt.scramble_raw_data){ //useful as a check to see if binning is actually doing anything more than reducing resolution on the covariance matrix
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


template<typename DataSeriesType>
inline DataSeriesType getDataInRange(const DataSeriesType &data, const int tmin, const int tmax){
  return DataSeriesType(tmax - tmin + 1, [&](const int i){ return data[tmin + i]; });
}

template<typename Out, typename In>
Out pconvert(const In &in){ Out out(in.size()); for(int i=0;i<in.size();i++) out(i) = in(i); return out; }

template<typename Out, typename In>
jackknifeDistribution<Out> pconvert(const jackknifeDistribution<In> &in){ 
  jackknifeDistribution<Out> out(in.size());
  for(int s=0;s<in.size();s++) out.sample(s) = pconvert<Out,In>(in.sample(s));
  return out;
}

template<typename Out, typename In>
bootstrapDistribution<Out> pconvert(const bootstrapDistribution<In> &in){ 
  bootstrapDistribution<Out> out(in.getInitializer());
  for(int s=0;s<in.size();s++) out.sample(s) = pconvert<Out,In>(in.sample(s));
  out.best() = pconvert<Out,In>(in.best());
  return out;
}
  


#endif
