#pragma once

std::vector<int> nonZeroSourceTimeSlices(const rawDataDistributionMatrix &M, const int conf){
  const int Lt = M.size();
  std::vector<int> nonzero_slices;
  for(int tsrc=0;tsrc<Lt;tsrc++){
    bool allzero = true;
    for(int tsep=0;tsep<Lt;tsep++)
      if(M(tsrc,tsep).sample(conf) != 0. ){ allzero = false; break; }
    if(!allzero) nonzero_slices.push_back(tsrc);
  }
  return nonzero_slices;
}
//assign multiplicities randomly using an RNG seeded from the trajectory index
//int conf : the config index
//traj : the trajectory index (not the config index but the actual trajectory index, used for initializing a unique rng)
std::vector<int> randomMultiplicity(const rawDataDistributionMatrix &M, const int traj, const int conf, const int nexact_timeslices){
  std::vector<int> nonzero_slices = nonZeroSourceTimeSlices(M,conf);
  int Lt = M.size();
  int n_nonzero_slices = nonzero_slices.size();
  int rem = nexact_timeslices - n_nonzero_slices;

  std::vector<int> multiplicity(Lt,0);
  
  //assign mult=1 to all nonzero slices
  for(int tnonzero: nonzero_slices)
    multiplicity[tnonzero] = 1;

  //randomly assign remainder
  RNGstore rng;
  rng.initialize(traj); //seed based on traj to ensure reproducibility across all data on the configuration
  for(int i=0;i<rem;i++){
    int idx = int( floor(uniformRandom<double>(0, n_nonzero_slices, rng)) ); //uniform random is exclusive of the upper bound on the range
    ++multiplicity[ nonzero_slices[idx] ];
  }
  return multiplicity;
}
//multiplicity[conf][t]   for conf \in 0..nsample-1
rawDataDistributionVector computeAMAcorrection(const rawDataDistributionMatrix &sloppy, const rawDataDistributionMatrix &exact, const std::vector< std::vector<int> > &multiplicity, int nexact_timeslices){
  const int Lt = exact.size();
  const int nsample = sloppy(0,0).size();
  rawDataDistributionVector out(Lt,rawDataDistributionD(nsample,0.));

  for(int conf=0;conf<nsample;conf++){
    std::cout << "Conf " << conf << " src_timeslice:multiplicity ";
    int sum_multiplicity = 0; //record the sum of multiplicities to ensure they add up to nexact_timeslices

    std::vector<double> avg_corr(Lt,0.);

    for(int tsrc=0;tsrc < Lt; tsrc++){
      if(multiplicity[conf][tsrc] > 0){
	std::cout << tsrc << ":" << multiplicity[conf][tsrc] << " ";
	sum_multiplicity += multiplicity[conf][tsrc];

	for(int tsep=0;tsep<Lt;tsep++)
	  avg_corr[tsep] += multiplicity[conf][tsrc] * (  exact(tsrc,tsep).sample(conf) - sloppy(tsrc,tsep).sample(conf)  ) / nexact_timeslices;	
      }
    }
    assert( sum_multiplicity == nexact_timeslices );
    
    for(int tsep=0;tsep<Lt;tsep++)
      out[tsep].sample(conf) = out[tsep].sample(conf) + avg_corr[tsep];
    
    std::cout << " for total multiplicity " << sum_multiplicity << std::endl;
  }
  return out;
}

rawDataDistributionMatrix timeReflect(const rawDataDistributionMatrix &m){
  //Boundary is *at* 0. 0->Lt=0, 1->Lt-1, 2->Lt-2 ... Lt-1 -> 1    
  const int Lt = m.size();
  rawDataDistributionMatrix out(Lt);
  for(int tsrc=0;tsrc<Lt;tsrc++)
    for(int tsep=0;tsep<Lt;tsep++){
      //int trefl = tsep == 0 ? 0 : Lt-tsep; 
      int trefl = (Lt - tsep) % Lt;
      out(tsrc,trefl) = m(tsrc,tsep);
    }
  return out;
}
rawDataDistributionVector sourceTimeSliceAverage(const rawDataDistributionMatrix &m){
  const int Lt = m.size();
  const int nsample = m(0,0).size();
  rawDataDistributionVector out(Lt,rawDataDistributionD(nsample));
  for(int tsrc=0;tsrc<Lt;tsrc++)
    for(int tsep=0;tsep<Lt;tsep++)
      out(tsep) = out(tsep) + m(tsrc,tsep);
    
  return out/double(Lt);
}

template<typename jackknifeTimeSeriesType>
jackknifeTimeSeriesType resampleVector(const rawDataDistributionVector &data, const int Lt){
  jackknifeTimeSeriesType out;
  if(data.size() == 0) return out;
  else if(data.size() != Lt) error_exit(std::cout << "resample called on data vector of size " << data.size() << ". Expected 0 or Lt=" << Lt << std::endl);

  out.resize(Lt);
  for(int t=0;t<Lt;t++){
    out.value(t).resample(data[t]);
    out.coord(t) = t;
  }
  return out;
}
