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

rawDataDistributionVector computeAMAcorrection(const rawDataDistributionMatrix &sloppy, const rawDataDistributionMatrix &exact){
  const int Lt = exact.size();
  const int nsample = sloppy(0,0).size();
  rawDataDistributionVector out(Lt,rawDataDistributionD(nsample,0.));

  for(int conf=0;conf<nsample;conf++){
    std::vector<int> nonzero_slices = nonZeroSourceTimeSlices(exact, conf);
    std::cout << "On conf " << conf << " exact is nonzero on tsrc={";
    for(int i=0;i<nonzero_slices.size();i++){
      const int tsrc = nonzero_slices[i];
      for(int tsep=0;tsep<Lt;tsep++)
	out[tsep].sample(conf) = out[tsep].sample(conf) + exact(tsrc,tsep).sample(conf) - sloppy(tsrc,tsep).sample(conf);      
      std::cout << tsrc << " ";
    }
    for(int tsep=0;tsep<Lt;tsep++)
      out[tsep].sample(conf) = out[tsep].sample(conf) / double(nonzero_slices.size());    
    std::cout << "}\n";
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
