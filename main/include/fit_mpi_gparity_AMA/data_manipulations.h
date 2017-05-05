#ifndef _FIT_MPI_GPARITY_AMA_DATA_MANIPULATIONS_H
#define _FIT_MPI_GPARITY_AMA_DATA_MANIPULATIONS_H

std::vector<int> nonZeroSourceTimeSlices(const distributionMatrix &M, const int conf){
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

distributionVector computeAMAcorrection(const distributionMatrix &sloppy, const distributionMatrix &exact){
  const int Lt = exact.size();
  const int nsample = sloppy(0,0).size();
  distributionVector out(Lt,distributionD(nsample,0.));

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

distributionMatrix timeReflect(const distributionMatrix &m){
  const int Lt = m.size();
  distributionMatrix out(Lt);
  for(int tsrc=0;tsrc<Lt;tsrc++)
    for(int tsep=0;tsep<Lt;tsep++)
      out(tsrc,Lt-tsep-1) = m(tsrc,tsep);
  return out;
}
distributionVector sourceTimeSliceAverage(const distributionMatrix &m){
  const int Lt = m.size();
  const int nsample = m(0,0).size();
  distributionVector out(Lt,distributionD(nsample));
  for(int tsrc=0;tsrc<Lt;tsrc++)
    for(int tsep=0;tsep<Lt;tsep++)
      out(tsep) = out(tsep) + m(tsrc,tsep);
    
  return out/double(Lt);
}

class filterCoordTrange{
  double t_min;
  double t_max;
  bool _invert;
public:
  filterCoordTrange(const double _t_min, const double _t_max): t_min(_t_min), t_max(_t_max), _invert(false){}

  void invert(){ _invert = !_invert; } 
  
  template<typename T>
  bool accept(const Coord &x, const T &y) const{
    bool cond = x.t >= t_min && x.t <= t_max;
    return _invert ? !cond : cond;
  }
};


#endif
