#ifndef _FIT_MPI_GPARITY_AMA_DATA_MANIPULATIONS_H
#define _FIT_MPI_GPARITY_AMA_DATA_MANIPULATIONS_H

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

rawDataDistributionVector readCombine(const Args &args, const int type_idx){
  const int ntraj = (args.traj_lessthan - args.traj_start)/args.traj_inc;
  assert(ntraj > 0);

  TwoPointFunction const & fargs = args.data[type_idx];
  DataType type = fargs.type;

  rawDataDistributionVector corrected[2];
  int FF=0, BB=1;
  
  for(int fb=0;fb<2;fb++){
    bool include_data = fb == FF ? fargs.FF_data.include_data : fargs.BB_data.include_data;
    if(!include_data) continue;

    const SloppyExact &se = fb == FF ? fargs.FF_data : fargs.BB_data;
    
    rawDataDistributionMatrix exact_data(args.Lt, rawDataDistributionD(ntraj));
    rawDataDistributionMatrix sloppy_data(args.Lt, rawDataDistributionD(ntraj));

#pragma omp parallel for
    for(int i=0;i<ntraj;i++){
      const int c = args.traj_start + i*args.traj_inc;
      read(exact_data, sloppy_data, se, c, i);
    }

    if(fb == BB){
      exact_data = timeReflect(exact_data);
      sloppy_data = timeReflect(sloppy_data);
      if(type == AP_LW_data || type == AP_WW_data){ //sinh-form data pick up - sign under reflection
	exact_data = -exact_data;
	sloppy_data = -sloppy_data;
      }
    }
    
    rawDataDistributionVector sloppy_avg;
    sloppy_avg = sourceTimeSliceAverage(sloppy_data);

    rawDataDistributionVector correction = computeAMAcorrection(sloppy_data, exact_data);

    corrected[fb] = sloppy_avg + correction;

    std::string nm = (fb == FF ? "FF" : "BB");

    std::cout << type << " " << nm << " sloppy data:\n";
    for(int t=0;t<args.Lt;t++) std::cout << t << " " << sloppy_avg[t] << std::endl;

    std::cout << type << " " << nm << " corrected data:\n";
    for(int t=0;t<args.Lt;t++) std::cout << t << " " << corrected[fb][t] << std::endl;
  }
  rawDataDistributionVector out;
  if(fargs.FF_data.include_data && fargs.BB_data.include_data) return (corrected[FF] + corrected[BB])/2.;
  else if(fargs.FF_data.include_data) return corrected[FF];
  else return corrected[BB];
}



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



#endif
