#ifndef _FIT_MPI_GPARITY_AMA_READ_DATA_H
#define _FIT_MPI_GPARITY_AMA_READ_DATA_H

void read(rawDataDistributionMatrix &exact, rawDataDistributionMatrix &sloppy, const SloppyExact &args, const int traj, const int sample_idx){
  if(args.include_data){
    read(exact, sample_idx, args.exact_fmt, traj, args.reim);
    read(sloppy, sample_idx, args.sloppy_fmt, traj, args.reim);
  }
}

rawDataDistributionVector readCombine(const Args &args, const int type_idx){
  const int ntraj = (args.traj_lessthan - args.traj_start)/args.traj_inc;
  assert(ntraj > 0);

  TwoPointFunction const & fargs = args.data[type_idx];
  const std::string& type = fargs.type;
  FitFuncType fitfunc = fargs.fitfunc;
  
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
      if(fitfunc == FSinh){ //sinh-form data pick up - sign under reflection
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

#endif
