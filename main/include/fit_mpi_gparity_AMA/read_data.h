#ifndef _FIT_MPI_GPARITY_AMA_READ_DATA_H
#define _FIT_MPI_GPARITY_AMA_READ_DATA_H

typedef NumericSquareMatrix<rawDataDistribution<double> > rawDataDistributionMatrix;
typedef NumericVector<rawDataDistribution<double> > rawDataDistributionVector;

void read(rawDataDistributionMatrix &into, const int sample, const std::string &fmt, const int traj, const ReIm reim){
  const int Lt = into.size();

  const std::string file = subsIdx(fmt,traj);
  std::ifstream is(file.c_str());
  if(!is.good()){ std::cout << "Could not open file \"" << file << "\"\n"; std::cout.flush(); exit(-1); }

  const int nelems = Lt*Lt;
  int i,j;
  double re, im;
  for(int e=0;e<nelems;e++){
    if(!(is >> i >> j)) error_exit(std::cout << "configData failed to read indices for file " << file << "\n");
    if(!(is >> re >> im)) error_exit(std::cout << "configData failed to read data for file " << file << "\n");
    into(i,j).sample(sample) = reim == Real ? re : im;
  }
  if(is.fail() || is.bad()){ std::cout << "Error reading file \"" << file << "\"\n"; std::cout.flush(); exit(-1); }
  is.close();
}

void read(rawDataDistributionMatrix &exact, rawDataDistributionMatrix &sloppy, const SloppyExact &args, const int traj, const int sample_idx){
  if(args.include_data){
    read(exact, sample_idx, args.exact_fmt, traj, args.reim);
    read(sloppy, sample_idx, args.sloppy_fmt, traj, args.reim);
  }
}





#endif
