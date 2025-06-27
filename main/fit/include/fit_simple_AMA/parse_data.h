#pragma once
#include "data_manipulations.h"

GENERATE_ENUM_AND_PARSER(ReIm, (Real)(Imaginary) );

GENERATE_ENUM_AND_PARSER(AMAparserType, (ParserComplexLabeled)(ParserReal) );

struct AMAparser{
  int Lt;

  AMAparser(int Lt): Lt(Lt){}
  virtual void parse(rawDataDistributionMatrix &into, std::istream &is, const int sample, const std::string &file) const = 0;  
};

struct AMAparserComplexLabeled: public AMAparser{ //Expect Lt lines of format  <tsrc> <tsep> <re> <im>
  ReIm reim;
  AMAparserComplexLabeled(int Lt, const ReIm reim): AMAparser(Lt), reim(reim){}
  void parse(rawDataDistributionMatrix &into, std::istream &is, const int sample, const std::string &file) const{
    const int nelems = Lt*Lt;
    int i,j;
    double re, im;
    for(int e=0;e<nelems;e++){
      if(!(is >> i >> j)) error_exit(std::cout << "AMAparserComplexLabeled failed to read indices for file " << file << "\n");
      if(!(is >> re >> im)) error_exit(std::cout << "AMAparserComplexLabeled failed to read data for file " << file << "\n");
      into(i,j).sample(sample) = reim == ReIm::Real ? re : im;
    }
    if(is.fail() || is.bad()){ std::cout << "Error reading file \"" << file << "\"\n"; std::cout.flush(); exit(-1); }
  }
};
  
struct AMAparserReal: public AMAparser{ //Expect Lt lines of format  <re tsep=0> <re tsep=1> ...
  AMAparserReal(int Lt): AMAparser(Lt){}
  void parse(rawDataDistributionMatrix &into, std::istream &is, const int sample, const std::string &file) const{
    double re;
    for(int tsrc=0;tsrc<Lt;tsrc++){
      for(int tsep=0;tsep<Lt;tsep++){
	if(!(is >> re)) error_exit(std::cout << "AMAparserReal failed to read data for file " << file << " tsrc " << tsrc << " tsep " << tsep << std::endl);
	into(tsrc,tsep).sample(sample) = re;	
      }
    }
    if(is.fail() || is.bad()){ std::cout << "Error reading file \"" << file << "\"\n"; std::cout.flush(); exit(-1); }

    // std::cout << "AMAparserReal read " << file << std::endl;
    // for(int tsrc=0;tsrc<Lt;tsrc++){
    //   for(int tsep=0;tsep<Lt;tsep++){
    // 	std::cout << tsrc << " " << tsep << " " << into(tsrc,tsep) << std::endl;
    //   }
    // }
  }
};
  
  
AMAparser* AMAparserFactory(const AMAparserType p, int Lt, const ReIm reim){
  switch(p){
  case AMAparserType::ParserComplexLabeled:
    return new AMAparserComplexLabeled(Lt,reim);
  case AMAparserType::ParserReal:
    return new AMAparserReal(Lt);
  default:
    error_exit(std::cout << "AMAparserFactory: Unknown parser " << p << std::endl);
  };
}


void read(rawDataDistributionMatrix &into, const int sample, const std::string &fmt, const int traj, const ReIm reim, const AMAparserType parser){
  const int Lt = into.size();

  const std::string file = subsIdx(fmt,traj);
  std::ifstream is(file.c_str());
  if(!is.good()){ std::cout << "Could not open file \"" << file << "\"\n"; std::cout.flush(); exit(-1); }

  std::unique_ptr<AMAparser> p(AMAparserFactory(parser,Lt,reim));
  p->parse(into, is, sample, file);
  is.close();
}

std::vector<int> readMultiplicity(int traj, const std::string &fmt, int Lt, int nexact_timeslice){
  subStringReplace repl(fmt,{subStringSpecify("%d")});
  std::string filename = repl.replace({ std::to_string(traj) });

  std::cout << "Parsing multiplicity " << filename << std::endl;
  std::ifstream is(filename.c_str());
  if(is.good()){
    std::vector<int> out(Lt);
    int sum = 0;
    for(int t=0;t<Lt;t++){
      is >> out[t];
      assert(!is.fail());
      sum += out[t];
    }
    if(sum != nexact_timeslice){
      error_exit(std::cout << "readMultiplicity sum of multiplicities " << sum << " does not add to nexact_timeslice=" << nexact_timeslice << std::endl);
    }
    return out;
  }else{
    error_exit(std::cout << "readMultiplicity failed to read file " << filename << std::endl);
  }
}

rawDataDistributionVector readData(const int traj_start, const int traj_inc, const int traj_lessthan,
				   const std::string &sloppy_fmt, const std::string &exact_fmt, 
				   const bool have_multiplicity, const std::string &mult_file_fmt,
				   const int nexact_timeslice,
				   const ReIm reim, const int Lt, AMAparserType parser){
  const int ntraj = (traj_lessthan - traj_start)/traj_inc;
  assert(ntraj > 0);

  rawDataDistributionVector corrected;
    
  rawDataDistributionMatrix exact_data(Lt, rawDataDistributionD(ntraj));
  rawDataDistributionMatrix sloppy_data(Lt, rawDataDistributionD(ntraj));
  std::vector<std::vector<int> > multiplicity(ntraj);

#pragma omp parallel for
  for(int i=0;i<ntraj;i++){
    const int c = traj_start + i*traj_inc;
    read(exact_data, i, exact_fmt, c, reim, parser);
    read(sloppy_data, i, sloppy_fmt, c, reim, parser);

    if(have_multiplicity) multiplicity[i] = readMultiplicity(c,mult_file_fmt,Lt,nexact_timeslice);
    else multiplicity[i] = randomMultiplicity(exact_data, c, i, nexact_timeslice);
  }
  
  rawDataDistributionVector sloppy_avg = sourceTimeSliceAverage(sloppy_data);
  rawDataDistributionVector correction = computeAMAcorrection(sloppy_data, exact_data, multiplicity, nexact_timeslice);
  
  corrected = sloppy_avg + correction;

  std::cout << "Sloppy data:\n";
  for(int tsrc=0;tsrc<Lt;tsrc++)
    for(int t=0;t<Lt;t++)
      std::cout << tsrc << " " << t << " " << sloppy_data(tsrc,t) << std::endl;

  std::cout << "Sloppy timeslice averaged data:\n";
  for(int t=0;t<Lt;t++) std::cout << t << " " << sloppy_avg[t] << std::endl;
  
  std::cout << "Corrected timeslice averaged data:\n";
  for(int t=0;t<Lt;t++) std::cout << t << " " << corrected[t] << std::endl;
  
  return corrected;
}
