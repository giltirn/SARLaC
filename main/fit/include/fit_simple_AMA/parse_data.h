#pragma once

GENERATE_ENUM_AND_PARSER(ReIm, (Real)(Imaginary) );

typedef NumericSquareMatrix<rawDataDistribution<double> > rawDataDistributionMatrix;
typedef NumericVector<rawDataDistribution<double> > rawDataDistributionVector;

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
