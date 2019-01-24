#include<fit/fit_wrapper/fit_wrapper_freeze.h>
#include<physics/luscherzeta.h>
#include<physics/delta_0_pheno.h>

using namespace CPSfit;

//Parse a space-separated set of integers
std::vector<int> parseIntVector(const std::string &from){
  std::stringstream ss;
  ss << from;
  int i;
  std::vector<int> out;
  while(ss >> i){
    out.push_back(i);
  }
  return out;
}
template<typename T, int N>
std::array<T,N> parseArray(const std::string &from){
  std::stringstream ss;
  ss << from;
  T i;
  std::array<T,N> out;
  int idx=0;
  while(ss >> i){  
    assert(idx != N);
    out[idx++] = i;
  }
  return out;
}



int main(const int argc, const char* argv[]){
#ifdef HAVE_HDF5
 
  if(argc < 8){
    std::cout << "Usage: <exe> <Epipi file> <Epipi idx> <Epi file> <Epi idx> <twists> <L> <output file>  [options]\n";
    std::cout << "Where:\n";
    std::cout << "\t<... idx> are array indices for the files in question. If the file contains a multi-dimensional array the separate indices should be space separated and enclosed in quotation marks, e.g. \"0 1\"\n";
    std::cout << "\t<twists> should be a space separated set of 3 integers in quotation marks with value 0 for periodic directions and 1 for antiperiodic/G-parity, eg. \"1 1 1\"\n";
    std::cout << "Options:\n";
    std::cout << "\t-comoving <d>    Epipi has non-zero total momentum. <d> should be \\vec P_CM in units of 2pi/L\n";
    std::cout << "\t-schenk <I> Print the Schenk values (3 curves A,B,C). <I> is the isospin\n";
    std::cout << "\t-colangelo <I>  Print the Colangelo value. <I> is the isospin\n";
    std::cout << "\t-write_schenk <filename> / -write_colangelo <filename>   Write the Schenk / Colangelo numbers to <filename>. Use in conjunction with the above\n";
    exit(0);
  }

  std::string Epipi_file = argv[1];
  std::vector<int> Epipi_idx = parseIntVector(argv[2]);
  std::string Epi_file = argv[3];
  std::vector<int> Epi_idx = parseIntVector(argv[4]);
  std::array<int,3> twists = parseArray<int,3>(argv[5]);
  int L = strToAny<int>(argv[6]);
  std::string delta_file = argv[7];

  if(twists.size() != 3) error_exit(std::cout << "Could not parse \"" << argv[5] << "\" into a vector of size 3, got " << twists << std::endl);

  std::array<double,3> d = {0.,0.,0.}; //Pcm in units of 2pi/L

  int isospin = 0;
  bool print_schenk = false;
  bool print_colangelo = false;
  bool write_schenk = false;
  bool write_colangelo = false;
  std::string write_schenk_file;
  std::string write_colangelo_file;
  
  int ii=8;
  while(ii < argc){
    if(std::string(argv[ii]) == "-comoving"){
      d = parseArray<double,3>(argv[ii+1]);
      ii+=2;
    }else if(std::string(argv[ii]) == "-schenk"){
      print_schenk = true;
      isospin = strToAny<int>(argv[ii+1]);
      ii+=2;
    }else if(std::string(argv[ii]) == "-colangelo"){
      print_colangelo = true;
      isospin = strToAny<int>(argv[ii+1]);
      ii+=2;
    }else if(std::string(argv[ii]) == "-write_schenk"){
      write_schenk = true;
      write_schenk_file = argv[ii+1];
      ii+=2;
    }else if(std::string(argv[ii]) == "-write_colangelo"){
      write_colangelo = true;
      write_colangelo_file = argv[ii+1];
      ii+=2;
    }else{
      error_exit(std::cout << "Unrecognized argument: " << argv[ii] << std::endl);      
    }
  }


  jackknifeDistribution<double> Epipi;
  readHDF5file(Epipi, Epipi_file,Epipi_idx);

  jackknifeDistribution<double> Epi;
  readHDF5file(Epi, Epi_file,Epi_idx);

  assert(Epipi.size() == Epi.size());
  int nsample = Epipi.size();

  int ntwist = twists[0] + twists[1] + twists[2];

  jackknifeDistribution<double> mpi = Epi;

  if(ntwist > 0){
    double p2 = ntwist * pow(M_PI/L,2);
    mpi = sqrt( Epi*Epi - p2 );
  }
  
  std::cout << "Epipi = " << Epipi << std::endl;
  std::cout << "Epi = " << Epi << std::endl;
  std::cout << "mpi = " << mpi << std::endl;

  LuscherZeta zeta(twists,d);
  jackknifeDistribution<double> delta(nsample);
#pragma omp parallel for
  for(int s=0;s<nsample;s++)
    delta.sample(s) = phaseShiftZ(zeta, Epipi.sample(s),mpi.sample(s),L);

  std::cout << "delta = " << delta << std::endl;

  writeParamsStandard(delta, delta_file);


  if(print_schenk || print_colangelo){
    const GSLvector &d = zeta.getd();
    double _2pidL = 2*M_PI/L;
    double Pcm2 = _2pidL*_2pidL* d.norm2();
    
    jackknifeDistribution<double> gamma(nsample, [&](const int s){ return Epipi.sample(s)/sqrt(Epipi.sample(s)*Epipi.sample(s) - Pcm2);  });    
    jackknifeDistribution<double> Epipi_CM = Epipi/gamma;
    jackknifeDistribution<double> S = Epipi_CM*Epipi_CM;
    
    if(print_schenk){
      char curves[3] = {'A','B','C'};
      std::vector<jackknifeDistribution<double> > delta(3);
      for(int c=0;c<3;c++){
	delta[c] = jackknifeDistribution<double>(nsample, [&](const int s){ return PhenoCurveSchenk::compute(S.sample(s), isospin, mpi.sample(s), curves[c]);});
	delta[c] = delta[c] * 180./M_PI;
      }
      std::cout << "Schenk value for I=" << isospin << " = " << delta[0] << " (A) " << delta[1] << " (B) " << delta[2] << " (C)" << std::endl;      

      if(write_schenk) writeParamsStandard(delta, write_schenk_file);
    }
    if(print_colangelo){
      jackknifeDistribution<double> delta(nsample, [&](const int s){ return PhenoCurveColangelo::compute(S.sample(s), isospin, mpi.sample(s));});
      delta = delta * 180./M_PI;

      std::cout << "Colangelo value for I=" << isospin << " = " << delta << std::endl;

      if(write_colangelo) writeParamsStandard(delta, write_colangelo_file);
    }
  }
#else
  error_exit(std::cout << "Require HDF5\n");
#endif
  return 0;
}
