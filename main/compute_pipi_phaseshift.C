#include<fit/fit_wrapper/fit_wrapper_freeze.h>
#include<physics/luscherzeta.h>

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


int main(const int argc, const char* argv[]){
#ifdef HAVE_HDF5

  assert(argc == 8);
  std::string Epipi_file = argv[1];
  std::vector<int> Epipi_idx = parseIntVector(argv[2]);
  std::string Epi_file = argv[3];
  std::vector<int> Epi_idx = parseIntVector(argv[4]);
  std::vector<int> twists = parseIntVector(argv[5]);
  int L = strToAny<int>(argv[6]);
  std::string delta_file = argv[7];

  if(twists.size() != 3) error_exit(std::cout << "Could not parse \"" << argv[5] << "\" into a vector of size 3, got " << twists << std::endl);

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

  LuscherZeta zeta(twists[0],twists[1],twists[2]);
  jackknifeDistribution<double> delta(nsample, [&](const int s){ return phaseShiftZ(Epipi.sample(s),mpi.sample(s),L,zeta); });

  std::cout << "delta = " << delta << std::endl;

  writeParamsStandard(delta, delta_file);

#else
  error_exit(std::cout << "Require HDF5\n");
#endif
  return 0;
}
