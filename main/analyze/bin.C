//Read a rawDataDistribution and bin the value

#include<distribution.h>
#include<fit/fit_wrapper/fit_wrapper_freeze.h>
#include<common.h>
#include<utils.h>

using namespace CPSfit;

std::vector<int> parseIndices(const std::string &ind){
  std::vector<int> indices;
  std::istringstream buffer(ind);
  std::string line;
  while(std::getline(buffer, line, ',')){
    int i;
    std::stringstream ss; ss << line; ss >> i;
    indices.push_back(i);
  }
  return indices;
}

int main(const int argc, const char** argv){
  if(argc != 5){
    std::cout << "Usage bin <infile> <index string> <bin size> <outfile>" << std::endl;
    return 0;
  }

  std::string infile = argv[1];
  std::string idx_str = argv[2];
  int bin_size = strToAny<int>(argv[3]);
  std::string outfile = argv[4];
  
  assert(fileExists(infile));
  assert(infile != outfile);
  std::vector<int> idx = parseIndices(idx_str);
  
  rawDataDistributionD A;
  readHDF5file(A,infile,idx);

  std::cout << "Read " << A << std::endl;

  A = A.bin(bin_size, true);

  std::cout << "Binned value " << A << std::endl;

  writeParamsStandard(A, outfile);
  std::cout << "Done" << std::endl;
  
  return 0;
}
