#include<distribution.h>
#include<utils.h>
#include<common.h>
#include<fit/fit_wrapper/fit_wrapper_freeze.h>

//Take a regular jackknife distribution and boost it to a "simple" superjackknife distribution - i.e. just a regular jackknife class instance but setup as a superjackknife

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
  if(argc != 6){
    std::cout << "Usage: <exe> <infile> <index string (comma separated)> <num total samples> <sample offset of input in output> <outfile>" << std::endl;
    exit(0);
  } 

  std::string infile = argv[1];
  std::vector<int> infile_idx = parseIndices(argv[2]);
  int nsamples_total = strToAny<int>(argv[3]);
  int offset = strToAny<int>(argv[4]);
  std::string outfile = argv[5];

  assert(fileExists(infile));

  jackknifeDistributionD in;
  readHDF5file(in, infile, infile_idx);

  std::cout << "Input data: " << in << std::endl;

  assert(in.size() < nsamples_total);
  assert(offset + in.size() <= nsamples_total);

  double c = in.mean();
  jackknifeDistributionD out(nsamples_total, c);
  for(int i=offset;i<offset + in.size();i++){
    out.sample(i) = in.sample(i-offset);
  }
  std::cout << "Boosted data: " << out << std::endl;

  writeParamsStandard(out,outfile);

  std::cout << "Done" << std::endl;
  return 0;
}
