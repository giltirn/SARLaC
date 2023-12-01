//Perform a binary operation on a pair of distributions

#include<config.h>

#ifdef HAVE_HDF5
#include<iostream>
#include<cassert>
#include<sstream>
#include<array>

#include<distribution.h>
#include<parser.h>
#include<common.h>

using namespace CPSfit;

std::vector<int> parseElemStr(const std::string &e){
  std::istringstream buffer(e);
  std::string line;
  std::vector<std::string> split;
  while (std::getline(buffer, line, ',')){
    split.push_back(line);
  }

  std::vector<int> out(split.size());
  for(int i=0;i<split.size();i++){
    std::stringstream ss; ss << split[i]; ss >> out[i];
  }
  return out;
}
  
  
struct Args{
  std::string filename[2];
  std::vector<int> elem[2];  
  std::string operation;

  DistributionTypeEnum type;
  int vector_depth[2];

  Args(const int argc, const char* argv[], const int start = 1){
    filename[0] = argv[start];
    elem[0] = parseElemStr(argv[start+1]);
    filename[1] = argv[start+2];
    elem[1] = parseElemStr(argv[start+3]);
    operation = argv[start+4];

    DistributionTypeEnum types[2];
    getTypeInfo(types[0],vector_depth[0],filename[0]);
    getTypeInfo(types[1],vector_depth[1],filename[1]);
    
    if(types[0] != types[1]) error_exit(std::cout << "Error: distribution types must be the same\n");

    type = types[0];

    for(int i=0;i<2;i++)
      if(vector_depth[i] != elem[i].size()) error_exit(std::cout << "Error: Expected a list of " << vector_depth[i] << " comma separated ints for the elem index of distribution " << i << std::endl);    
  }
};
  
template<typename Dist>
Dist getData(const std::string &filename, const std::vector<int> &elem, const int vector_depth){
  if(vector_depth == 1){
    std::vector<Dist> out;
    readParamsStandard(out,filename);
    return std::move(out[elem[0]]);
  }else if(vector_depth == 2){
    std::vector<std::vector<Dist> > out;
    readParamsStandard(out,filename);
    return std::move(out[elem[0]][elem[1]]);
  }else error_exit(std::cout << "Unsupported vector depth " << vector_depth << std::endl);
}

template<typename Dist>
struct Correlation{
  static void apply(const Dist &a, const Dist &b){
    auto covab = Dist::covariance(a,b);
    auto covaa = Dist::covariance(a,a);
    auto covbb = Dist::covariance(b,b);
    
    auto corr = covab/sqrt(covaa)/sqrt(covbb);

    std::cout << "Correlation of " << a << " , " << b << " : " << corr << std::endl;
  }
};
template<>
struct Correlation<superJackknifeDistribution<double> >{
  typedef superJackknifeDistribution<double> Dist;
  static void apply(const Dist &a, const Dist &b){
    const superJackknifeLayout &layout = a.getLayout();
    assert(layout == b.getLayout());

    for(int i=0;i<layout.nEnsembles();i++){
      std::cout << "For sub-ensemble " << layout.ensTag(i) << " ";
      Correlation<jackknifeDistribution<double> >::apply(a.getEnsembleJackknife(i), b.getEnsembleJackknife(i));
    }
  }
};

template<typename Dist>
void run(const Args &args){
  Dist a = getData<Dist>(args.filename[0],args.elem[0],args.vector_depth[0]);
  Dist b = getData<Dist>(args.filename[1],args.elem[1],args.vector_depth[1]);

  if(args.operation == "correlation"){
    Correlation<Dist>::apply(a,b);
  }else error_exit(std::cout << "Unsupported operation " << args.operation << std::endl);  
}

  
int main(const int argc, const char* argv[]){
  assert(argc >= 6);
  Args args(argc,argv,1);

  switch(args.type){
  case DistributionTypeEnum::Jackknife:
    run<jackknifeDistributionD>(args); break;
  case DistributionTypeEnum::JackknifeC:
    run<jackknifeCdistributionD>(args); break;
  case DistributionTypeEnum::Raw:
    run<rawDataDistributionD>(args); break;
  case DistributionTypeEnum::SuperJackknife:
    run<superJackknifeDistribution<double> >(args); break;
  default:
    error_exit(std::cout << "Unsupported distribution type " << args.type << std::endl);
  }
    
  return 0;
}


#else

int main(const int argc, const char* argv[]){
  std::cout << "Require HDF5\n";
  return 1;
}

#endif
