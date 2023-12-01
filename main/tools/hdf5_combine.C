#include<distribution.h>
#include<utils.h>

using namespace SARLaC;


//Append conventionally-stored distributions from fit parameter files to an output file, either new or existing. Can be used to create new parameter files from multiple fit outputs for, say, freezing parameters
//Expects output files to be vector of distributions (not vector-vector); input file can be either

template<typename DistributionType>
void runD(const std::string &infile, const std::string &outfile, const int vector_depth_in, bool outfile_exists, const int argc, const char** argv){
  //Extract the data we want from the input file
  std::vector<DistributionType> out;
  if(outfile_exists)
    readParamsStandard(out, outfile);

  if(vector_depth_in == 1){
    std::vector<DistributionType> in;
    readParamsStandard(in, infile);
    for(int i=3;i<argc;i++){
      int in_idx = strToAny<int>(argv[i]);
      assert(in_idx >=0 && in_idx < in.size());
      out.push_back(in[in_idx]);
    }
  }else if(vector_depth_in == 2){
    std::vector<std::vector<DistributionType> > in;
    readParamsStandard(in, infile);
    for(int i=3;i<argc;i+=2){
      int in_idx0 = strToAny<int>(argv[i]);
      int in_idx1 = strToAny<int>(argv[i+1]);
      assert(in_idx0 >=0 && in_idx0 < in.size());
      assert(in_idx1 >=0 && in_idx1 < in[in_idx0].size());
      out.push_back(in[in_idx0][in_idx1]);
    }
  }
  writeParamsStandard(out, outfile);
}

void run(DistributionTypeEnum type, const std::string &infile, const std::string &outfile, const int vector_depth_in, bool outfile_exists, const int argc, const char** argv){
  switch(type){
  case DistributionTypeEnum::Raw:
    runD<rawDataDistribution<double> >(infile,outfile,vector_depth_in,outfile_exists,argc,argv);  break;
  case DistributionTypeEnum::Jackknife:
    runD<jackknifeDistribution<double> >(infile,outfile,vector_depth_in,outfile_exists,argc,argv);  break;
  case DistributionTypeEnum::JackknifeC:
    runD<jackknifeCdistribution<double> >(infile,outfile,vector_depth_in,outfile_exists,argc,argv);  break;
  case DistributionTypeEnum::DoubleJackknife:
    runD<doubleJackknifeDistribution<double> >(infile,outfile,vector_depth_in,outfile_exists,argc,argv);  break;
  case DistributionTypeEnum::SuperJackknife:
    runD<superJackknifeDistribution<double> >(infile,outfile,vector_depth_in,outfile_exists,argc,argv);  break;    
  case DistributionTypeEnum::Bootstrap:
    runD<bootstrapDistribution<double> >(infile,outfile,vector_depth_in,outfile_exists,argc,argv);  break;
  default:
    error_exit(std::cout << "run(...) unknown type " << type << std::endl);
  }
}



int main(const int argc, const char** argv){
  assert(argc >= 4);
  std::string infile = argv[1];
  std::string outfile = argv[2];
  
  assert(fileExists(infile));
  bool outfile_exists = fileExists(outfile);

  DistributionTypeEnum type_in;
  int vector_depth_in;
  getTypeInfo(type_in,vector_depth_in,infile);
  
  if(outfile_exists){
    DistributionTypeEnum type_out;
    int vector_depth_out;
    getTypeInfo(type_out,vector_depth_out,outfile);
    assert(vector_depth_out == 1);
    assert(type_out == type_in);
  }
  run(type_in,infile, outfile, vector_depth_in, outfile_exists, argc, argv);

  return 0;
}
