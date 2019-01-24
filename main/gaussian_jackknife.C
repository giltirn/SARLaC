//Generate a jackknife distribution with a chosen central value and error using a random gaussian

#include<common.h>
#include<utils.h>
#include<random.h>

using namespace CPSfit;

int main(const int argc, const char* argv[]){
  assert(argc >= 5);
  int i=1;
  double cen = strToAny<double>(argv[i++]);
  double std_err = strToAny<double>(argv[i++]);
  int nsample = strToAny<int>(argv[i++]);
  std::string outfile = argv[i++];

  double cen_tol = 1e-4;
  double err_tol = 1e-3;
  int ntry = 10000;
  
  int seed = 4967;

  while(i<argc){
    std::string arg(argv[i]);
    if(arg == "-cen_tol"){ 
      cen_tol = strToAny<double>(argv[i+1]);
      i+=2;
    }else if(arg == "-err_tol"){ 
      err_tol = strToAny<double>(argv[i+1]);
      i+=2;
    }else if(arg == "-ntry"){ 
      ntry = strToAny<int>(argv[i+1]);
      i+=2;
    }else if(arg == "-seed"){ 
      seed = strToAny<int>(argv[i+1]);
      i+=2;
    }else{
      error_exit(std::cout << "Unrecognized argument: " << arg << std::endl);
    }
  }
  
  RNG.initialize(seed);

  jackknifeDistributionD out(nsample);
  
  for(int i=0;i<ntry;i++){
    gaussianRandom(out, cen, std_err / sqrt((double)nsample-1.));

    double c = out.mean();
    double e = out.standardError();

    double ce = fabs( (c - cen)/c );
    double ee = fabs( (e - std_err)/e );

    std::cout << out << " " << ce << " " << ee << std::endl;
    if(ce <= cen_tol && ee <= err_tol){
      writeParamsStandard(out, outfile);
      return 0;
    }
  }
  std::cout << "Error: Could not generate distribution with desired tolerance" << std::endl;
  return 1;
}
