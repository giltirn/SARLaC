//Generate a jackknife distribution with a chosen central value and error using a random gaussian

#include<common.h>
#include<utils.h>
#include<random.h>

using namespace CPSfit;


jackknifeDistributionD randomJackknife(double cen, double err, double err_tol, int N, int iter_max = 1000, int iter = 0){
  double dev = err / sqrt(double(N-1));
  jackknifeDistributionD out(N);
  for(int i=0;i<N;i++)
    out.sample(i) = gaussianRandom<double>(cen, dev);

  //Shift slightly to correct mean
  double delta = cen - out.mean();
  for(int i=0;i<N;i++)
    out.sample(i) += delta;

  //Make it recursive to get the error correct to within tol
  if( abs(  (out.standardError() - err) / err  ) > err_tol ){
    if(iter >= iter_max){
      throw "Error: Could not generate distribution with desired tolerance";
    }
    return randomJackknife(cen, err, err_tol, N, iter_max, iter+1);
  }
  return out;
}



int main(const int argc, const char* argv[]){
  assert(argc >= 5);
  int i=1;
  double cen = strToAny<double>(argv[i++]);
  double std_err = strToAny<double>(argv[i++]);
  int nsample = strToAny<int>(argv[i++]);
  std::string outfile = argv[i++];

  double err_tol = 1e-3;
  int ntry = 10000;
  
  int seed = 4967;

  while(i<argc){
    std::string arg(argv[i]);
    if(arg == "-cen_tol"){ 
      std::cout << "Warning: cen_tol is deprecated" << std::endl;
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

  jackknifeDistributionD out(nsample, cen);
  if(std_err != 0.){
    try{
      out = randomJackknife(cen, std_err, err_tol, nsample, ntry);
    }catch(const char* e){
      std::cerr << "Error: " << e << std::endl;
      return 1;
    }
    std::cout << out << std::endl;

    writeParamsStandard(out, outfile);
    return 0;
  }else{
    writeParamsStandard(out, outfile);
    return 0;
  }
}
