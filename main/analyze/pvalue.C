//Compute the p-value from a hdf5-stored chi^2 or chi^2/dof

#include<common.h>
#include<utils.h>
#include<fit/pvalue.h>

using namespace SARLaC;

int main(const int argc, const char* argv[]){
  assert(argc >= 3);
  int i=1;
  const std::string file = argv[i++];
  double dof = strToAny<double>(argv[i++]);

  bool is_chisq_per_dof = false;
  bool write_out = false;
  std::string outfile;

  enum WhichDistribution { ChiSquare, Tsquare }; //choose between \chi^2 or Hotelling's T^2 distribution
  WhichDistribution dist = ChiSquare;

  int Nfix = -1;

  while(i<argc){
    std::string arg(argv[i]);
    if(arg == "-chisq_per_dof"){ //input distribution is chi^2/dof, not chi^2
      is_chisq_per_dof = true;
      i++;
    }else if(arg == "-o"){ //output pvalue to file
      write_out = true;
      outfile = argv[i+1];
      i+=2;
    }else if(arg == "-infer_dof"){ //provide chi^2 and chi^2/dof as separate files, from which we infer dof
      jackknifeDistributionD chisq, chisq_per_dof;
      readParamsStandard(chisq, argv[i+1]);
      readParamsStandard(chisq_per_dof, argv[i+2]);
      i+=3;
      dof = chisq.best()/chisq_per_dof.best();
      std::cout << "Inferred dof = " << dof << std::endl;
    }else if(arg == "-distribution"){
      std::string dtype(argv[i+1]); 
      if(dtype == "ChiSquare"){
	dist = ChiSquare;
      }else if(dtype == "Tsquare"){
	dist = Tsquare;
      }else error_exit(std::cout << "Unknown distribution " << dtype << std::endl);
      i+=2;
    }else if(arg == "-override_N"){ //fix the value for the number of samples N used for the T^2 distribution - defaults to the number of jackknife samples
      Nfix = strToAny<int>(argv[i+1]);
      assert(Nfix > 0);
      i+=2;
    }else{
      error_exit(std::cout << "Unrecognized argument: " << arg << std::endl);
    }
  }
  jackknifeDistributionD chisq;
  readParamsStandard(chisq, file);

  if(is_chisq_per_dof) chisq = chisq * double(dof);
  
  int nsample = chisq.size();

  int N = Nfix == -1 ? nsample : Nfix;

  if(dist == Tsquare) std::cout << "Using T^2 distribution with N=" << N << " (Nsample = " << nsample << ")" << std::endl;

  jackknifeDistributionD pvalue(nsample, [&](const int s){ 
      return dist == ChiSquare ? 
	chiSquareDistribution::pvalue(dof, chisq.sample(s)) :
	TsquareDistribution::pvalue(chisq.sample(s), dof, N-1);	
    });
  
  std::string dnm = dist == ChiSquare ? "chi^2" : "T^2";
  
  std::cout << "chi^2 = " << chisq << " dof = " << dof << " p(" << dnm << ") = " << pvalue << std::endl;

  if(write_out) writeParamsStandard(pvalue, outfile);

  return 0;
}
