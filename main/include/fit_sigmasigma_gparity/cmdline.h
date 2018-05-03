#ifndef _SIGMA_CMDLINE_H_
#define _SIGMA_CMDLINE_H_

struct SigmaCMDline{
  bool load_guess;
  std::string guess_file;

  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;
 
  SigmaCMDline(){
    load_guess = false;
    load_frozen_fit_params= false;
  }
  SigmaCMDline(const int argc, const char** argv, const int begin = 0): SigmaCMDline(){
    setup(argc,argv,begin);
  }
  
  void setup(const int argc, const char** argv, const int begin = 0){
    const int sz = argc-begin;
    std::vector<std::string> sargv(sz);
    for(int i=begin; i<argc; i++) sargv[i-begin] = std::string(argv[i]);

    int i = 0;
    while(i<sz){
      if(sargv[i] == "-load_guess"){
	load_guess = true;
	guess_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-nthread"){
	omp_set_num_threads(strToAny<int>(sargv[i+1]));
	i+=2;
      }else if(sargv[i] == "-load_frozen_fit_params"){
	load_frozen_fit_params = true;
	load_frozen_fit_params_file = sargv[i+1];
	i+=2;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }

  void transfer(SigmaFitArgs &fargs) const{
    fargs.load_guess = load_guess;
    fargs.guess_file = guess_file;
    fargs.load_frozen_fit_params = load_frozen_fit_params;
    fargs.load_frozen_fit_params_file = load_frozen_fit_params_file;
  }
};


#endif
