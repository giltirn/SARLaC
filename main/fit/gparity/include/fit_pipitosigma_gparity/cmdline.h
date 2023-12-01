#ifndef _PIPI_TO_SIGMA_CMDLINE_H_
#define _PIPI_TO_SIGMA_CMDLINE_H_

struct PiPiToSigmaCMDline{
  bool load_guess;
  std::string guess_file;

  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;
 
  bool force_disconn_tstep_src; //the code computes the disconnected component for all tsrc, but this option can be used to constrain the number of source timeslices to observe the effect
  int disconn_tstep_src;

  bool use_disconn_complex_prod; //use Re( pipi_bubble * sigma_bubble )  [original strategy] for disconnected piece rather than Re ( pipi_bubble ) * Re ( sigma_bubble )

  PiPiToSigmaCMDline(){
    load_guess = false;
    load_frozen_fit_params= false;
    force_disconn_tstep_src = false;
    use_disconn_complex_prod = false;
  }
  PiPiToSigmaCMDline(const int argc, const char** argv, const int begin = 0): PiPiToSigmaCMDline(){
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
      }else if(sargv[i] == "-force_disconn_tstep_src"){
	force_disconn_tstep_src = true;
	disconn_tstep_src = strToAny<int>(sargv[i+1]);
	i+=2;
      }else if(sargv[i] == "-use_disconn_complex_prod"){
	use_disconn_complex_prod = true;
	i++;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }

  void transfer(PiPiToSigmaFitArgs &fargs) const{
    fargs.load_guess = load_guess;
    fargs.guess_file = guess_file;
    fargs.load_frozen_fit_params = load_frozen_fit_params;
    fargs.load_frozen_fit_params_file = load_frozen_fit_params_file;
  }
};


#endif
