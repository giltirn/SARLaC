#ifndef _PIPI_SIGMA_SIM_FIT_CMDLINE_H_
#define _PIPI_SIGMA_SIM_FIT_CMDLINE_H_

struct PiPiSigmaSimCMDline{
  bool load_guess;
  std::string guess_file;

  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;
 
  bool include_pipi_2pt;
  bool include_pipi_to_sigma;
  bool include_sigma_2pt;

  bool save_checkpoint;
  std::string save_checkpoint_file;
  
  bool load_checkpoint;
  std::string load_checkpoint_file;
  

  PiPiSigmaSimCMDline(){
    load_guess = false;
    load_frozen_fit_params= false;
    include_pipi_2pt = true;
    include_pipi_to_sigma = true;
    include_sigma_2pt = true;
    save_checkpoint = false;
    load_checkpoint = false;
  }
  PiPiSigmaSimCMDline(const int argc, const char** argv, const int begin = 0): PiPiSigmaSimCMDline(){
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
	if(load_frozen_fit_params_file == "TEMPLATE"){
	  std::cout << "Saving frozen fit params template file to freeze_template.args" << std::endl;
	  FreezeParams fp;
	  std::ofstream of("freeze_template.args");
	  of << fp;
	  of.close();
	  exit(0);
	}
	i+=2;
      }else if(sargv[i] == "-exclude_pipi_2pt"){
	include_pipi_2pt = false;
	i++;
      }else if(sargv[i] == "-exclude_pipi_to_sigma"){
	include_pipi_to_sigma = false;
	i++;
      }else if(sargv[i] == "-exclude_sigma_2pt"){
	include_sigma_2pt = false;
	i++;
      }else if(sargv[i] == "-save_checkpoint"){
	save_checkpoint = true;
	save_checkpoint_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_checkpoint"){
	load_checkpoint = true;
	load_checkpoint_file = sargv[i+1];
	i+=2;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }

  void transfer(SimFitArgs &fargs) const{
    fargs.load_guess = load_guess;
    fargs.guess_file = guess_file;
    fargs.load_frozen_fit_params = load_frozen_fit_params;
    fargs.load_frozen_fit_params_file = load_frozen_fit_params_file;
  }
};


#endif
