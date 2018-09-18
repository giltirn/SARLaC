#ifndef _FIT_PIPI_GND_EXC_GPARITY_CMDLINE_H
#define _FIT_PIPI_GND_EXC_GPARITY_CMDLINE_H

struct CMDline{
  bool load_guess;
  std::string guess_file;

  bool load_hdf5_data_checkpoint;
  std::string load_hdf5_data_checkpoint_stub;
  bool save_hdf5_data_checkpoint;
  std::string save_hdf5_data_checkpoint_stub;

  bool load_combined_data;
  std::string load_combined_data_file;

  bool save_combined_data;
  std::string save_combined_data_file;

  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  CMDline(){
    load_guess = false;
    load_hdf5_data_checkpoint = false;
    save_hdf5_data_checkpoint = false;
    load_combined_data = false;
    save_combined_data = false;
    load_frozen_fit_params = false;
  }
  CMDline(const int argc, const char** argv, const int begin = 0): CMDline(){
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
      }else if(sargv[i] == "-load_hdf5_data_checkpoint"){
	load_hdf5_data_checkpoint = true;
	load_hdf5_data_checkpoint_stub = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_hdf5_data_checkpoint"){
	save_hdf5_data_checkpoint = true;
	save_hdf5_data_checkpoint_stub = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_combined_data"){ //load the double-jackknife data set previously generated
	load_combined_data = true;
	load_combined_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_combined_data"){ //save the double-jackknife data set previously generated
	save_combined_data = true;
	save_combined_data_file = sargv[i+1];
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
	}else if(!fileExists(load_frozen_fit_params_file)) error_exit(std::cout << "CMDline freeze data file " << load_frozen_fit_params_file << " does not exist!\n");
	i+=2;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};


#endif
