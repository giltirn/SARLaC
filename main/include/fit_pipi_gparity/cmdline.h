#ifndef _FIT_PIPI_GPARITY_CMDLINE_H
#define _FIT_PIPI_GPARITY_CMDLINE_H

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

  bool load_mlparams;
  std::string mlparams_file;

  CMDline(){
    load_guess = false;
    load_hdf5_data_checkpoint = false;
    save_hdf5_data_checkpoint = false;
    load_combined_data = false;
    save_combined_data = false;
    load_frozen_fit_params= false;
    load_mlparams = false;
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
	  FreezeParams templ;
	  std::ofstream of("freeze_template.dat");
	  of << templ;
	  of.close();
	  std::cout << "Wrote freeze template to freeze_template.dat\n";
	  exit(0);	 
	}
	i+=2;
      }else if(sargv[i] == "-load_mlparams"){
	load_mlparams = true;
	mlparams_file = sargv[i+1];
	if(mlparams_file == "TEMPLATE"){
	  MarquardtLevenbergParameters<double> templ;
	  std::ofstream of("mlparams_template.args");
	  of << templ;
	  of.close();
	  std::cout << "Wrote MLparams template to mlparams_template.args\n";
	  exit(0);	 
	}
	i+=2;
      }else if(sargv[i] == "-allow_bin_cropping"){ //when binning allow extra configs that don't fit a full bin to be cropped from the end of the ensemble
	rawDataDistributionOptions::binAllowCropByDefault() = true;
	i++;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }


  void Export(pipiFitOptions &to){
#define I(A) to.A = A
    I(load_frozen_fit_params);
    I(load_frozen_fit_params_file);
    I(load_guess);
    I(guess_file);
    I(load_mlparams);
    I(mlparams_file);
#undef I
  }


};


#endif
