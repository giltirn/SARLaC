#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_CMDLINE_H
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_CMDLINE_H

struct CMDline{
  bool load_guess;
  std::string guess_file;

  bool load_combined_data;
  std::string load_combined_data_file;

  bool save_combined_data;
  std::string save_combined_data_file;

  bool load_raw_data;
  std::string load_raw_data_file;

  bool save_raw_data;
  std::string save_raw_data_file;

  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  bool write_covariance_matrix;
  std::string write_covariance_matrix_file;

  bool save_guess_template;

  bool load_priors;
  std::string load_priors_file;

  bool load_mlparams;
  std::string mlparams_file;

  bool remove_samples_in_range;
  int remove_samples_in_range_start; //units are sample index, not trajectories!
  int remove_samples_in_range_lessthan;

  bool write_fit_data;

  CMDline(){
    load_guess = false;
    load_raw_data = false;
    save_raw_data = false;
    load_combined_data = false;
    save_combined_data = false;
    load_frozen_fit_params = false;
    save_guess_template = false;
    write_covariance_matrix = false;
    load_priors = false;
    load_mlparams = false;
    remove_samples_in_range = false;
    write_fit_data = false;
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
      }else if(sargv[i] == "-save_guess_template"){ //save template for guess file then exit
	save_guess_template = true;
	i++;
      }else if(sargv[i] == "-nthread"){
	omp_set_num_threads(strToAny<int>(sargv[i+1]));
	i+=2;
      }else if(sargv[i] == "-load_raw_data"){ //load the raw data pre vacuum subtraction/binning
	load_raw_data = true;
	load_raw_data_file = sargv[i+1];
	std::cout << "Enabled reading raw data from " << load_raw_data_file << std::endl;
	i+=2;
      }else if(sargv[i] == "-save_raw_data"){ //save the raw data pre vacuum subtraction/binning
	save_raw_data = true;
	save_raw_data_file = sargv[i+1];
	std::cout << "Enabled saving raw data to " << save_raw_data_file << std::endl;
	i+=2;
      }else if(sargv[i] == "-load_combined_data"){ //load the double-jackknife data set previously generated
	load_combined_data = true;
	load_combined_data_file = sargv[i+1];
	std::cout << "Enabled reading combined data from " << load_combined_data_file << std::endl;
	i+=2;
      }else if(sargv[i] == "-save_combined_data"){ //save the double-jackknife data set previously generated
	save_combined_data = true;
	save_combined_data_file = sargv[i+1];
	std::cout << "Enabled saving combined data to " << save_combined_data_file << std::endl;
	i+=2;
      }else if(sargv[i] == "-write_covariance_matrix"){
	write_covariance_matrix = true;
	write_covariance_matrix_file = sargv[i+1];
	std::cout << "Enabled saving covariance matrix to " << write_covariance_matrix_file << std::endl;
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
      }else if(sargv[i] == "-load_priors"){
	load_priors = true;
	load_priors_file = sargv[i+1];
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
      }else if(sargv[i] == "-allow_bin_cropping"){ //when #configs is not an exact multiple of bin size, allow discarding of excess configs
	rawDataDistributionOptions::binAllowCropByDefault() = true;
	i++;
      }else if(sargv[i] == "-remove_samples_in_range"){ //drop data from raw data being read. Only works if reading raw data from original files or checkpoint
	remove_samples_in_range = true;
	remove_samples_in_range_start = strToAny<int>(sargv[i+1]);
	remove_samples_in_range_lessthan = strToAny<int>(sargv[i+2]);
	i+=3;
      }else if(sargv[i] == "-write_fit_data"){
	write_fit_data=true;
	i++;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
  
  void exportOptions(fitOptions &opt){
#define COPYIT(A) opt.A = A
    COPYIT(load_frozen_fit_params);
    COPYIT(load_frozen_fit_params_file);
    COPYIT(write_covariance_matrix);
    COPYIT(write_covariance_matrix_file);
    COPYIT(load_priors);
    COPYIT(load_priors_file);
    COPYIT(load_mlparams);
    COPYIT(mlparams_file);
  }

};


#endif
