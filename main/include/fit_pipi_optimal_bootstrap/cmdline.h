#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_CMDLINE_BOOTSTRAP_H
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_CMDLINE_BOOTSTRAP_H

struct CMDline{
  bool load_guess;
  std::string guess_file;

  bool load_raw_data;
  std::string load_raw_data_file;

  bool save_raw_data;
  std::string save_raw_data_file;

  bool load_combined_data;
  std::string load_combined_data_file;

  bool save_combined_data;
  std::string save_combined_data_file;

  bool save_guess_template;

  bool load_minimizer_params;
  std::string minimizer_params_file;

  bool load_boot_resample_table;
  std::string load_boot_resample_table_file;

  bool save_boot_resample_table;
  std::string save_boot_resample_table_file;

  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  CMDline(){
    load_guess = false;
    load_raw_data = false;
    save_raw_data = false;
    load_combined_data = false;
    save_combined_data = false;
    save_guess_template = false;
    load_minimizer_params = false;
    load_boot_resample_table = false;
    save_boot_resample_table = false;
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
      }else if(sargv[i] == "-load_minimizer_params"){
	load_minimizer_params = true;
	minimizer_params_file = sargv[i+1];
	if(minimizer_params_file == "TEMPLATE"){
	  std::ofstream of("min_params_template.args");
	  MarquardtLevenbergParameters<double> templ;
	  of << templ;
	  of.close();
	  std::cout << "Wrote minimizer params template to min_params_template.args\n";

	  exit(0);	 
	}
	i+=2;
      }else if(sargv[i] == "-allow_bin_cropping"){ //when #configs is not an exact multiple of bin size, allow discarding of excess configs
	rawDataDistributionOptions::binAllowCropByDefault() = true;
	i++;
      }else if(sargv[i] == "-load_boot_resample_table"){ 
	load_boot_resample_table = true;
	load_boot_resample_table_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_boot_resample_table"){ 
	save_boot_resample_table = true;
	save_boot_resample_table_file = sargv[i+1];
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
