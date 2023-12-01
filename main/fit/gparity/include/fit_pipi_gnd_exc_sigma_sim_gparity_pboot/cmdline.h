#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_PBOOT_CMDLINE_H
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_PBOOT_CMDLINE_H

struct CMDline{
  bool load_guess;
  std::string guess_file;

  bool load_raw_data;
  std::string load_raw_data_file;

  bool save_raw_data;
  std::string save_raw_data_file;

  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  bool save_guess_template;

  bool load_priors;
  std::string load_priors_file;

  bool load_bounds;
  std::string load_bounds_file;

  bool load_minimizer_params;
  std::string minimizer_params_file;

  bool load_filters;
  std::string load_filters_file;

  bool load_boot_resample_table;
  std::string load_boot_resample_table_file;

  bool save_boot_resample_table;
  std::string save_boot_resample_table_file;

  CMDline(){
    load_guess = false;
    load_raw_data = false;
    save_raw_data = false;
    load_frozen_fit_params = false;
    save_guess_template = false;
    load_priors = false;
    load_bounds = false;
    load_minimizer_params = false;
    load_filters = false;
    load_boot_resample_table = false;
    save_boot_resample_table = false;
  }
  CMDline(const int argc, const char** argv, MinimizerType minimizer, const int begin = 0): CMDline(){
    setup(argc,argv,minimizer,begin);
  }
  
  void setup(const int argc, const char** argv, MinimizerType minimizer, const int begin = 0){
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

      }else if(sargv[i] == "-load_filters"){
	load_filters = true;
	load_filters_file = sargv[i+1];

	if(load_filters_file == "TEMPLATE"){
	  std::cout << "Saving filters template file to filters_template.args" << std::endl;
	  Filters fp;
	  std::ofstream of("filters_template.args");
	  of << fp;
	  of.close();
	  exit(0);
	}else if(!fileExists(load_filters_file)) error_exit(std::cout << "CMDline filters file " << load_filters_file << " does not exist!\n");
	i+=2;

      }else if(sargv[i] == "-load_priors"){
	load_priors = true;
	load_priors_file = sargv[i+1];
	i+=2;

      }else if(sargv[i] == "-load_bounds"){
	load_bounds = true;
	load_bounds_file = sargv[i+1];

	if(load_bounds_file == "TEMPLATE"){
	  std::cout << "Saving bounds template file to bounds_template.args" << std::endl;
	  BoundArgs fp;
	  std::ofstream of("bounds_template.args");
	  of << fp;
	  of.close();
	  exit(0);
	}else if(!fileExists(load_bounds_file)) error_exit(std::cout << "CMDline bounds file " << load_bounds_file << " does not exist!\n");
	i+=2;

      }else if(sargv[i] == "-load_minimizer_params"){
	load_minimizer_params = true;
	minimizer_params_file = sargv[i+1];
	if(minimizer_params_file == "TEMPLATE"){
	  std::ofstream of("min_params_template.args");
	  if(minimizer == MinimizerType::MarquardtLevenberg){
	    MarquardtLevenbergParameters<double> templ;
	    of << templ;
	  }else if(minimizer == MinimizerType::GSLtrs){
	    GSLtrsMinimizerParams templ;
	    of << templ;
	  }else if(minimizer == MinimizerType::GSLmultimin){
	    GSLmultidimMinimizerParams templ;
	    of << templ;
	  }else if(minimizer == MinimizerType::Minuit2){
#ifdef HAVE_MINUIT2
	    Minuit2minimizerParams templ;
	    of << templ;
#else
	    error_exit(std::cout << "Library not compiled with Minuit2\n");
#endif
	  }else assert(0);
	  
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
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
  
  void exportOptions(fitCentralOptions &opt) const{
#define COPYIT(A) opt.A = A
    COPYIT(load_frozen_fit_params);
    COPYIT(load_frozen_fit_params_file);
    COPYIT(load_priors);
    COPYIT(load_priors_file);
    COPYIT(load_bounds);
    COPYIT(load_bounds_file);
    COPYIT(load_minimizer_params);
    COPYIT(minimizer_params_file);
#undef COPYIT
  }

};


#endif
