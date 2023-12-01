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
  
  bool load_resampled_data;
  std::string load_resampled_data_file;

  bool save_resampled_data;
  std::string save_resampled_data_file;
  
  bool use_pipitosigma_disconn_complex_prod; //use Re( pipi_bubble * sigma_bubble )  [original strategy] for pipi->sigma disconnected piece rather than Re ( pipi_bubble ) * Re ( sigma_bubble )

  bool write_covariance_matrix;
  std::string write_covariance_matrix_file;

  bool analyze_population_distribution; //analyze the population distribution using the bootstrap method to ensure gaussianity. Requires either reading original
                                        //data files or loading the raw data checkpoint

  bool ledoit_wolf; //apply ledoit wolf procedure to correlation matrix

  bool load_mlparams;
  std::string mlparams_file;

  PiPiSigmaSimCMDline(){
    load_guess = false;
    load_frozen_fit_params= false;
    include_pipi_2pt = true;
    include_pipi_to_sigma = true;
    include_sigma_2pt = true;
    save_checkpoint = false;
    load_checkpoint = false;
    load_resampled_data = false;
    save_resampled_data = false;
    use_pipitosigma_disconn_complex_prod = false;
    write_covariance_matrix = false;
    analyze_population_distribution = false;
    ledoit_wolf = false;
    load_mlparams = false;
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
	}else if(!fileExists(load_frozen_fit_params_file)) error_exit(std::cout << "PiPiSigmaSimCMDline freeze data file " << load_frozen_fit_params_file << " does not exist!\n");
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
      }else if(sargv[i] == "-save_resampled_data"){
	save_resampled_data = true;
	save_resampled_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_resampled_data"){
	load_resampled_data = true;
	load_resampled_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-use_pipitosigma_disconn_complex_prod"){
	use_pipitosigma_disconn_complex_prod = true;
	i++;
      }else if(sargv[i] == "-write_covariance_matrix"){
	write_covariance_matrix = true;
	write_covariance_matrix_file = sargv[i+1];
	std::cout << "Enabled saving covariance matrix to " << write_covariance_matrix_file << std::endl;
	i+=2;
      }else if(sargv[i] == "-allow_bin_cropping"){ //when binning allow extra configs that don't fit a full bin to be cropped from the end of the ensemble
	rawDataDistributionOptions::binAllowCropByDefault() = true;
	i++;
      }else if(sargv[i] == "-analyze_population_distribution"){
	analyze_population_distribution = true;
	i++;
      }else if(sargv[i] == "-ledoit_wolf"){
	ledoit_wolf = true;
	i++;
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
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }

  void transfer(SimFitArgs &fargs) const{
#define COPYIT(A) fargs.A = A
    COPYIT(load_guess);
    COPYIT(guess_file);
    COPYIT(load_frozen_fit_params);
    COPYIT(load_frozen_fit_params_file);
    COPYIT(write_covariance_matrix);
    COPYIT(write_covariance_matrix_file);
    COPYIT(ledoit_wolf);
    COPYIT(load_mlparams);
    COPYIT(mlparams_file);
#undef COPYIT
  }
};


#endif
