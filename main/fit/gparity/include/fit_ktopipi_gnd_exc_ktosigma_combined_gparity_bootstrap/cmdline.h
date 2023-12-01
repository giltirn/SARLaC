#ifndef _FIT_KTOPIPI_KTOSIGMA_GPARITY_BOOTSTRAP_CMDLINE_H
#define _FIT_KTOPIPI_KTOSIGMA_GPARITY_BOOTSTRAP_CMDLINE_H

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE


struct CMDline{
#define DEF_ARGS(NM) \
  bool load_##NM##_data_checkpoint; \
  std::string load_##NM##_data_checkpoint_stub; \
  bool save_##NM##_data_checkpoint; \
  std::string save_##NM##_data_checkpoint_stub; \
  bool load_##NM##_amplitude_data; \
  std::string load_##NM##_amplitude_data_file; \
  bool save_##NM##_amplitude_data; \
  std::string save_##NM##_amplitude_data_file

  DEF_ARGS(ktopipi);
  DEF_ARGS(ktopipi_exc);
  DEF_ARGS(ktosigma);
#undef DEF_ARGS
 
  bool load_raw_data_container_checkpoint;
  std::string load_raw_data_container_checkpoint_file;

  bool save_raw_data_container_checkpoint;
  std::string save_raw_data_container_checkpoint_file;

  bool load_resampled_data_container_checkpoint;
  std::string load_resampled_data_container_checkpoint_file;

  bool save_resampled_data_container_checkpoint;
  std::string save_resampled_data_container_checkpoint_file;

  bool load_boot_resample_table;
  std::string load_boot_resample_table_file;

  bool save_boot_resample_table;
  std::string save_boot_resample_table_file;

  bool fronthalf_backhalf;

  bool remove_samples_in_range;
  int remove_samples_in_range_start; //units are sample index, not trajectories!
  int remove_samples_in_range_lessthan;

  bool write_alpha_and_pseudoscalar_matrix_elem;

  bool alpha_vary_plot;  
  std::string alpha_vary_plot_args;

  bool set_alpha_coeff; //change the coefficient of alpha from 1.0 to some other value to illustrate effect
  double alpha_coeff;

  bool prune_data;
  std::string prune_data_file;

  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  bool disable_vacuum_subtraction;

  bool exclude_V_diagram; //exclude disconnected diagrams

  CMDline(){
#define INIT_ARGS(NM)				\
    load_##NM##_data_checkpoint = false;	\
    save_##NM##_data_checkpoint = false;	\
    load_##NM##_amplitude_data = false;	\
    save_##NM##_amplitude_data = false

    INIT_ARGS(ktopipi);
    INIT_ARGS(ktopipi_exc);
    INIT_ARGS(ktosigma);
#undef INIT_ARGS

    load_raw_data_container_checkpoint = false;
    save_raw_data_container_checkpoint = false;

    load_resampled_data_container_checkpoint = false;
    save_resampled_data_container_checkpoint = false;

    load_boot_resample_table = false;
    save_boot_resample_table = false;

    fronthalf_backhalf = false;

    alpha_vary_plot = false;

    write_alpha_and_pseudoscalar_matrix_elem = false;

    set_alpha_coeff = false;

    prune_data = false;

    load_frozen_fit_params = false;

    disable_vacuum_subtraction = false;

    remove_samples_in_range = false;

    exclude_V_diagram = false;
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
      if(sargv[i] == "-nthread"){
	omp_set_num_threads(strToAny<int>(sargv[i+1]));
	i+=2;
      }else if(sargv[i] == "-allow_bin_cropping"){ //when #configs is not an exact multiple of bin size, allow discarding of excess configs
	rawDataDistributionOptions::binAllowCropByDefault() = true;
	i++;

#define PARSE_ARGS(NM) \
      }else if(sargv[i] == "-load_" #NM "_data_checkpoint"){ \
	load_##NM##_data_checkpoint = true;		     \
        load_##NM##_data_checkpoint_stub = sargv[i+1];     \
        i+=2;					     \
      }else if(sargv[i] == "-save_" #NM "_data_checkpoint"){ \
	save_##NM##_data_checkpoint = true;		     \
        save_##NM##_data_checkpoint_stub = sargv[i+1];     \
        i+=2;					     \
      }else if(sargv[i] == "-load_" #NM "_amplitude_data"){\
	load_##NM##_amplitude_data = true;\
	load_##NM##_amplitude_data_file = sargv[i+1];\
	i+=2;\
      }else if(sargv[i] == "-save_" #NM "_amplitude_data"){\
	save_##NM##_amplitude_data = true;\
	save_##NM##_amplitude_data_file = sargv[i+1];\
	i+=2

      PARSE_ARGS(ktopipi);
      PARSE_ARGS(ktopipi_exc);
      PARSE_ARGS(ktosigma);
#undef PARSE_ARGS
      }else if(sargv[i] == "-save_raw_data_container_checkpoint"){
	save_raw_data_container_checkpoint = true;
	save_raw_data_container_checkpoint_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_raw_data_container_checkpoint"){
	load_raw_data_container_checkpoint = true;
	load_raw_data_container_checkpoint_file = sargv[i+1];
	i+=2;      
      }else if(sargv[i] == "-save_resampled_data_container_checkpoint"){
	save_resampled_data_container_checkpoint = true;
	save_resampled_data_container_checkpoint_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_resampled_data_container_checkpoint"){
	load_resampled_data_container_checkpoint = true;
	load_resampled_data_container_checkpoint_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_boot_resample_table"){ 
	load_boot_resample_table = true;
	load_boot_resample_table_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_boot_resample_table"){ 
	save_boot_resample_table = true;
	save_boot_resample_table_file = sargv[i+1];
      	i+=2;
      }else if(sargv[i] == "-fronthalf_backhalf"){
	fronthalf_backhalf = true;
	i++;
      }else if(sargv[i] == "-remove_samples_in_range"){ //drop data from raw data being read. Only works if reading raw data from original files or checkpoint
	remove_samples_in_range = true;
	remove_samples_in_range_start = strToAny<int>(sargv[i+1]);
	remove_samples_in_range_lessthan = strToAny<int>(sargv[i+2]);
	i+=3;
      }else if(sargv[i] == "-write_alpha_and_pseudoscalar_matrix_elem"){
	write_alpha_and_pseudoscalar_matrix_elem = true;
	i++;
      }else if(sargv[i] == "-alpha_vary_plot"){
	alpha_vary_plot = true;
	alpha_vary_plot_args = sargv[i+1];
	i+=2;
	std::cout << "Doing alpha vary plot with args from " << alpha_vary_plot_args << std::endl;

	if(alpha_vary_plot_args == "TEMPLATE"){
	  AlphaVaryPlotArgs pargs;
	  std::ofstream of("alpha_vary_plot_arg.args"); of << pargs; of.close();	 
	  exit(0);
	}
      }else if(sargv[i] == "-set_alpha_coeff"){
	set_alpha_coeff = true;
	alpha_coeff = strToAny<double>(sargv[i+1]);
	i+=2;
      }else if(sargv[i] == "-prune_data"){
	prune_data = true;
	prune_data_file = sargv[i+1];
	i+=2;

	if(prune_data_file == "TEMPLATE"){
	  pruneArgs pargs;
	  std::ofstream of("prune_args_template.args"); of << pargs; of.close();
	  exit(0);
	}
      }else if(sargv[i] == "-disable_vacuum_subtraction"){
	disable_vacuum_subtraction = true;
	std::cout << "Not performing vacuum subtraction" << std::endl;
	i++;
      }else if(sargv[i] == "-load_frozen_fit_params"){
	load_frozen_fit_params = true;
	load_frozen_fit_params_file = sargv[i+1];
	i+=2;

	if(load_frozen_fit_params_file == "TEMPLATE"){
	  FreezeMatrixElements pargs;
	  std::ofstream of("freeze_template.args"); of << pargs; of.close();
	  exit(0);
	}
      }else if(sargv[i] == "-exclude_V_diagram"){
	exclude_V_diagram = true;
	std::cout << "Excluding disconnected diagrams" << std::endl;
	i++;      
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};



CPSFIT_END_NAMESPACE

#endif
