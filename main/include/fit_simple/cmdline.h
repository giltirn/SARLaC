#ifndef _FIT_SIMPLE_CMDLINE_H
#define _FIT_SIMPLE_CMDLINE_H

struct CMDline{
  bool load_guess;
  std::string guess_file;

  bool save_combined_data;
  std::string save_combined_data_file;

  bool load_combined_data;
  std::string load_combined_data_file;

  bool save_raw_data;
  std::string save_raw_data_file;

  bool load_raw_data;
  std::string load_raw_data_file;

  bool remove_samples_in_range;
  int remove_samples_in_range_start; //units are sample index, not trajectories!
  int remove_samples_in_range_lessthan;

  bool load_mlparams;
  std::string mlparams_file;

  bool scramble_raw_data; //randomly scramble samples of raw data (only works if reading/loading raw data distributions!)

  bool load_boot_resample_table;
  std::string load_boot_resample_table_file;

  bool save_boot_resample_table;
  std::string save_boot_resample_table_file;

  bool save_bootstrap_fit_result;
  std::string save_bootstrap_fit_result_file;

  bool print_fitdata_samples;

  CMDline(){
    load_guess = false;
    save_combined_data = false;
    load_combined_data = false;
    save_raw_data = false;
    load_raw_data = false;  
    remove_samples_in_range = false;
    load_mlparams = false;
    scramble_raw_data = false;
    load_boot_resample_table = false;
    save_boot_resample_table = false;
    save_bootstrap_fit_result = false;
    print_fitdata_samples = false;
  }
  CMDline(const int argc, const char** argv, const int begin = 0): CMDline(){
    setup(argc,argv,begin);
  }
  
  void setup(const int argc, const char** argv, const int begin = 0){
    const int sz = argc-begin;
    if(sz <= 0) return;

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
      }else if(sargv[i] == "-save_combined_data"){ //save the resampled data
	save_combined_data = true;
	save_combined_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_combined_data"){ //load the resampled data
	load_combined_data = true;
	load_combined_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_raw_data"){ //save the raw, unbinned data
	save_raw_data = true;
	save_raw_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_raw_data"){ //load the raw, unbinned data
	load_raw_data = true;
	load_raw_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-allow_bin_cropping"){ //when #configs is not an exact multiple of bin size, allow discarding of excess configs
	rawDataDistributionOptions::binAllowCropByDefault() = true;
	i++;
      }else if(sargv[i] == "-remove_samples_in_range"){ //drop data from raw data being read. Only works if reading raw data from original files or checkpoint
	remove_samples_in_range = true;
	remove_samples_in_range_start = strToAny<int>(sargv[i+1]);
	remove_samples_in_range_lessthan = strToAny<int>(sargv[i+2]);
	i+=3;
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
      }else if(sargv[i] == "-scramble_raw_data"){
	scramble_raw_data = true;
	i++;
      }else if(sargv[i] == "-load_boot_resample_table"){ 
	load_boot_resample_table = true;
	load_boot_resample_table_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_boot_resample_table"){ 
	save_boot_resample_table = true;
	save_boot_resample_table_file = sargv[i+1];
      	i+=2;
      }else if(sargv[i] == "-save_bootstrap_fit_result"){ //When bootstrap p-value is computed, store the fit to the (unrecentered) data
	save_bootstrap_fit_result = true;
	save_bootstrap_fit_result_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-print_fitdata_samples"){ //print the jackknife samples for data in the fit window
	print_fitdata_samples = true;
	i++;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};


#endif
