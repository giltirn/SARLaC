#ifndef _FIT_PIPI_GND_EXC_SIGMA_GEVP_GPARITY_CMDLINE_H_
#define _FIT_PIPI_GND_EXC_SIGMA_GEVP_GPARITY_CMDLINE_H_

struct CMDline{
  bool load_combined_data;
  std::string load_combined_data_file;

  bool save_combined_data;
  std::string save_combined_data_file;

  bool load_raw_data;
  std::string load_raw_data_file;

  bool save_raw_data;
  std::string save_raw_data_file;

  bool verbose_solver;
  
  bool subtract_from_data = false;
  std::string subtract_from_data_file;

  CMDline(){
    load_raw_data = false;
    save_raw_data = false;
    load_combined_data = false;
    save_combined_data = false;
    verbose_solver = false;
    subtract_from_data = false;
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
      }else if(sargv[i] == "-allow_bin_cropping"){ //when #configs is not an exact multiple of bin size, allow discarding of excess configs
	rawDataDistributionOptions::binAllowCropByDefault() = true;
	i++;
      }else if(sargv[i] == "-verbose_solver"){
	verbose_solver = true;
	i++;
      }else if(sargv[i] == "-subtract_from_data"){
	subtract_from_data = true;
	subtract_from_data_file = sargv[i+1];
	i+=2;

	if(subtract_from_data_file == "TEMPLATE"){
	  SubArgs templ;
	  std::ofstream of("subargs_template.args");
	  of << templ;
	  of.close();
	  std::cout << "Wrote SubArgs template to subargs_template.args\n";
	  exit(0);	 
	}

      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};


#endif
