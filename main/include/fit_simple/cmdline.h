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


  CMDline(){
    load_guess = false;
    save_combined_data = false;
    load_combined_data = false;
    save_raw_data = false;
    load_raw_data = false;  
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
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};


#endif
