#ifndef _PIPI_SIGMA_GEVP_FIT_CMDLINE_H_
#define _PIPI_SIGMA_GEVP_FIT_CMDLINE_H_

struct PiPiSigmaGEVPcmdLine{
  bool save_checkpoint;
  std::string save_checkpoint_file;
  
  bool load_checkpoint;
  std::string load_checkpoint_file;
  
  bool load_resampled_data;
  std::string load_resampled_data_file;

  bool save_resampled_data;
  std::string save_resampled_data_file;
  
  bool use_pipitosigma_disconn_complex_prod; //use Re( pipi_bubble * sigma_bubble )  [original strategy] for pipi->sigma disconnected piece rather than Re ( pipi_bubble ) * Re ( sigma_bubble )

  PiPiSigmaGEVPcmdLine(){
    save_checkpoint = false;
    load_checkpoint = false;
    load_resampled_data = false;
    save_resampled_data = false;
    use_pipitosigma_disconn_complex_prod = false;
  }
  PiPiSigmaGEVPcmdLine(const int argc, const char** argv, const int begin = 0): PiPiSigmaGEVPcmdLine(){
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
      }else if(sargv[i] == "-allow_bin_cropping"){ //when binning allow extra configs that don't fit a full bin to be cropped from the end of the ensemble
	rawDataDistributionOptions::binAllowCropByDefault() = true;
	i++;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};


#endif
