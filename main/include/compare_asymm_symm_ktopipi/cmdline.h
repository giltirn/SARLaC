#ifndef _COMPARE_ASYMM_SYMM_KTOPIPI_GPARITY_CMDLINE_H
#define _COMPARE_ASYMM_SYMM_KTOPIPI_GPARITY_CMDLINE_H

enum AsymmSymm { Asymmetric, Symmetric };

struct ComparisonCMDline{
  bool load_guess;
  std::string guess_file;

  bool load_data_checkpoint;
  std::string load_data_checkpoint_asymm_stub; //will append  _tsepkpi<VAL>.hdf5
  std::string load_data_checkpoint_symm_stub; //will append  _tsepkpi<VAL>.hdf5
  
  bool save_data_checkpoint;
  std::string save_data_checkpoint_asymm_stub;
  std::string save_data_checkpoint_symm_stub;

  bool load_amplitude_data;
  std::string load_amplitude_data_asymm_file;
  std::string load_amplitude_data_symm_file;
  
  bool save_amplitude_data;
  std::string save_amplitude_data_asymm_file;
  std::string save_amplitude_data_symm_file;
  
  ComparisonCMDline(){
    load_data_checkpoint = false;
    save_data_checkpoint = false;
    load_amplitude_data = false;
    save_amplitude_data = false;
  }
  ComparisonCMDline(const int argc, const char** argv, const int begin = 0): ComparisonCMDline(){
    setup(argc,argv,begin);
  }
  
  void setup(const int argc, const char** argv, const int begin = 0){
    const int sz = argc-begin;
    std::vector<std::string> sargv(sz);
    for(int i=begin; i<argc; i++) sargv[i-begin] = std::string(argv[i]);

    int i = 0;
    while(i<sz){
      if(sargv[i] == "-load_data_checkpoint"){
	load_data_checkpoint = true;
	load_data_checkpoint_asymm_stub = sargv[i+1];
	load_data_checkpoint_symm_stub = sargv[i+2];
	i+=3;
      }else if(sargv[i] == "-save_data_checkpoint"){
	save_data_checkpoint = true;
	save_data_checkpoint_asymm_stub = sargv[i+1];
	save_data_checkpoint_symm_stub = sargv[i+2];
	i+=3;
      }else if(sargv[i] == "-load_amplitude_data"){
	load_amplitude_data = true;
	load_amplitude_data_asymm_file = sargv[i+1];
	load_amplitude_data_symm_file = sargv[i+2];
	i+=3;
      }else if(sargv[i] == "-save_amplitude_data"){
	save_amplitude_data = true;
	save_amplitude_data_asymm_file = sargv[i+1];
	save_amplitude_data_symm_file = sargv[i+2];
	i+=3;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }

  CMDline toCMDline(const AsymmSymm type) const{
    CMDline out;
    if(load_data_checkpoint){
      out.load_data_checkpoint = true;
      out.load_data_checkpoint_stub = type == Asymmetric ? load_data_checkpoint_asymm_stub : load_data_checkpoint_symm_stub;
    }
    if(save_data_checkpoint){
      out.save_data_checkpoint = true;
      out.save_data_checkpoint_stub = type == Asymmetric ? save_data_checkpoint_asymm_stub : save_data_checkpoint_symm_stub;
    }
    if(load_amplitude_data){
      out.load_amplitude_data = true;
      out.load_amplitude_data_file = type == Asymmetric ? load_amplitude_data_asymm_file : load_amplitude_data_symm_file;
    }
    if(save_amplitude_data){
      out.save_amplitude_data = true;
      out.save_amplitude_data_file = type == Asymmetric ? save_amplitude_data_asymm_file : save_amplitude_data_symm_file;
    }
    out.use_symmetric_quark_momenta = type == Symmetric;
    return out;
  }
  
};


#endif
