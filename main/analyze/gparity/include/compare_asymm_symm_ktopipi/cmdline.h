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
  
  std::string symmetric_quark_momenta_figure_file_extension;

  ComparisonCMDline(){
    load_data_checkpoint = false;
    save_data_checkpoint = false;
    load_amplitude_data = false;
    save_amplitude_data = false;
    symmetric_quark_momenta_figure_file_extension = "_symm";
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
      }else if(sargv[i] == "-symmetric_quark_momenta_figure_file_extension"){
	symmetric_quark_momenta_figure_file_extension = sargv[i+1];
	i+=2;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }

  readKtoPiPiDataOptions getReadOptions(const AsymmSymm type) const{
    readKtoPiPiDataOptions opt;
    opt.load_data_checkpoint = load_data_checkpoint;
    opt.save_data_checkpoint = save_data_checkpoint;

    if(type == Asymmetric){
      opt.load_data_checkpoint_stub = load_data_checkpoint_asymm_stub;
      opt.save_data_checkpoint_stub = save_data_checkpoint_asymm_stub;
    }else{
      opt.load_data_checkpoint_stub = load_data_checkpoint_symm_stub;
      opt.save_data_checkpoint_stub = save_data_checkpoint_symm_stub;
    }
    return opt;
  }
};


#endif
