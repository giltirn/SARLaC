#ifndef _COMPARE_ASYMM_SYMM_PIPI_GPARITY_CMDLINE_H
#define _COMPARE_ASYMM_SYMM_PIPI_GPARITY_CMDLINE_H

enum AsymmSymm { Asymmetric, Symmetric };

struct ComparisonCMDline{
  bool load_hdf5_data_checkpoint;
  std::string load_hdf5_data_checkpoint_asymm_stub;
  std::string load_hdf5_data_checkpoint_symm_stub;
  
  bool save_hdf5_data_checkpoint;
  std::string save_hdf5_data_checkpoint_asymm_stub;
  std::string save_hdf5_data_checkpoint_symm_stub;

  bool load_combined_data;
  std::string load_combined_data_asymm_file;
  std::string load_combined_data_symm_file;

  bool save_combined_data;
  std::string save_combined_data_asymm_file;
  std::string save_combined_data_symm_file;
  
  ComparisonCMDline(){
    load_hdf5_data_checkpoint = false;
    save_hdf5_data_checkpoint = false;
    load_combined_data = false;
    save_combined_data = false;
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
      if(sargv[i] == "-load_hdf5_data_checkpoint"){
	load_hdf5_data_checkpoint = true;
	load_hdf5_data_checkpoint_asymm_stub = sargv[i+1];
	load_hdf5_data_checkpoint_symm_stub = sargv[i+2];
	i+=3;
      }else if(sargv[i] == "-save_hdf5_data_checkpoint"){
	save_hdf5_data_checkpoint = true;
	save_hdf5_data_checkpoint_asymm_stub = sargv[i+1];
	save_hdf5_data_checkpoint_symm_stub = sargv[i+2];
	i+=3;
      }else if(sargv[i] == "-load_combined_data"){ //load the double-jackknife data set previously generated
	load_combined_data = true;
	load_combined_data_asymm_file = sargv[i+1];
	load_combined_data_symm_file = sargv[i+2];
	i+=3;
      }else if(sargv[i] == "-save_combined_data"){ //save the double-jackknife data set previously generated
	save_combined_data = true;
	save_combined_data_asymm_file = sargv[i+1];
	save_combined_data_symm_file = sargv[i+2];
	i+=3;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }

  CMDline toCMDline(const AsymmSymm type) const{
    CMDline out;
    if(load_hdf5_data_checkpoint){
      out.load_hdf5_data_checkpoint = true;
      out.load_hdf5_data_checkpoint_stub = type == Asymmetric ? load_hdf5_data_checkpoint_asymm_stub : load_hdf5_data_checkpoint_symm_stub;
    }
    if(save_hdf5_data_checkpoint){
      out.save_hdf5_data_checkpoint = true;
      out.save_hdf5_data_checkpoint_stub = type == Asymmetric ? save_hdf5_data_checkpoint_asymm_stub : save_hdf5_data_checkpoint_symm_stub;
    }
    if(load_combined_data){
      out.load_combined_data = true;
      out.load_combined_data_file = type == Asymmetric ? load_combined_data_asymm_file : load_combined_data_symm_file;
    }
    if(save_combined_data){
      out.save_combined_data = true;
      out.save_combined_data_file = type == Asymmetric ? save_combined_data_asymm_file : save_combined_data_symm_file;
    }
    out.use_symmetric_quark_momenta = type == Asymmetric ? false : true;
    return out;    
  }
  
};


#endif
