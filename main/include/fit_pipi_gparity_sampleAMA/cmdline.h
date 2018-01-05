#ifndef _FIT_PIPI_GPARITY_SAMPLEAMA_CMDLINE_H
#define _FIT_PIPI_GPARITY_SAMPLEAMA_CMDLINE_H

enum SloppyExact {Sloppy, Exact};

struct CMDlineSampleAMA{
  bool load_guess;
  std::string guess_file;

  bool load_hdf5_data_checkpoint;
  std::string load_hdf5_data_checkpoint_file_sloppy_S;
  std::string load_hdf5_data_checkpoint_file_sloppy_C;
  std::string load_hdf5_data_checkpoint_file_exact_C;

  bool save_hdf5_data_checkpoint;
  std::string save_hdf5_data_checkpoint_file_sloppy_S;
  std::string save_hdf5_data_checkpoint_file_sloppy_C;
  std::string save_hdf5_data_checkpoint_file_exact_C;

  bool load_combined_data;
  std::string load_combined_data_file;

  bool save_combined_data;
  std::string save_combined_data_file;

  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  
  CMDlineSampleAMA(){
    load_guess = false;
    load_hdf5_data_checkpoint = false;
    save_hdf5_data_checkpoint = false;
    load_combined_data = false;
    save_combined_data = false;
    load_frozen_fit_params= false;
  }
  CMDlineSampleAMA(const int argc, const char** argv, const int begin = 0): CMDlineSampleAMA(){
    setup(argc,argv,begin);
  }
  
  CMDline toCMDline(const SloppyExact &se, const char ens) const{
    CMDline out;
    out.load_guess = load_guess;
    out.guess_file = guess_file;
    out.load_combined_data = load_combined_data;
    out.load_combined_data_file = load_combined_data_file;
    out.save_combined_data = save_combined_data;
    out.save_combined_data_file = save_combined_data_file;
    out.load_frozen_fit_params = load_frozen_fit_params;
    out.load_frozen_fit_params_file = load_frozen_fit_params_file;

    out.load_hdf5_data_checkpoint = load_hdf5_data_checkpoint;
    out.save_hdf5_data_checkpoint = save_hdf5_data_checkpoint;
    if(se == Sloppy){
      out.use_symmetric_quark_momenta = false;
      out.load_hdf5_data_checkpoint_file = ens == 'S' ? load_hdf5_data_checkpoint_file_sloppy_S : load_hdf5_data_checkpoint_file_sloppy_C;
      out.save_hdf5_data_checkpoint_file = ens == 'S' ? save_hdf5_data_checkpoint_file_sloppy_S : save_hdf5_data_checkpoint_file_sloppy_C;
    }else{
      assert(ens == 'C');
      out.use_symmetric_quark_momenta = true;
      out.load_hdf5_data_checkpoint_file = load_hdf5_data_checkpoint_file_exact_C;
      out.save_hdf5_data_checkpoint_file = save_hdf5_data_checkpoint_file_exact_C;
    }
    return out;
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
      }else if(sargv[i] == "-load_hdf5_data_checkpoint"){
	load_hdf5_data_checkpoint = true;
	load_hdf5_data_checkpoint_file_sloppy_S = sargv[i+1];
	load_hdf5_data_checkpoint_file_sloppy_C = sargv[i+2];
	load_hdf5_data_checkpoint_file_exact_C = sargv[i+3];
	i+=4;
      }else if(sargv[i] == "-save_hdf5_data_checkpoint"){
	save_hdf5_data_checkpoint = true;
	save_hdf5_data_checkpoint_file_sloppy_S = sargv[i+1];
	save_hdf5_data_checkpoint_file_sloppy_C = sargv[i+2];
	save_hdf5_data_checkpoint_file_exact_C = sargv[i+3];
	i+=4;
      }else if(sargv[i] == "-load_combined_data"){ //load the double-jackknife data set previously generated
	load_combined_data = true;
	load_combined_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_combined_data"){ //save the double-jackknife data set previously generated
	save_combined_data = true;
	save_combined_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_frozen_fit_params"){
	load_frozen_fit_params = true;
	load_frozen_fit_params_file = sargv[i+1];
	i+=2;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};


#endif
