#ifndef _FIT_PIPI_GPARITY_SAMPLEAMA_CMDLINE_H
#define _FIT_PIPI_GPARITY_SAMPLEAMA_CMDLINE_H

enum SloppyExact {Sloppy, Exact};

struct CMDlineSampleAMA{
  bool load_guess;
  std::string guess_file;

  bool load_hdf5_data_checkpoint_sloppy_S;
  std::string load_hdf5_data_checkpoint_stub_sloppy_S;
  bool load_hdf5_data_checkpoint_sloppy_C;
  std::string load_hdf5_data_checkpoint_stub_sloppy_C;
  bool load_hdf5_data_checkpoint_exact_C;
  std::string load_hdf5_data_checkpoint_stub_exact_C;

  bool load_checkpoint(std::string &filename, const SloppyExact se, const char ens) const{
    if(se == Sloppy && ens == 'S'){ filename = load_hdf5_data_checkpoint_stub_sloppy_S; return load_hdf5_data_checkpoint_sloppy_S; }
    if(se == Sloppy && ens == 'C'){ filename = load_hdf5_data_checkpoint_stub_sloppy_C; return load_hdf5_data_checkpoint_sloppy_C; }
    if(se == Exact && ens == 'C'){ filename = load_hdf5_data_checkpoint_stub_exact_C; return load_hdf5_data_checkpoint_exact_C; }
    assert(0);
  }

  bool save_hdf5_data_checkpoint_sloppy_S;
  std::string save_hdf5_data_checkpoint_stub_sloppy_S;
  bool save_hdf5_data_checkpoint_sloppy_C;
  std::string save_hdf5_data_checkpoint_stub_sloppy_C;
  bool save_hdf5_data_checkpoint_exact_C;
  std::string save_hdf5_data_checkpoint_stub_exact_C;

  bool save_checkpoint(std::string &filename, const SloppyExact se, const char ens) const{
    if(se == Sloppy && ens == 'S'){ filename = save_hdf5_data_checkpoint_stub_sloppy_S; return save_hdf5_data_checkpoint_sloppy_S; }
    if(se == Sloppy && ens == 'C'){ filename = save_hdf5_data_checkpoint_stub_sloppy_C; return save_hdf5_data_checkpoint_sloppy_C; }
    if(se == Exact && ens == 'C'){ filename = save_hdf5_data_checkpoint_stub_exact_C; return save_hdf5_data_checkpoint_exact_C; }
    assert(0);
  }


  bool load_combined_data;
  std::string load_combined_data_file;

  bool save_combined_data;
  std::string save_combined_data_file;

  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  
  CMDlineSampleAMA(){
    load_guess = false;
    load_hdf5_data_checkpoint_sloppy_S = false;
    load_hdf5_data_checkpoint_sloppy_C = false;
    load_hdf5_data_checkpoint_exact_C = false;
    save_hdf5_data_checkpoint_sloppy_S = false;
    save_hdf5_data_checkpoint_sloppy_C = false;
    save_hdf5_data_checkpoint_exact_C = false;
    load_combined_data = false;
    save_combined_data = false;
    load_frozen_fit_params= false;
  }
  CMDlineSampleAMA(const int argc, const char** argv, const int begin = 0): CMDlineSampleAMA(){
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
      }else if(sargv[i] == "-load_hdf5_data_checkpoint_sloppy_S"){
	load_hdf5_data_checkpoint_sloppy_S = true;
	load_hdf5_data_checkpoint_stub_sloppy_S = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_hdf5_data_checkpoint_sloppy_C"){
	load_hdf5_data_checkpoint_sloppy_C = true;
	load_hdf5_data_checkpoint_stub_sloppy_C = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_hdf5_data_checkpoint_exact_C"){
	load_hdf5_data_checkpoint_exact_C = true;
	load_hdf5_data_checkpoint_stub_exact_C = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_hdf5_data_checkpoint_sloppy_S"){
	save_hdf5_data_checkpoint_sloppy_S = true;
	save_hdf5_data_checkpoint_stub_sloppy_S = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_hdf5_data_checkpoint_sloppy_C"){
	save_hdf5_data_checkpoint_sloppy_C = true;
	save_hdf5_data_checkpoint_stub_sloppy_C = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_hdf5_data_checkpoint_exact_C"){
	save_hdf5_data_checkpoint_exact_C = true;
	save_hdf5_data_checkpoint_stub_exact_C = sargv[i+1];
	i+=2;
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
