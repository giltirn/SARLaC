#ifndef _FIT_KTOSIGMA_GPARITY_CMDLINE_H
#define _FIT_KTOSIGMA_GPARITY_CMDLINE_H

#include<config.h>
#include<utils/macros.h>

#include "args.h"

SARLAC_START_NAMESPACE

struct CMDline{
  bool load_guess;
  std::string guess_file;

  bool load_data_checkpoint;
  std::string load_data_checkpoint_stub; //will append  _tsepkpi<VAL>.hdf5
  bool save_data_checkpoint;
  std::string save_data_checkpoint_stub;

  bool load_amplitude_data;
  std::string load_amplitude_data_file;
  bool save_amplitude_data;
  std::string save_amplitude_data_file;
  
  bool load_freeze_data;
  std::string freeze_data;

  bool use_scratch; //save memory at peak times by storing to disk and reloading later
  std::string use_scratch_stub;

  bool use_existing_scratch_files; //if scratch files from a previous run exist, use these

  CMDline(){
    load_guess = false;
    load_data_checkpoint = false;
    save_data_checkpoint = false;
    load_amplitude_data = false;
    save_amplitude_data = false;
    load_freeze_data = false;
    use_scratch = false;
    use_existing_scratch_files = false;
    use_scratch_stub = "scratch";
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
      if(sargv[i] == "-load_guess"){
	load_guess = true;
	guess_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-nthread"){
	omp_set_num_threads(strToAny<int>(sargv[i+1]));
	i+=2;
      }else if(sargv[i] == "-load_data_checkpoint"){
	load_data_checkpoint = true;
	load_data_checkpoint_stub = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_data_checkpoint"){
	save_data_checkpoint = true;
	save_data_checkpoint_stub = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_amplitude_data"){
	load_amplitude_data = true;
	load_amplitude_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_amplitude_data"){
	save_amplitude_data = true;
	save_amplitude_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-use_scratch"){
	use_scratch = true;
	i++;
      }else if(sargv[i] == "-set_scratch_stub"){
	use_scratch_stub = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-use_existing_scratch_files"){
	use_existing_scratch_files = true;
	i++;
      }else if(sargv[i] == "-freeze"){
	load_freeze_data = true;
	freeze_data = sargv[i+1];
	if(freeze_data == "TEMPLATE"){
	  KtoPiPiFreezeParams p;
	  std::ofstream of("freeze_template.args");
	  of << p;
	  of.close();
	  std::cout << "Wrote freeze template argument file to freeze_template.args\n";
	  exit(0);
	}else{
	  std::cout << "Testing freeze params file " << freeze_data << " for correctness" << std::endl;
	  KtoPiPiFreezeParams fparams;
	  parse(fparams,freeze_data);
	}
	i+=2;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};


SARLAC_END_NAMESPACE

#endif
