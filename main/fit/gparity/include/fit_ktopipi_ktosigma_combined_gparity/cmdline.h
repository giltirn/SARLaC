#ifndef _FIT_KTOPIPI_KTOSIGMA_GPARITY_CMDLINE_H
#define _FIT_KTOPIPI_KTOSIGMA_GPARITY_CMDLINE_H

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE


struct CMDline{
  bool load_ktopipi_data_checkpoint;
  std::string load_ktopipi_data_checkpoint_stub; //will append  _tsepkpi<VAL>.hdf5
  bool save_ktopipi_data_checkpoint;
  std::string save_ktopipi_data_checkpoint_stub;

  bool load_ktosigma_data_checkpoint;
  std::string load_ktosigma_data_checkpoint_stub;
  bool save_ktosigma_data_checkpoint;
  std::string save_ktosigma_data_checkpoint_stub;

  bool load_ktopipi_amplitude_data;
  std::string load_ktopipi_amplitude_data_file;
  bool save_ktopipi_amplitude_data;
  std::string save_ktopipi_amplitude_data_file;

  bool load_ktosigma_amplitude_data;
  std::string load_ktosigma_amplitude_data_file;
  bool save_ktosigma_amplitude_data;
  std::string save_ktosigma_amplitude_data_file;


  bool tianle_compare;
  std::string tianle_compare_file;

  bool skip_frzparam_load; //for testing, just read the data and exit without loading the frozen parameters from the 2pt fits
  
  bool bind_M1; //for simultaneous fit, choose whether the M1 matrix element is shared "bound" between the two operators or allowed to differ

  CMDline(){
    load_ktopipi_data_checkpoint = false;
    save_ktopipi_data_checkpoint = false;
    load_ktopipi_amplitude_data = false;
    save_ktopipi_amplitude_data = false;

    load_ktosigma_data_checkpoint = false;
    save_ktosigma_data_checkpoint = false;
    load_ktosigma_amplitude_data = false;
    save_ktosigma_amplitude_data = false;
    tianle_compare = false;
    skip_frzparam_load = false;
    bind_M1 = true;
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


      }else if(sargv[i] == "-load_ktopipi_data_checkpoint"){
	load_ktopipi_data_checkpoint = true;
	load_ktopipi_data_checkpoint_stub = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_ktopipi_data_checkpoint"){
	save_ktopipi_data_checkpoint = true;
	save_ktopipi_data_checkpoint_stub = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_ktopipi_amplitude_data"){
	load_ktopipi_amplitude_data = true;
	load_ktopipi_amplitude_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_ktopipi_amplitude_data"){
	save_ktopipi_amplitude_data = true;
	save_ktopipi_amplitude_data_file = sargv[i+1];
	i+=2;


      }else if(sargv[i] == "-load_ktosigma_data_checkpoint"){
	load_ktosigma_data_checkpoint = true;
	load_ktosigma_data_checkpoint_stub = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_ktosigma_data_checkpoint"){
	save_ktosigma_data_checkpoint = true;
	save_ktosigma_data_checkpoint_stub = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_ktosigma_amplitude_data"){
	load_ktosigma_amplitude_data = true;
	load_ktosigma_amplitude_data_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_ktosigma_amplitude_data"){
	save_ktosigma_amplitude_data = true;
	save_ktosigma_amplitude_data_file = sargv[i+1];
	i+=2;

	
      }else if(sargv[i] == "-tianle_compare"){
	tianle_compare = true;
	tianle_compare_file = sargv[i+1];
	i+=2;

      }else if(sargv[i] == "-skip_frzparam_load"){
	skip_frzparam_load = true;
	i++;

      }else if(sargv[i] == "-unbind_M1"){
	bind_M1 = false;
	i++;

      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};


SARLAC_END_NAMESPACE

#endif
