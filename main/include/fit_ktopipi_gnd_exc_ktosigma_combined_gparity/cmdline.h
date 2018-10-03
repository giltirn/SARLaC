#ifndef _FIT_KTOPIPI_KTOSIGMA_GPARITY_CMDLINE_H
#define _FIT_KTOPIPI_KTOSIGMA_GPARITY_CMDLINE_H

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE


struct CMDline{
#define DEF_ARGS(NM) \
  bool load_##NM##_data_checkpoint; \
  std::string load_##NM##_data_checkpoint_stub; \
  bool save_##NM##_data_checkpoint; \
  std::string save_##NM##_data_checkpoint_stub; \
  bool load_##NM##_amplitude_data; \
  std::string load_##NM##_amplitude_data_file; \
  bool save_##NM##_amplitude_data; \
  std::string save_##NM##_amplitude_data_file

  DEF_ARGS(ktopipi);
  DEF_ARGS(ktopipi_exc);
  DEF_ARGS(ktosigma);
#undef DEF_ARGS
 
  CMDline(){
#define INIT_ARGS(NM)				\
    load_##NM##_data_checkpoint = false;	\
    save_##NM##_data_checkpoint = false;	\
    load_##NM##_amplitude_data = false;	\
    save_##NM##_amplitude_data = false

    INIT_ARGS(ktopipi);
    INIT_ARGS(ktopipi_exc);
    INIT_ARGS(ktosigma);
#undef INIT_ARGS
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

#define PARSE_ARGS(NM) \
      }else if(sargv[i] == "-load_" #NM "_data_checkpoint"){ \
	load_##NM##_data_checkpoint = true;		     \
        load_##NM##_data_checkpoint_stub = sargv[i+1];     \
        i+=2;					     \
      }else if(sargv[i] == "-save_" #NM "_data_checkpoint"){ \
	save_##NM##_data_checkpoint = true;		     \
        save_##NM##_data_checkpoint_stub = sargv[i+1];     \
        i+=2;					     \
      }else if(sargv[i] == "-load_" #NM "_amplitude_data"){\
	load_##NM##_amplitude_data = true;\
	load_##NM##_amplitude_data_file = sargv[i+1];\
	i+=2;\
      }else if(sargv[i] == "-save_" #NM "_amplitude_data"){\
	save_##NM##_amplitude_data = true;\
	save_##NM##_amplitude_data_file = sargv[i+1];\
	i+=2

      PARSE_ARGS(ktopipi);
      PARSE_ARGS(ktopipi_exc);
      PARSE_ARGS(ktosigma);
#undef PARSE_ARGS

      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};


CPSFIT_END_NAMESPACE

#endif
