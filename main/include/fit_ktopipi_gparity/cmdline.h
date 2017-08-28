#ifndef _FIT_KTOPIPI_GPARITY_CMDLINE_H
#define _FIT_KTOPIPI_GPARITY_CMDLINE_H

#include<fit_ktopipi_gparity/args.h>

struct CMDline{
  bool load_guess;
  std::string guess_file;

  bool load_data_checkpoint;
  std::string load_data_checkpoint_stub; //will append  _tsepkpi<VAL>.hdf5
  bool save_data_checkpoint;
  std::string save_data_checkpoint_stub;

  bool load_freeze_data;
  std::string freeze_data;
  
  CMDline(){
    load_guess = false;
    load_data_checkpoint = false;
    save_data_checkpoint = false;
    load_freeze_data = false;
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
      }else if(sargv[i] == "-freeze"){
	load_freeze_data = true;
	freeze_data = sargv[i+1];
	if(freeze_data == "TEMPLATE"){
	  FreezeParams p;
	  std::ofstream of("freeze_template.args");
	  of << p;
	  of.close();
	  std::cout << "Wrote freeze template argument file to freeze_template.args\n";
	  exit(0);
	}
	i+=2;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};


#endif
