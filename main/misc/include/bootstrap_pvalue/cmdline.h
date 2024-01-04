#pragma once

struct CMDline{
  bool write_data;
  bool exit_after_preanalysis;

  int seed; //for main rng
  int seed_thr; //for thread rng

  bool bootstrap_resid_diagonalize; //do diagonalization of residuals prior to resampling

  CMDline(){
    write_data = false;
    exit_after_preanalysis = false;
    bootstrap_resid_diagonalize = false;

    seed = 1234;
    seed_thr = 5678;
  }
  CMDline(const int argc, const char** argv, const int begin = 0): CMDline(){
    setup(argc,argv,begin);
  }
  
  void setup(const int argc, const char** argv, const int begin = 0){
    const int sz = argc-begin;
    if(sz <= 0) return;

    std::vector<std::string> sargv(sz);
    for(int i=begin; i<argc; i++) sargv[i-begin] = std::string(argv[i]);

    int i = 0;
    while(i<sz){
      if(sargv[i] == "-write_data"){
	write_data = true;
	i++;
      }else if(sargv[i] == "-exit_after_preanalysis"){
	exit_after_preanalysis = true;;
	i++;
      }else if(sargv[i] == "-seed"){
	seed = strToAny<int>(sargv[i+1]);
	seed_thr = strToAny<int>(sargv[i+2]);
	std::cout << "Set seeds to " << seed << " " << seed_thr << std::endl;
	i+=3;
      }else if(sargv[i] == "-bootstrap_resid_diagonalize"){
	bootstrap_resid_diagonalize = true;
	i++;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};
