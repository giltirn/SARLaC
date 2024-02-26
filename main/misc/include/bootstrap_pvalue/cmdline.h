#pragma once

struct CMDline{
  bool exit_after_preanalysis;

  int seed; //for main rng
  int seed_thr; //for thread rng

  bool recenter_orig_ens; //this should be on, but can be disabled to demonstrate its effect

  bool bootstrap_resid_diagonalize; //do diagonalization of residuals prior to resampling
  bool bootstrap_resid_diagonalize_evals; //do diagonalization of residuals and normalize by rooted-inverse-evals prior to resampling

  bool output_bootstrap_q2sorted_rtable;  //for the primary bootstrap analysis, write out the rtable sorted by q^2 in ascending order, with the value of q^2 as the first entry on each line
  CMDline(){
    exit_after_preanalysis = false;
    bootstrap_resid_diagonalize = false;
    bootstrap_resid_diagonalize_evals = false;
    recenter_orig_ens = true;
    output_bootstrap_q2sorted_rtable = false;

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
      if(sargv[i] == "-exit_after_preanalysis"){
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
      }else if(sargv[i] == "-bootstrap_resid_diagonalize_evals"){
	bootstrap_resid_diagonalize_evals = true;
	i++;
      }else if(sargv[i] == "-no_recenter_orig_ens"){
	recenter_orig_ens = false;
	i++;
      }else if(sargv[i] == "-output_bootstrap_q2sorted_rtable"){
	output_bootstrap_q2sorted_rtable = true;
	i++;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};
