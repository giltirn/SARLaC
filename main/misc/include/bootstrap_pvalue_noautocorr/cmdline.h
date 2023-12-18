#pragma once

struct CMDline{
  bool write_data;
  bool exit_after_preanalysis;

  CMDline(){
    write_data = false;
    exit_after_preanalysis = false;
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
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};
