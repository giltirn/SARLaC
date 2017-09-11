#ifndef _FIT_PIPI_GPARITY_CMDLINE_H
#define _FIT_PIPI_GPARITY_CMDLINE_H

struct CMDline{
  bool load_guess;
  std::string guess_file;

  bool load_data_checkpoint;
  std::string load_data_checkpoint_file;
  bool save_data_checkpoint;
  std::string save_data_checkpoint_file;

  bool load_text_data_checkpoint;
  std::string load_text_data_checkpoint_file;
  bool save_text_data_checkpoint;
  std::string save_text_data_checkpoint_file;

  bool load_hdf5_data_checkpoint;
  std::string load_hdf5_data_checkpoint_file;
  bool save_hdf5_data_checkpoint;
  std::string save_hdf5_data_checkpoint_file;

  
  CMDline(){
    load_guess = false;
    load_data_checkpoint = false;
    save_data_checkpoint = false;
    load_text_data_checkpoint = false;
    save_text_data_checkpoint = false;
    load_hdf5_data_checkpoint = false;
    save_hdf5_data_checkpoint = false;
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
	load_data_checkpoint_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_data_checkpoint"){
	save_data_checkpoint = true;
	save_data_checkpoint_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_text_data_checkpoint"){
	load_text_data_checkpoint = true;
	load_text_data_checkpoint_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_text_data_checkpoint"){
	save_text_data_checkpoint = true;
	save_text_data_checkpoint_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_hdf5_data_checkpoint"){
	load_hdf5_data_checkpoint = true;
	load_hdf5_data_checkpoint_file = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_hdf5_data_checkpoint"){
	save_hdf5_data_checkpoint = true;
	save_hdf5_data_checkpoint_file = sargv[i+1];
	i+=2;	
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
    assert(!(load_data_checkpoint && load_text_data_checkpoint));
  }
};


#endif
