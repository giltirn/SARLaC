#ifndef _FIT_KTOPIPI_GPARITY_SAMPLEAMA_CMDLINE_H
#define _FIT_KTOPIPI_GPARITY_SAMPLEAMA_CMDLINE_H

#include<fit_ktopipi_gparity/args.h>

enum SloppyExact {Sloppy, Exact};

struct SampleAMAcmdLine{
  bool load_guess;
  std::string guess_file;

  bool load_data_checkpoint_sloppy_S;
  std::string load_data_checkpoint_stub_sloppy_S; //will append  _tsepkpi<VAL>.hdf5
  bool load_data_checkpoint_sloppy_C;
  std::string load_data_checkpoint_stub_sloppy_C;
  bool load_data_checkpoint_exact_C;
  std::string load_data_checkpoint_stub_exact_C;

  bool save_data_checkpoint_sloppy_S;
  std::string save_data_checkpoint_stub_sloppy_S;
  bool save_data_checkpoint_sloppy_C;
  std::string save_data_checkpoint_stub_sloppy_C;
  bool save_data_checkpoint_exact_C;
  std::string save_data_checkpoint_stub_exact_C;

  bool load_amplitude_data;
  std::string load_amplitude_data_file;

  bool save_amplitude_data;
  std::string save_amplitude_data_file;
  
  bool load_freeze_data;
  std::string freeze_data;

  bool use_scratch; //save memory at peak times by storing to disk and reloading later
  std::string use_scratch_stub;
  bool use_existing_scratch_files; //if scratch files from a previous run exist, use these

  bool checkpoint_and_exit; //don't use the date, just load and checkpoint

  std::string symmetric_quark_momenta_figure_file_extension;

  bool SAMAexpand; //append N extra superjackknife samples to data
  int SAMAexpandN;

  bool plot_only;

  CMDline toCMDline(const char ens, const SloppyExact se) const{
    CMDline out;
    out.load_guess = load_guess;
    out.guess_file = guess_file;

    out.load_amplitude_data = load_amplitude_data;
    out.load_amplitude_data_file = load_amplitude_data_file;
    out.save_amplitude_data = save_amplitude_data;
    out.save_amplitude_data_file = save_amplitude_data_file;
    out.load_freeze_data = load_freeze_data;
    out.freeze_data = freeze_data;
    out.use_scratch = use_scratch;
    out.use_existing_scratch_files = use_existing_scratch_files;
    out.symmetric_quark_momenta_figure_file_extension = symmetric_quark_momenta_figure_file_extension;

    if(ens == 'S'){
      assert(se == Sloppy);
      out.load_data_checkpoint = load_data_checkpoint_sloppy_S;
      out.save_data_checkpoint = save_data_checkpoint_sloppy_S;
      out.load_data_checkpoint_stub = load_data_checkpoint_stub_sloppy_S;
      out.save_data_checkpoint_stub = save_data_checkpoint_stub_sloppy_S;
      out.use_scratch_stub = use_scratch_stub + "_sloppy_S";
      out.use_symmetric_quark_momenta = false;
    }else{
      out.load_data_checkpoint = se == Sloppy ? load_data_checkpoint_sloppy_C : load_data_checkpoint_exact_C;     
      out.save_data_checkpoint = se == Sloppy ? save_data_checkpoint_sloppy_C : save_data_checkpoint_exact_C;     
      out.load_data_checkpoint_stub = se == Sloppy ? load_data_checkpoint_stub_sloppy_C : load_data_checkpoint_stub_exact_C;
      out.save_data_checkpoint_stub = se == Sloppy ? save_data_checkpoint_stub_sloppy_C : save_data_checkpoint_stub_exact_C;
      out.use_scratch_stub = use_scratch_stub + (se == Sloppy ? "_sloppy_C" : "_exact_C");
      out.use_symmetric_quark_momenta = (se == Exact);
    }
    return out;
  }


  SampleAMAcmdLine(){
    load_guess = false;
    load_data_checkpoint_sloppy_S = false;
    load_data_checkpoint_sloppy_C = false;
    load_data_checkpoint_exact_C = false;
    save_data_checkpoint_sloppy_S = false;
    save_data_checkpoint_sloppy_C = false;
    save_data_checkpoint_exact_C = false;
    load_amplitude_data = false;
    save_amplitude_data = false;
    load_freeze_data = false;
    use_scratch = false;
    use_existing_scratch_files = false;
    use_scratch_stub = "scratch";
    checkpoint_and_exit = false;
    symmetric_quark_momenta_figure_file_extension = "_symm";
    SAMAexpand = false;
    plot_only = false;
  }
  SampleAMAcmdLine(const int argc, const char** argv, const int begin = 0): SampleAMAcmdLine(){
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
      }else if(sargv[i] == "-load_data_checkpoint_sloppy_S"){
	load_data_checkpoint_sloppy_S = true;
	load_data_checkpoint_stub_sloppy_S = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_data_checkpoint_sloppy_C"){
	load_data_checkpoint_sloppy_C = true;
	load_data_checkpoint_stub_sloppy_C = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-load_data_checkpoint_exact_C"){
	load_data_checkpoint_exact_C = true;
	load_data_checkpoint_stub_exact_C = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_data_checkpoint_sloppy_S"){
	save_data_checkpoint_sloppy_S = true;
	save_data_checkpoint_stub_sloppy_S = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_data_checkpoint_sloppy_C"){
	save_data_checkpoint_sloppy_C = true;
	save_data_checkpoint_stub_sloppy_C = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-save_data_checkpoint_exact_C"){
	save_data_checkpoint_exact_C = true;
	save_data_checkpoint_stub_exact_C = sargv[i+1];
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
      }else if(sargv[i] == "-superjack_maxmem"){ //in GB
	std::stringstream ss; ss << sargv[i+1];
	size_t bytes;  ss >> bytes;
	const size_t k = 1024;
	bytes *= k*k*k;					  
	constrainedMemoryManager::maxSize() = bytes;
	i+=2;
      }else if(sargv[i] == "-checkpoint_and_exit"){
	checkpoint_and_exit = true;
	i++;
      }else if(sargv[i] == "-symmetric_quark_momenta_figure_file_extension"){
	symmetric_quark_momenta_figure_file_extension = sargv[i+1];
	i+=2;
      }else if(sargv[i] == "-SAMAexpand"){
	SAMAexpand = true;
	SAMAexpandN = strToAny<int>(sargv[i+1]);
	i+=2;
      }else if(sargv[i] == "-plot_only"){
	plot_only = true;
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


#endif
