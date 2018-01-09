#ifndef _KTOPIPI_FIT_GPARITY_SCRATCH_H__
#define _KTOPIPI_FIT_GPARITY_SCRATCH_H__

class scratch{
  std::vector<std::string> scratch_files;
  const Args &_args;
  const CMDline &_cmdline;
public:

  scratch(const Args &args, const CMDline &cmdline): scratch_files(args.tsep_k_pi.size()), _args(args), _cmdline(cmdline){
    if(cmdline.use_scratch){
#ifndef HAVE_HDF5
      error_exit("getData: scratch usage requires HDF5\n");
#endif
      for(int tsep_k_pi_idx=0;tsep_k_pi_idx<args.tsep_k_pi.size();tsep_k_pi_idx++){
	std::ostringstream f; f << cmdline.use_scratch_stub << "_" << tsep_k_pi_idx;
	scratch_files[tsep_k_pi_idx] = f.str();
      }
    }
  }

  void writeScratch(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
		    std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
		    const int tsep_k_pi_idx){
    if(_cmdline.use_scratch){
      printMem("Pre-scratch write");
#ifdef HAVE_HDF5
      std::cout << "Scratch writing to " << scratch_files[tsep_k_pi_idx] << std::endl;
      boost::timer::cpu_timer timer;
      
      timer.start();
      HDF5writer writer(scratch_files[tsep_k_pi_idx]);
      write(writer, A0_all_j, "A0_all_j");
      write(writer, A0_all_dj, "A0_all_dj");
      timer.stop();
      std::cout << "Scratch write time " << timer.elapsed().wall/1e9 << "s\n";

      //Clear the data
      for(int q=0;q<10;q++){
	A0_all_j[q].clear(); A0_all_dj[q].clear();
      }
#endif
      printMem("Post-scratch write");
    }
  }

  void reloadScratch(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j,
		     std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj){
    if(_cmdline.use_scratch){
      for(int tsep_k_pi_idx=0;tsep_k_pi_idx<_args.tsep_k_pi.size();tsep_k_pi_idx++){
#ifdef HAVE_HDF5
	printMem("Pre scratch read");
	std::cout << "Scratch reading from " << scratch_files[tsep_k_pi_idx] << std::endl;
	boost::timer::cpu_timer timer;
	
	timer.start();
	HDF5reader reader(scratch_files[tsep_k_pi_idx]);
	std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > tmp_A0_all_j;
	std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > tmp_A0_all_dj;
	read(reader, tmp_A0_all_j, "A0_all_j");
	read(reader, tmp_A0_all_dj, "A0_all_dj");
	timer.stop();
	std::cout << "Scratch read time " << timer.elapsed().wall/1e9 << "s\n";

	for(int q=0;q<10;q++){
	  for(int i=0;i<tmp_A0_all_j[q].size();i++){
	    A0_all_j[q].push_back(tmp_A0_all_j[q][i]);
	    A0_all_dj[q].push_back(tmp_A0_all_dj[q][i]);
	  }
	}
	printMem("Post scratch read");
#endif
      }
    }
  }

  //We can optionally use scratch files on disk already to avoid reads
  inline bool doSkipLoad(const int tsep_k_pi_idx){
    return _cmdline.use_scratch && _cmdline.use_existing_scratch_files && fileExists(scratch_files[tsep_k_pi_idx]);
  }

};


#endif
