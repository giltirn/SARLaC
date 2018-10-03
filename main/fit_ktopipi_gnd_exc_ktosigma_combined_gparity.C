#include<ktopipi_common/ktopipi_common.h>

#include<fit_ktopipi_ktosigma_combined_gparity/simfit_generic.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/args.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/cmdline.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/simfit_plot.h>
#include<fit_ktopipi_gnd_exc_ktosigma_combined_gparity/simfit.h>

using namespace CPSfit;


int main(const int argc, const char* argv[]){
  printMem("Beginning of execution");
  
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  const std::string arg_file = argv[1];
  parse(args, arg_file);

  CMDline cmdline(argc,argv,2); 

  int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc/args.bin_size;

  simultaneousFitBase* fitter = args.fitfunc == SimFitFunction::TwoState ? (simultaneousFitBase*)(new simultaneousFit3state) : (simultaneousFitBase*)(new simultaneousFit3state);

  fitter->load2ptFitParams(args.input_params, nsample); 

  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > ktopipi_A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > ktopipi_A0_all_dj(10);

#define COPYOPTS(INTO, NM) \
  INTO.load_amplitude_data = cmdline.load_##NM##_amplitude_data; \
  INTO.load_amplitude_data_file = cmdline.load_##NM##_amplitude_data_file; \
  INTO.save_amplitude_data = cmdline.save_##NM##_amplitude_data; \
  INTO.save_amplitude_data_file = cmdline.save_##NM##_amplitude_data_file; \
  INTO.read_opts.load_data_checkpoint = cmdline.load_##NM##_data_checkpoint; \
  INTO.read_opts.load_data_checkpoint_stub = cmdline.load_##NM##_data_checkpoint_stub; \
  INTO.read_opts.save_data_checkpoint = cmdline.save_##NM##_data_checkpoint; \
  INTO.read_opts.save_data_checkpoint_stub = cmdline.save_##NM##_data_checkpoint_stub


  readKtoPiPiAllDataOptions read_opts;
  COPYOPTS(read_opts, ktopipi);

  getData(ktopipi_A0_all_j, ktopipi_A0_all_dj, args.tsep_k_pi, args.data_dir, 
	  args.ktopipi_type_file_fmt, args.ktopipi_type1_pimom_proj, 
	  args.pipi_bubble_file_fmt, args.pipi_bubble_pimom_proj,
	  args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size, args.Lt, args.tsep_pipi, read_opts);


  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > ktopipi_exc_A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > ktopipi_exc_A0_all_dj(10);

  COPYOPTS(read_opts, ktopipi_exc);
    
  getData(ktopipi_exc_A0_all_j, ktopipi_exc_A0_all_dj, args.tsep_k_pi, args.data_dir, 
	  args.ktopipi_exc_type_file_fmt, args.ktopipi_exc_type1_pimom_proj, 
	  args.pipi_bubble_file_fmt, args.pipi_exc_bubble_pimom_proj,
	  args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size, args.Lt, args.tsep_pipi, read_opts);

  const std::vector<std::pair<threeMomentum, double> > sigma_bub_quarkmom_proj = { { {1,1,1}, 1./8 }, { {-1,-1,-1}, 1./8 },
										{ {-3,1,1}, 1./8 }, { {3,-1,-1}, 1./8 },
										{ {1,-3,1}, 1./8 }, { {-1,3,-1}, 1./8 },
										{ {1,1,-3}, 1./8 }, { {-1,-1,3}, 1./8 }       };


  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > ktosigma_A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > ktosigma_A0_all_dj(10);

  COPYOPTS(read_opts, ktosigma);

  getKtoSigmaData(ktosigma_A0_all_j, ktosigma_A0_all_dj, args.tsep_k_sigma, args.data_dir, args.ktosigma_type_file_fmt,  args.sigma_bubble_file_fmt, 
		  sigma_bub_quarkmom_proj, args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size, args.Lt, read_opts);


  fitter->fit(ktopipi_A0_all_j, ktopipi_A0_all_dj, 
	     ktopipi_exc_A0_all_j, ktopipi_exc_A0_all_dj, 
	     ktosigma_A0_all_j, ktosigma_A0_all_dj,
	     args.Lt, args.tmin_k_op, args.tmin_op_snk, args.correlated);

  std::cout << "Done" << std::endl;
  
  return 0;
}
