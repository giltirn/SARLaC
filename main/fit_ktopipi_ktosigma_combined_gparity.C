#include<ktopipi_common/ktopipi_common.h>

#include<fit_ktopipi_ktosigma_combined_gparity/simfit_generic.h>
#include<fit_ktopipi_ktosigma_combined_gparity/args.h>
#include<fit_ktopipi_ktosigma_combined_gparity/cmdline.h>
#include<fit_ktopipi_ktosigma_combined_gparity/projfit.h>
#include<fit_ktopipi_ktosigma_combined_gparity/simfit.h>

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

  jackknifeDistributionD mK;
  jackknifeDistributionD cK;
  {
    std::vector<jackknifeDistributionD> p;
    readParamsStandard(p,  args.input_params.kaon2pt_fit_result);
    mK = p[args.input_params.idx_mK];
    cK = sqrt( p[args.input_params.idx_cK] );
  }

  int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc/args.bin_size;
  jackknifeDistributionD E0, E1;
  NumericSquareMatrix<jackknifeDistributionD> coeffs(2); //row = (0=pipi, 1=sigma)  col = (0=gnd state, 1=exc state)
  {
    double scale = sqrt(args.input_params.pipi_sigma_sim_fit_Ascale);
    std::vector<jackknifeDistributionD> p;
    readParamsStandard(p, args.input_params.pipi_sigma_sim_fit_result);
    for(int i=0;i<p.size();i++) assert(p[i].size() == nsample);
    coeffs(0,0) = p[args.input_params.idx_coeff_pipi_state0] * scale;
    coeffs(0,1) = p[args.input_params.idx_coeff_pipi_state1] * scale;
    coeffs(1,0) = p[args.input_params.idx_coeff_sigma_state0] * scale;
    coeffs(1,1) = p[args.input_params.idx_coeff_sigma_state1] * scale;
    E0 = p[args.input_params.idx_E0];
    E1 = p[args.input_params.idx_E1];
  }

  std::cout << "cK = " << cK << std::endl;
  std::cout << "mK = " << mK << std::endl;
  std::cout << "E0 = " << E0 << std::endl;
  std::cout << "E1 = " << E1 << std::endl;

  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > ktopipi_A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > ktopipi_A0_all_dj(10);

  readKtoPiPiAllDataOptions kpp_opts;
  kpp_opts.load_amplitude_data = cmdline.load_ktopipi_amplitude_data;
  kpp_opts.load_amplitude_data_file = cmdline.load_ktopipi_amplitude_data_file; 
  kpp_opts.save_amplitude_data = cmdline.save_ktopipi_amplitude_data;
  kpp_opts.save_amplitude_data_file = cmdline.save_ktopipi_amplitude_data_file;
  kpp_opts.read_opts.load_data_checkpoint = cmdline.load_ktopipi_data_checkpoint;
  kpp_opts.read_opts.load_data_checkpoint_stub = cmdline.load_ktopipi_data_checkpoint_stub;
  kpp_opts.read_opts.save_data_checkpoint = cmdline.save_ktopipi_data_checkpoint;
  kpp_opts.read_opts.save_data_checkpoint_stub = cmdline.save_ktopipi_data_checkpoint_stub;


  getData(ktopipi_A0_all_j, ktopipi_A0_all_dj, args.tsep_k_pi, args.data_dir, 
	  args.ktopipi_type_file_fmt, args.ktopipi_type1_pimom_proj, 
	  args.pipi_bubble_file_fmt, args.pipi_bubble_pimom_proj,
	  args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size, args.Lt, args.tsep_pipi, kpp_opts);


  const std::vector<std::pair<threeMomentum, double> > sigma_bub_quarkmom_proj = { { {1,1,1}, 1./8 }, { {-1,-1,-1}, 1./8 },
										{ {-3,1,1}, 1./8 }, { {3,-1,-1}, 1./8 },
										{ {1,-3,1}, 1./8 }, { {-1,3,-1}, 1./8 },
										{ {1,1,-3}, 1./8 }, { {-1,-1,3}, 1./8 }       };


  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > ktosigma_A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > ktosigma_A0_all_dj(10);

  readKtoPiPiAllDataOptions ks_opts;
  ks_opts.load_amplitude_data = cmdline.load_ktosigma_amplitude_data;
  ks_opts.load_amplitude_data_file = cmdline.load_ktosigma_amplitude_data_file; //"data_ktosigma.hdf5";
  ks_opts.save_amplitude_data = cmdline.save_ktosigma_amplitude_data;
  ks_opts.save_amplitude_data_file = cmdline.save_ktosigma_amplitude_data_file;
  ks_opts.read_opts.load_data_checkpoint = cmdline.load_ktosigma_data_checkpoint;
  ks_opts.read_opts.load_data_checkpoint_stub = cmdline.load_ktosigma_data_checkpoint_stub;
  ks_opts.read_opts.save_data_checkpoint = cmdline.save_ktosigma_data_checkpoint;
  ks_opts.read_opts.save_data_checkpoint_stub = cmdline.save_ktosigma_data_checkpoint_stub;


  getKtoSigmaData(ktosigma_A0_all_j, ktosigma_A0_all_dj, args.tsep_k_sigma, args.data_dir, args.ktosigma_type_file_fmt,  args.sigma_bubble_file_fmt, 
		  sigma_bub_quarkmom_proj, args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size, args.Lt, ks_opts);



  if(!args.do_simfit){
    projectionFit(ktopipi_A0_all_j, ktopipi_A0_all_dj, ktosigma_A0_all_j, ktosigma_A0_all_dj,
		  mK, cK, E0, E1, coeffs, args.Lt, args.tmin_k_op, args.tmin_op_snk, args.correlated);

  }else{
    simultaneousFit(ktopipi_A0_all_j, ktopipi_A0_all_dj, ktosigma_A0_all_j, ktosigma_A0_all_dj,
		    mK, cK, E0, E1, coeffs, args.Lt, args.tmin_k_op, args.tmin_op_snk, args.correlated);
  }

  std::cout << "Done" << std::endl;
  
  return 0;
}
