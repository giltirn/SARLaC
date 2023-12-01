#include<ktopipi_common/ktopipi_common.h>

#include<fit_ktosigma_gparity/args.h>
#include<fit_ktosigma_gparity/cmdline.h>

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

  const std::vector<std::pair<threeMomentum, double> > bubble_quarkmom_proj = { { {1,1,1}, 1./8 }, { {-1,-1,-1}, 1./8 },
										{ {-3,1,1}, 1./8 }, { {3,-1,-1}, 1./8 },
										{ {1,-3,1}, 1./8 }, { {-1,3,-1}, 1./8 },
										{ {1,1,-3}, 1./8 }, { {-1,-1,3}, 1./8 }       };
  readKtoPiPiAllDataOptions read_opt;
  read_opt.import(cmdline);

  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > A0_all_dj(10);
  getKtoSigmaData(A0_all_j, A0_all_dj, args.tsep_k_sigma, args.data_dir, args.type_file_fmt,  args.bubble_file_fmt, 
		  bubble_quarkmom_proj, args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size, args.Lt, read_opt);
  
  printMem("Prior to fitting");

  FitKtoPiPiOptions opt;
  opt.load_freeze_data = cmdline.load_freeze_data;
  opt.freeze_data = cmdline.freeze_data;
  
  fitAndPlot(A0_all_j,A0_all_dj, args.Lt, args.tmin_k_op, args.tmin_op_sigma, args.fitfunc, args.correlated, opt);
  
  std::cout << "Done" << std::endl;
  
  return 0;
}
