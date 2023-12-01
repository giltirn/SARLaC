#include<ktopipi_common/ktopipi_common.h>

using namespace SARLaC;

#include <fit_ktopipi_gparity/cmdline.h>
#include <fit_ktopipi_gparity/args.h>

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
  
  //Prepare the data
  readKtoPiPiAllDataOptions read_opt;
  read_opt.import(cmdline);

  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > A0_all_dj(10);
  printMem("Prior to getData");

  getData(A0_all_j, A0_all_dj, args.tsep_k_pi, args.data_dir, 
	  args.type_file_fmt, args.type1_pimom_proj, 
	  args.bubble_file_fmt, args.bubble_pimom_proj,
	  args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size, args.Lt, args.tsep_pipi, read_opt);
  
  printMem("Prior to fitting");

  FitKtoPiPiOptions opt;
  opt.load_freeze_data = cmdline.load_freeze_data;
  opt.freeze_data = cmdline.freeze_data;  

  fitAndPlot(A0_all_j,A0_all_dj, args.Lt, args.tmin_k_op, args.tmin_op_pi, args.fitfunc, args.correlated, opt);
  
  std::cout << "Done" << std::endl;
  
  return 0;
}

