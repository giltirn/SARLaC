#include<ktopipi_common/ktopipi_common.h>
#include<ktopipi_sampleAMA_common/ktopipi_sampleAMA_common.h>

using namespace CPSfit;

#include <fit_ktopipi_gparity_sampleAMA/cmdline.h>
#include <fit_ktopipi_gparity_sampleAMA/args.h>
#include <fit_ktopipi_gparity_sampleAMA/fit_sama_expand.h>

int main(const int argc, const char* argv[]){
  printMem("Beginning of execution");
  
  SampleAMAargs args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  const std::string arg_file = argv[1];
  parse(args, arg_file);
  
  SampleAMAcmdLine cmdline(argc,argv,2);
  freezeCheck<>(cmdline);

  readKtoPiPiAllDataSampleAMAoptions data_opt;
  data_opt.importGlobalOptions(cmdline);
  data_opt.read_opts = cmdline.getSampleAMAreadOptions();

  if(cmdline.plot_only){
    std::cout << "Plotting results and exiting\n";
    assert(cmdline.load_amplitude_data);
    PlotOnlyInputs pi(args.Lt, args.tmin_k_op, args.tmin_op_pi, cmdline.load_amplitude_data_file);

    fitfuncCall<PlotOnlyInputs,PlotOnlyCall>(args.fitfunc,pi);
    std::cout << "Done" << std::endl;
    return 0;
  }

  //Prepare the data
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > A0_all_dj(10);
  printMem("Prior to getData");
    
  if(cmdline.checkpoint_and_exit){
    checkpointRawOnly(args.tsep_k_pi, 
		      args.data_dir_S, args.traj_start_S, args.traj_lessthan_S,
		      args.data_dir_C, args.traj_start_C, args.traj_lessthan_C,
		      args.traj_inc, args.bin_size, args.Lt, args.tsep_pipi, data_opt.read_opts);
    return 0;
  }

  getDataSampleAMA(A0_all_j, A0_all_dj, args.tsep_k_pi, 
		   args.data_dir_S, args.traj_start_S, args.traj_lessthan_S,
		   args.data_dir_C, args.traj_start_C, args.traj_lessthan_C,
		   args.traj_inc, args.bin_size, args.Lt, args.tsep_pipi, data_opt);

  printMem("Prior to fitting");
  if(cmdline.SAMAexpand) fitAndPlotSAMAexpand(A0_all_j,A0_all_dj,args,cmdline);
  else fitAndPlot(A0_all_j,A0_all_dj, args.Lt, args.tmin_k_op, args.tmin_op_pi, args.fitfunc, args.correlated, cmdline.load_freeze_data, cmdline.freeze_data);
  
  std::cout << "Done" << std::endl;
  
  return 0;
}

