#if 1

//Code needs to be updated and fixed

int main(void){
  return 0;
}


#else

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
    
  std::vector<std::string> data_file_fmt_sloppy =
    { "traj_<TRAJ>_type1_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>_mom<MOM>",
      "traj_<TRAJ>_type2_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>",
      "traj_<TRAJ>_type3_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>",
      "traj_<TRAJ>_type4" };

  std::vector<std::string> data_file_fmt_exact(data_file_fmt_sloppy);
  for(int i=0;i<data_file_fmt_exact.size();i++) data_file_fmt_exact[i] = data_file_fmt_exact[i] + cmdline.symmetric_quark_momenta_figure_file_extension;

  std::vector<std::pair<threeMomentum, double> > type1_pimom_proj = {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  { {-1,1,1}, 3.0/8.0 }, { {1,-1,-1}, 3.0/8.0 }  };

  std::string bubble_file_fmt_sloppy = "traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_mom<PB>";
  std::string bubble_file_fmt_exact =  bubble_file_fmt_sloppy + "_symm";

  std::vector<std::pair<threeMomentum, double> > bubble_pimom_proj =  {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  
								         { {-1,1,1}, 1.0/8.0 }, { {1,-1,-1}, 1.0/8.0 }, 
								         { {1,-1,1}, 1.0/8.0 }, { {-1,1,-1}, 1.0/8.0 }, 
								         { {1,1,-1}, 1.0/8.0 }, { {-1,-1,1}, 1.0/8.0 } };

  if(cmdline.checkpoint_and_exit){
    checkpointRawOnly(args.tsep_k_pi, 
		      bubble_file_fmt_sloppy, bubble_file_fmt_exact, bubble_pimom_proj,
		      data_file_fmt_sloppy, data_file_fmt_exact, type1_pimom_proj,
		      args.data_dir_S, args.traj_start_S, args.traj_lessthan_S,
		      args.data_dir_C, args.traj_start_C, args.traj_lessthan_C,
		      args.traj_inc, args.bin_size, args.Lt, args.tsep_pipi, data_opt.read_opts);
    return 0;
  }

  getDataSampleAMA(A0_all_j, A0_all_dj, args.tsep_k_pi, 
		   bubble_file_fmt_sloppy, bubble_file_fmt_exact, bubble_pimom_proj,
		   data_file_fmt_sloppy, data_file_fmt_exact, type1_pimom_proj,
		   args.data_dir_S, args.traj_start_S, args.traj_lessthan_S,
		   args.data_dir_C, args.traj_start_C, args.traj_lessthan_C,
		   args.traj_inc, args.bin_size, args.Lt, args.tsep_pipi, data_opt);

  printMem("Prior to fitting");
  if(cmdline.SAMAexpand) fitAndPlotSAMAexpand(A0_all_j,A0_all_dj,args,cmdline);
  else{
    FitKtoPiPiOptions opt;
    opt.load_freeze_data = cmdline.load_freeze_data;
    opt.freeze_data = cmdline.freeze_data;    
  
    fitAndPlot(A0_all_j,A0_all_dj, args.Lt, args.tmin_k_op, args.tmin_op_pi, args.fitfunc, args.correlated, opt);
  }
  std::cout << "Done" << std::endl;
  
  return 0;
}

#endif
