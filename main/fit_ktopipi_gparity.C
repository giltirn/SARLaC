#include<ktopipi_common/ktopipi_common.h>

using namespace CPSfit;

#include <fit_ktopipi_gparity/cmdline.h>
#include <fit_ktopipi_gparity/args.h>


// inline std::string typeFile(const int traj, const int type, const int tsep_k_pi, const int tsep_pipi, const std::string &data_dir, const bool use_symmetric_quark_momenta, const std::string &symmetric_quark_momenta_figure_file_extension, const threeMomentum &mom = {0,0,0}){
//   std::ostringstream os;
//   os << data_dir << "/traj_" << traj << "_type" << type;
//   if(type != 4) os << "_deltat_" << tsep_k_pi << "_sep_" << tsep_pipi;
//   if(type == 1) os << "_mom" << momStr(mom);
//   if(use_symmetric_quark_momenta) os << symmetric_quark_momenta_figure_file_extension;
//   return os.

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

  std::vector<std::string> data_file_fmt =
    { "traj_<TRAJ>_type1_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>_mom<MOM>",
      "traj_<TRAJ>_type2_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>",
      "traj_<TRAJ>_type3_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>",
      "traj_<TRAJ>_type4" };

  std::vector<std::pair<threeMomentum, double> > type1_pimom_proj = {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  { {-1,1,1}, 3.0/8.0 }, { {1,-1,-1}, 3.0/8.0 }  };

  std::string bubble_file_fmt = "traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_mom<PB>";
  std::vector<std::pair<threeMomentum, double> > bubble_pimom_proj =  {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  
								         { {-1,1,1}, 1.0/8.0 }, { {1,-1,-1}, 1.0/8.0 }, 
								         { {1,-1,1}, 1.0/8.0 }, { {-1,1,-1}, 1.0/8.0 }, 
								         { {1,1,-1}, 1.0/8.0 }, { {-1,-1,1}, 1.0/8.0 } };

  getData(A0_all_j, A0_all_dj, args.tsep_k_pi, args.data_dir, 
	  data_file_fmt, type1_pimom_proj, 
	  bubble_file_fmt, bubble_pimom_proj,
	  args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size, args.Lt, args.tsep_pipi, read_opt);
  
  printMem("Prior to fitting");
  fitAndPlot(A0_all_j,A0_all_dj, args.Lt, args.tmin_k_op, args.tmin_op_pi, args.fitfunc, args.correlated, cmdline.load_freeze_data, cmdline.freeze_data);
  
  std::cout << "Done" << std::endl;
  
  return 0;
}

