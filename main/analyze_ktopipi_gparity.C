#include<distribution.h>
#include<common_defs.h>
#include<superjackknife.h>
#include<alpha_s.h>
#include<parser.h>

#include<analyze_ktopipi_gparity/params.h>
#include<analyze_ktopipi_gparity/LLfactor.h>
#include<analyze_ktopipi_gparity/renormalization.h>
#include<analyze_ktopipi_gparity/WilsonCoeffs.h>
#include<analyze_ktopipi_gparity/utils.h>
#include<analyze_ktopipi_gparity/read_data.h>
#include<analyze_ktopipi_gparity/main.h>


int main(const int argc, const char* argv[]){
  typedef superJackknifeDistribution<double> superJackD;  

  if(argc != 2){
    std::cout << "To use ./analyze_ktopipi_gparity <params_file>\n";
    std::cout << "Writing template to template.args\n";
    Args tp;
    std::ofstream of("template.args");
    of << tp;
    of.close();    
    exit(0);
  }
  Args args;
  parse(args, argv[1]);

  //Load the data and convert into common superJackknife format
  superJackknifeData data(args);  

  //Compute the Lellouch-Luscher factor
  superJackD F_sj, delta_0_sj;
  computePhaseShiftAndLLfactor(delta_0_sj, F_sj, data.Epipi_sj, data.mpi_sj, data.mK_sj, data.ainv_sj, args);

  //Compute the MSbar, infinite-volume matrix elements
  NumericTensor<superJackD,1> M_MSbar_std_sj = computePhysicalMSbarMatrixElements(data.M_lat_sj,data.ainv_sj,F_sj,data.NPR_sj,args);

  //Compute A0
  superJackD ReA0_sj, ImA0_sj;
  computeA0(ReA0_sj, ImA0_sj, M_MSbar_std_sj, args, "");

  //Compute epsilon'
  computeEpsilon(ReA0_sj, ImA0_sj,
		 data.ReA2_lat_sj, data.ImA2_lat_sj,
		 delta_0_sj, data.delta_2_sj,
		 data.ReA0_expt_sj, data.ReA2_expt_sj,
		 data.omega_expt_sj, data.mod_eps_sj,
		 args, "");
  
  std::cout << "Done\n";
  return 0;
}
