#include<distribution.h>
#include<common.h>
#include<physics.h>
#include<parser.h>

using namespace CPSfit;

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

  //Load the Wilson coefficients
  WilsonCoeffs wilson_coeffs(args.wilson_coeffs_file);
  
  //Compute the MSbar matching factors
  MSbarConvert MSbar(args.renormalization.mu,args.renormalization.scheme);

  std::pair<NumericTensor<double,1>, NumericTensor<double,1> > RI_Wilson_coeffs;
  std::pair<NumericTensor<superJackD,1>, NumericTensor<superJackD,1> > lat_Wilson_coeffs;
  
  computeRIandLatticeWilsonCoefficients(RI_Wilson_coeffs,lat_Wilson_coeffs,wilson_coeffs,MSbar,data.NPR_sj,"");
  
  //Compute the Lellouch-Luscher factor
  superJackD F_sj, delta_0_sj;
  computePhaseShiftAndLLfactor(delta_0_sj, F_sj, data.Epipi_sj, data.mpi_sj, data.mK_sj, data.ainv_sj, args);

  {
    std::cout << "\n\n---------------------------------------------------------------\n";
    std::cout << "Starting computation using standard method\n";
    
    //Compute the MSbar, infinite-volume matrix elements
    NumericTensor<superJackD,1> M_MSbar_std_sj = computePhysicalMSbarMatrixElements(data.M_lat_sj,data.ainv_sj,F_sj,data.NPR_sj,MSbar,args);

    computeMatrixElementRelations(M_MSbar_std_sj, "");
    
    //Compute A0
    superJackD ReA0_sj, ImA0_sj;
    computeA0(ReA0_sj, ImA0_sj, M_MSbar_std_sj, wilson_coeffs, args, "");
    
    //Compute epsilon'
    computeEpsilon(ReA0_sj, ImA0_sj,
		   data.ReA2_lat_sj, data.ImA2_lat_sj,
		   delta_0_sj, data.delta_2_sj,
		   data.ReA0_expt_sj, data.ReA2_expt_sj,
		   data.omega_expt_sj, data.mod_eps_sj,
		   args, "");


  }
  
  for(int i=0;i<7;i++){
    std::cout << "\n\n---------------------------------------------------------------\n";
    std::cout << "Starting computation using Re(A0) to eliminate Q'" << chiralBasisIdx(i) << std::endl;

    std::ostringstream nm; nm << "elim" << chiralBasisIdx(i) << "_";
        
    //Compute the MSbar, infinite-volume matrix elements
    NumericTensor<superJackD,1> M_MSbar_std_sj = computePhysicalMSbarMatrixElementsUsingReA0expt(data.M_lat_sj,data.ainv_sj,F_sj,data.NPR_sj,MSbar,lat_Wilson_coeffs.first,data.ReA0_expt_sj,i,args);

    computeMatrixElementRelations(M_MSbar_std_sj, nm.str());
    
    //Compute A0
    superJackD ReA0_sj, ImA0_sj;
    computeA0(ReA0_sj, ImA0_sj, M_MSbar_std_sj, wilson_coeffs, args, nm.str());
    
    //Compute epsilon'
    computeEpsilon(ReA0_sj, ImA0_sj,
		   data.ReA2_lat_sj, data.ImA2_lat_sj,
		   delta_0_sj, data.delta_2_sj,
		   data.ReA0_expt_sj, data.ReA2_expt_sj,
		   data.omega_expt_sj, data.mod_eps_sj,
		   args, nm.str());
  }

    
  
  std::cout << "Done\n";
  return 0;
}
