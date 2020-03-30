#include<distribution.h>
#include<common.h>
#include<physics.h>
#include<parser.h>
#include<fit.h>
#include<random.h>

using namespace CPSfit;

#include<analyze_ktopipi_gparity/params.h>
#include<analyze_ktopipi_gparity/LLfactor.h>
#include<analyze_ktopipi_gparity/renormalization.h>
#include<analyze_ktopipi_gparity/WilsonCoeffs.h>
#include<analyze_ktopipi_gparity/utils.h>
#include<analyze_ktopipi_gparity/read_data.h>
#include<analyze_ktopipi_gparity/convert_chiral.h>
#include<analyze_ktopipi_gparity/main.h>

typedef superMultiDistribution<double> superMultiD;  

jackknifeDistributionD randomJackknife(double cen, double err, double err_tol, int N, int iter_max = 1000, int iter = 0){
  jackknifeDistributionD out(N,cen);
  if(err == 0.0) return out;
  
  double dev = err / sqrt(double(N-1));
  for(int i=0;i<N;i++)
    out.sample(i) = gaussianRandom<double>(cen, dev);

  //Shift slightly to correct mean
  double delta = cen - out.mean();
  for(int i=0;i<N;i++)
    out.sample(i) += delta;

  //Make it recursive to get the error correct to within tol
  if( abs(  (out.standardError() - err) / err  ) > err_tol ){
    if(iter >= iter_max){
      throw "Error: Could not generate distribution with desired tolerance";
    }
    return randomJackknife(cen, err, err_tol, N, iter_max, iter+1);
  }
  return out;
}


//Must have added REA0_SYS and IMA0_SYS entries to supermulti layout
void addA0syserr(superMultiD &ReA0_sj, superMultiD &ImA0_sj,
		 const double reA0_sysfrac, const double imA0_sysfrac){
  const superMultiLayout &layout = ReA0_sj.getLayout();
  int reA0_sys_ensidx = layout.ensIdx("REA0_SYS");
  int imA0_sys_ensidx = layout.ensIdx("IMA0_SYS");
  assert( reA0_sys_ensidx !=0 && imA0_sys_ensidx !=0);
  
  int reA0_sys_sz = layout.nSamplesEns(reA0_sys_ensidx);
  int imA0_sys_sz = layout.nSamplesEns(imA0_sys_ensidx);
  
  double reA0_sys = fabs(reA0_sysfrac * ReA0_sj.best());
  double imA0_sys = fabs(imA0_sysfrac * ImA0_sj.best());
  
  std::cout << "ReA0 prior to including systematic error of size " << reA0_sys << " : " << ReA0_sj << std::endl;
  jackknifeDistributionD reA0_sys_j = randomJackknife(ReA0_sj.best(), reA0_sys, 1e-3, reA0_sys_sz, 1e6);
  std::cout << "Generated " << reA0_sys_j << std::endl;
  ReA0_sj.setEnsembleDistribution(reA0_sys_ensidx, reA0_sys_j);
  std::cout << "ReA0 after including systematic error : " << ReA0_sj << std::endl;

  std::cout << "ImA0 prior to including systematic error of size " << imA0_sys << " : " << ImA0_sj << std::endl;
  jackknifeDistributionD imA0_sys_j = randomJackknife(ImA0_sj.best(), imA0_sys, 1e-3, imA0_sys_sz, 1e6);
  std::cout << "Generated " << imA0_sys_j << std::endl;
  ImA0_sj.setEnsembleDistribution(imA0_sys_ensidx, imA0_sys_j);
  std::cout << "ImA0 after including systematic error : " << ImA0_sj << std::endl;
}

				  
  


  

int main(const int argc, const char* argv[]){
  RNG.initialize(1234);

  if(argc < 2){
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

  bool zero_all_but_NPR_errs = false;
  bool add_A0_sys = false;
  double reA0_sys_frac, imA0_sys_frac;

  int arg = 2;
  while(arg < argc){
    std::string sargv(argv[arg]);
    if(sargv == "-zero_all_but_NPR_errs"){
      zero_all_but_NPR_errs = true;
      arg++;
    }else if(sargv == "-add_A0_sys"){ //include A0 systematic in epsilon' determination
      add_A0_sys = true;
      reA0_sys_frac = strToAny<double>(argv[arg+1]);
      imA0_sys_frac = strToAny<double>(argv[arg+2]);
      arg+=3;
    }else{
      error_exit(std::cout << "Unknown option " << sargv << std::endl);
    }
  }

  //Load the data and convert into common superMulti format
  superMultiData data(args, add_A0_sys ? 50 : 0, add_A0_sys ? 50 : 0);  

  if(zero_all_but_NPR_errs) data.zeroAllButNPRerrs();

  //Setup the input parameters to the perturbation theory
  PerturbativeVariables pv(args.perturbative_inputs);

  //Compute the Wilson coefficients
  WilsonCoeffs wilson_coeffs(args.renormalization.mu, pv, args.constants.tau);
  wilson_coeffs.writeLatex("MSbar_Wilson_coeffs.dat");

  //Compute the MSbar matching factors
  MSbarConvert MSbar(args.renormalization.mu,args.renormalization.scheme,pv);
  
  {
    //Compute lattice -> MSbar renormalization matrix
    NumericTensor<superMultiD,2> NPRlattoMSbar(data.NPR_sj);
    for(int i=0;i<7;i++)
      for(int j=0;j<7;j++){
	NPRlattoMSbar({i,j}) = MSbar.getConversionMatrix(Chiral)({i,0}) * data.NPR_sj({0,j});
	for(int k=1;k<7;k++) NPRlattoMSbar({i,j}) = NPRlattoMSbar({i,j}) + MSbar.getConversionMatrix(Chiral)({i,k}) * data.NPR_sj({k,j});
      }
    writeParamsStandard(NPRlattoMSbar, "NPR_lat_to_MSbar.hdf5");
  }

  std::pair<NumericTensor<double,1>, NumericTensor<double,1> > RI_Wilson_coeffs;
  std::pair<NumericTensor<superMultiD,1>, NumericTensor<superMultiD,1> > lat_Wilson_coeffs;
  
  computeRIandLatticeWilsonCoefficients(RI_Wilson_coeffs,lat_Wilson_coeffs,wilson_coeffs,MSbar,data.NPR_sj,"");
  
  //Compute the Lellouch-Luscher factor
  superMultiD F_sj, delta_0_sj;
  computePhaseShiftAndLLfactor(delta_0_sj, F_sj, data.Epipi_sj, data.mpi_sj, data.mK_sj, data.ainv_sj, args);

  checkLatticeWilsonCoefficients(lat_Wilson_coeffs,data.M_lat_sj,data.ainv_sj,F_sj,args);


  NumericTensor<superMultiD,1> M_unrenorm_chiral_sj = computePhysicalUnrenormalizedMatrixElementsChiral(data.M_lat_sj,data.ainv_sj,F_sj,args);


  {
    std::cout << "\n\n---------------------------------------------------------------\n";
    std::cout << "Starting computation using standard method\n";
    
    //Compute the MSbar, infinite-volume matrix elements
    NumericTensor<superMultiD,1> M_MSbar_std_sj = computePhysicalMSbarMatrixElements(M_unrenorm_chiral_sj,data.NPR_sj,MSbar,args);

    computeMatrixElementRelations(M_MSbar_std_sj, "");
    
    //Compute A0
    superMultiD ReA0_sj, ImA0_sj;
    computeA0(ReA0_sj, ImA0_sj, M_MSbar_std_sj, wilson_coeffs, args, "");
    
    if(add_A0_sys)
      addA0syserr(ReA0_sj, ImA0_sj, reA0_sys_frac, imA0_sys_frac);

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
    NumericTensor<superMultiD,1> M_MSbar_std_sj = computePhysicalMSbarMatrixElementsUsingReA0expt(M_unrenorm_chiral_sj,data.NPR_sj,MSbar,lat_Wilson_coeffs.first,data.ReA0_expt_sj,i,args,nm.str());

    computeMatrixElementRelations(M_MSbar_std_sj, nm.str());
    
    //Compute A0
    superMultiD ReA0_sj, ImA0_sj;
    computeA0(ReA0_sj, ImA0_sj, M_MSbar_std_sj, wilson_coeffs, args, nm.str());
    
    if(add_A0_sys)
      addA0syserr(ReA0_sj, ImA0_sj, reA0_sys_frac, imA0_sys_frac);

    //Compute epsilon'
    computeEpsilon(ReA0_sj, ImA0_sj,
		   data.ReA2_lat_sj, data.ImA2_lat_sj,
		   delta_0_sj, data.delta_2_sj,
		   data.ReA0_expt_sj, data.ReA2_expt_sj,
		   data.omega_expt_sj, data.mod_eps_sj,
		   args, nm.str());
  }

  {
    std::cout << "\n\n---------------------------------------------------------------\n";
    std::cout << "Starting computation using Re(A0) to optimize statistical error" << std::endl;

    std::string nm = "opt_wReA0expt_";
        
    //Compute the MSbar, infinite-volume matrix elements
    NumericTensor<superMultiD,1> M_MSbar_std_sj = computePhysicalMSbarMatrixElementsUsingReA0exptOptimal(M_unrenorm_chiral_sj,data.NPR_sj,MSbar,
													 lat_Wilson_coeffs.first,lat_Wilson_coeffs.second,
													 data.ReA0_expt_sj,args,nm);

    computeMatrixElementRelations(M_MSbar_std_sj, nm);
    
    //Compute A0
    superMultiD ReA0_sj, ImA0_sj;
    computeA0(ReA0_sj, ImA0_sj, M_MSbar_std_sj, wilson_coeffs, args, nm);
    
    if(add_A0_sys)
      addA0syserr(ReA0_sj, ImA0_sj, reA0_sys_frac, imA0_sys_frac);


    //Compute epsilon'
    computeEpsilon(ReA0_sj, ImA0_sj,
		   data.ReA2_lat_sj, data.ImA2_lat_sj,
		   delta_0_sj, data.delta_2_sj,
		   data.ReA0_expt_sj, data.ReA2_expt_sj,
		   data.omega_expt_sj, data.mod_eps_sj,
		   args, nm);
  }
  
  std::cout << "Done\n";
  return 0;
}
