#include<array>
#include<vector>
#include<complex>
#include<iostream>
#include<fstream>

#include <boost/timer/timer.hpp>

#include<utils.h>
#include<distribution.h>
#include<common_defs.h>
#include<numeric_tensors.h>
#include<correlationfunction.h>
#include<fit_wrapper.h>
#include <parser.h>

#include <fit_pipi_gparity/data_containers.h>
#include <fit_pipi_gparity/mom_data_containers.h>
#include <fit_pipi_gparity/read_data.h>

#include <fit_ktopipi_gparity/args.h>
#include <fit_ktopipi_gparity/data_containers.h>
#include <fit_ktopipi_gparity/read_data.h>
#include <fit_ktopipi_gparity/compute_amplitude.h>
#include <fit_ktopipi_gparity/fitfunc.h>
#include <fit_ktopipi_gparity/main.h>


int main(const int argc, const char* argv[]){
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  {
    const std::string arg_file = argv[1];
    std::ifstream arg_f(arg_file.c_str());
    assert(arg_f.good());
    arg_f >> args;
    assert(!arg_f.bad() && !arg_f.fail());
    arg_f.close();
    std::cout << "Read arguments: \n" << args << std::endl;
  }
  assert( (args.traj_lessthan - args.traj_start) % args.traj_inc == 0 );  
  const int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc;
  
  typedef FitKtoPiPi FitFunc;
  std::vector<typename FitFunc::Params> guess(10);

  //Read the bubble data
  NumericTensor<rawDataDistributionD,1> bubble = readA2projectedBubble(args.traj_start,args.traj_inc,args.traj_lessthan,args.tsep_pipi,args.Lt,args.data_dir);
  NumericTensor<doubleJackknifeDistributionD,1> bubble_dj = bubble.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());

  //Read and prepare the amplitude data for fitting
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_fit(10);
  std::vector<NumericVector<jackknifeDistributionD> > sigma_fit(10);
  for(int tsep_k_pi_idx=0;tsep_k_pi_idx<args.tsep_k_pi.size();tsep_k_pi_idx++)
    getData(A0_fit,sigma_fit,bubble,bubble_dj,tsep_k_pi_idx,args);
  
  std::cout << "Including " << A0_fit[0].size() << " data points in fit\n";
  for(int q=0;q<10;q++){
    std::cout << "For Q" << q+1 << std::endl;
    for(int i=0;i<A0_fit[q].size();i++)
      std::cout << A0_fit[q].coord(i) << " : " << A0_fit[q].value(i) << "  with sigma=" << sigma_fit[q](i) << std::endl;
  }

  //Perform the fit
  typedef typename composeFitPolicy<amplitudeDataCoord, FitFunc, frozenFitFuncPolicy, uncorrelatedFitPolicy>::type FitPolicies;
  FitFunc fitfunc(args.AKscale, args.Apipiscale);

  for(int q=0;q<10;q++){  
    fitter<FitPolicies> fitter;
    fitter.importFitFunc(fitfunc);
    fitter.importCostFunctionParameters(sigma_fit[q]);  

    jackknifeDistribution<FitFunc::Params> freeze(nsample, FitFunc::Params(1,0.36,1,0.36,0));
    fitter.freeze({0,1,2,3}, freeze);
  
    jackknifeDistribution<FitFunc::Params> params(nsample, guess[q]);
    jackknifeDistributionD chisq;
    jackknifeDistributionD chisq_per_dof;
    fitter.fit(params, chisq, chisq_per_dof, A0_fit[q]);

    distributionPrint<decltype(params)>::printer(new ktopipiParamsPrinter<FitFunc>);

    std::cout << "Params: " << params << std::endl;
    std::cout << "Chisq: " << chisq << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  }

    
  
  return 0;
}

