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
#include<parser.h>
#include<hdf5_serialize.h>

#include <fit_pipi_gparity/data_containers.h>
#include <fit_pipi_gparity/mom_data_containers.h>
#include <fit_pipi_gparity/read_data.h>

#include <fit_ktopipi_gparity/cmdline.h>
#include <fit_ktopipi_gparity/args.h>
#include <fit_ktopipi_gparity/data_containers.h>
#include <fit_ktopipi_gparity/read_data.h>
#include <fit_ktopipi_gparity/compute_amplitude.h>
#include <fit_ktopipi_gparity/fitfunc.h>
#include <fit_ktopipi_gparity/plot.h>
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

  CMDline cmdline(argc,argv,2);
  
  typedef FitKtoPiPi FitFunc;
  std::vector<typename FitFunc::Params> guess(10);

  //Prepare the data
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > A0_all_dj(10);
  getData(A0_all_j, A0_all_dj,args,cmdline);

  //Compute fit weights from double jackknife
  std::vector<NumericVector<jackknifeDistributionD> > sigma_all_j(10);
  getSigma(sigma_all_j, A0_all_dj);  
  
  //Pull out data in fit range
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_fit(10);
  std::vector<NumericVector<jackknifeDistributionD> > sigma_fit(10);
  getFitData(A0_fit,sigma_fit,A0_all_j,sigma_all_j,args);
  
  std::cout << "Including " << A0_fit[0].size() << " data points in fit\n";
  for(int q=0;q<10;q++){
    std::cout << "For Q" << q+1 << std::endl;
    for(int i=0;i<A0_fit[q].size();i++)
      std::cout << A0_fit[q].coord(i) << " : " << A0_fit[q].value(i) << "  with sigma=" << sigma_fit[q](i) << std::endl;
  }

  //Perform the fit
  typedef typename composeFitPolicy<amplitudeDataCoord, FitFunc, frozenFitFuncPolicy, uncorrelatedFitPolicy>::type FitPolicies;
  FitFunc fitfunc;
  std::vector<jackknifeDistribution<FitFunc::Params> > fit_params(10);
  
  for(int q=0;q<10;q++){
    std::cout << "Starting fit for Q=" << q+1 << std::endl;
    fitter<FitPolicies> fitter;
    fitter.importFitFunc(fitfunc);
    fitter.importCostFunctionParameters(sigma_fit[q]);  

    readFrozenParams(fitter, q+1, cmdline, nsample);
  
    jackknifeDistribution<FitFunc::Params> &params = fit_params[q];
    params = jackknifeDistribution<FitFunc::Params>(nsample, guess[q]);
    jackknifeDistributionD chisq;
    jackknifeDistributionD chisq_per_dof;
    fitter.fit(params, chisq, chisq_per_dof, A0_fit[q]);

    distributionPrint<jackknifeDistribution<FitFunc::Params> >::printer(new ktopipiParamsPrinter<FitFunc>);

    std::cout << "Q" << q << " Params: " << params << std::endl;
    std::cout << "Q" << q << " Chisq: " << chisq << std::endl;
    std::cout << "Q" << q << " Chisq/dof: " << chisq_per_dof << std::endl;
  }
  extractMdata<FitFunc> extractor(fit_params);
  plotErrorWeightedData(A0_all_j,extractor,args);

#ifdef HAVE_HDF5
  {
    std::vector<jackknifeDistributionD> M(10);
    for(int q=0;q<10;q++) M[q] = jackknifeDistributionD(nsample, [&](const int s){ return extractor.getMfit(q,s); } );
    HDF5writer writer("result_M.hdf5");
    write(writer, M , "value");
  }
#endif
  
  return 0;
}

