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
#include <fit_ktopipi_gparity/fit.h>
#include <fit_ktopipi_gparity/main.h>

int main(const int argc, const char* argv[]){
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  const std::string arg_file = argv[1];
  parse(args, arg_file);
  
  assert( (args.traj_lessthan - args.traj_start) % args.traj_inc == 0 );  
  const int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc;

  CMDline cmdline(argc,argv,2);
  
  //Prepare the data
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > A0_all_dj(10);
  getData(A0_all_j, A0_all_dj,args,cmdline);

  //Extract the data we are going to fit
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_fit_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > A0_fit_dj(10);
  getFitData(A0_fit_j,A0_fit_dj,A0_all_j,A0_all_dj,args);
  
  std::cout << "Including " << A0_fit_j[0].size() << " data points in fit\n";
  for(int q=0;q<10;q++){
    std::cout << "For Q" << q+1 << std::endl;
    for(int i=0;i<A0_fit_j[q].size();i++)
      std::cout << A0_fit_j[q].coord(i) << " : " << A0_fit_j[q].value(i) << std::endl;
  }

  typedef FitKtoPiPi FitFunc;
  std::vector<jackknifeDistribution<typename FitFunc::Params> > fit_params = fit<FitFunc>(A0_fit_j,A0_fit_dj,args,cmdline);

  extractMdata<FitFunc> extractor(fit_params);
  plotErrorWeightedData(A0_all_j,extractor,args);

#ifdef HAVE_HDF5
  writeParamsStandard(fit_params, "params.hdf5");
#endif

  std::cout << "Done" << std::endl;
  
  return 0;
}

