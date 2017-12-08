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

#include <fit_ktopipi_gparity/freeze.h>
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
  
  assert( (args.traj_lessthan - args.traj_start) % args.traj_inc == 0 );  
  const int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc;

  CMDline cmdline(argc,argv,2);
  
  //Prepare the data
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > A0_all_dj(10);
  printMem("Prior to getData");
  getData(A0_all_j, A0_all_dj,args,cmdline);

  printMem("Prior to fitting");
  fitAndPlot(A0_all_j,A0_all_dj,args,cmdline);

  std::cout << "Done" << std::endl;
  
  return 0;
}

