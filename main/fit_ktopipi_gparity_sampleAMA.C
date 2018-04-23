#include<array>
#include<vector>
#include<complex>
#include<iostream>
#include<fstream>

#include <boost/timer/timer.hpp>

#include<utils.h>
#include<distribution.h>
#include<common.h>
#include<tensors.h>
#include<data_series.h>
#include<fit.h>
#include<parser.h>
#include<serialize.h>
#include<plot.h>
#include<containers.h>

using namespace CPSfit;

#include <fit_pipi_gparity/threemomentum.h>
#include <fit_pipi_gparity/data_containers.h>
#include <fit_pipi_gparity/mom_project.h>
#include <fit_pipi_gparity/mom_data_containers.h>
#include <fit_pipi_gparity/read_data.h>

#include <fit_ktopipi_gparity/utils.h>
#include <fit_ktopipi_gparity/freeze.h>
#include <fit_ktopipi_gparity/cmdline.h>
#include <fit_ktopipi_gparity/args.h>
#include <fit_ktopipi_gparity/data_containers.h>
#include <fit_ktopipi_gparity/read_data.h>
#include <fit_ktopipi_gparity/compute_amplitude.h>
#include <fit_ktopipi_gparity/amplitude_data.h>
#include <fit_ktopipi_gparity/fitfunc.h>
#include <fit_ktopipi_gparity/plot.h>
#include <fit_ktopipi_gparity/fit.h>
#include <fit_ktopipi_gparity/scratch.h>
#include <fit_ktopipi_gparity/main.h>

#include <fit_simple_sampleAMA/sampleAMA_resample.h>
#include <fit_ktopipi_gparity_sampleAMA/cmdline.h>
#include <fit_ktopipi_gparity_sampleAMA/args.h>
#include <fit_ktopipi_gparity_sampleAMA/resample_correct.h>
#include <fit_ktopipi_gparity_sampleAMA/data_structs.h>
#include <fit_ktopipi_gparity_sampleAMA/resample_average_typedata.h>
#include <fit_ktopipi_gparity_sampleAMA/alpha_vac_sub.h>
#include <fit_ktopipi_gparity_sampleAMA/main.h>

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
  
  allInputs inputs(args,cmdline);

  //Prepare the data
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > A0_all_dj(10);
  printMem("Prior to getData");
    
  if(cmdline.checkpoint_and_exit){
    checkpointRawOnly(inputs);
    return 0;
  }

  getDataSampleAMA(A0_all_j, A0_all_dj, inputs);

  printMem("Prior to fitting");
  fitAndPlot(A0_all_j,A0_all_dj,inputs.args_S,inputs.cmdline_sloppy_S);
  
  std::cout << "Done" << std::endl;
  
  return 0;
}

