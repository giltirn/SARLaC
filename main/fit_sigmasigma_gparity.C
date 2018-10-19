
#include<pipi_common/pipi_common.h>

using namespace CPSfit;

#include <fit_sigmasigma_gparity/fit.h>
#include <fit_sigmasigma_gparity/args.h>
#include <fit_sigmasigma_gparity/cmdline.h>

int main(const int argc, const char* argv[]){
  SigmaArgs args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  {
    const std::string arg_file = argv[1];
    parse(args,arg_file);
    std::cout << "Read arguments: \n" << args << std::endl;
  }

  SigmaCMDline cmdline(argc,argv,2);
 
  figureData raw_data;
  readSigmaSigma(raw_data, args.data_dir, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);

  sigmaSelfContraction raw_self;
  readSigmaSelf(raw_self, args.data_dir, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);

  rawCorrelationFunction correlator_raw = sourceAverage(raw_data);

  auto correlator_j = binResample<jackknifeCorrelationFunction>(correlator_raw, args.bin_size);
  auto correlator_dj = binResample<doubleJackCorrelationFunction>(correlator_raw, args.bin_size);

  const int nsample = correlator_dj.value(0).size();

  if(args.do_vacuum_subtraction){
    correlator_j = correlator_j - computeSigmaVacSub<jackknifeCorrelationFunction>(raw_self, args.bin_size);
    correlator_dj = correlator_dj - computeSigmaVacSub<doubleJackCorrelationFunction>(raw_self, args.bin_size);
  }
  
  correlator_j = fold(correlator_j, 0);
  correlator_dj = fold(correlator_dj, 0);

  //Filter out the data that is to be fitted
  filterXrange<double> trange(args.t_min,args.t_max);
  
  doubleJackCorrelationFunction correlator_dj_inrange;
  jackknifeCorrelationFunction correlator_j_inrange;
  for(int d=0;d<correlator_dj.size();d++)
    if(trange.accept(correlator_dj.coord(d),correlator_dj.value(d) )){
      correlator_dj_inrange.push_back(correlator_dj[d]);
      correlator_j_inrange.push_back(correlator_j[d]);
    }
  
  //Perform the fit
  SigmaFitArgs fargs;
  args.transfer(fargs);
  cmdline.transfer(fargs);
  std::pair<jackknifeDistributionD,jackknifeDistributionD> Epipi_and_const = fit_sigma(correlator_j_inrange,correlator_dj_inrange,fargs);
  
  std::cout << "Done\n";
  return 0;
}

