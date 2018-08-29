#include <fstream>
#include <algorithm>
#include <sstream>
#include <boost/timer/timer.hpp>

#include <fit.h>
#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>
#include <parser.h>
#include <common.h>

#include <pipi_common/pipi_common.h>

using namespace CPSfit;

#include <fit_sigmasigma_gparity/read_data.h>

#include <fit_pipitosigma_gparity/read_data.h>
#include <fit_pipitosigma_gparity/resampled_correlator.h>
#include <fit_pipitosigma_gparity/fit.h>
#include <fit_pipitosigma_gparity/args.h>
#include <fit_pipitosigma_gparity/cmdline.h>


int main(const int argc, const char* argv[]){
  PiPiToSigmaArgs args;
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
  
  PiPiToSigmaCMDline cmdline(argc,argv,2);

  //Read sigma bubble as complex
  sigmaSelfContractionZ sigma_self_data_Z;
  readSigmaSelf(sigma_self_data_Z, args.data_dir, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);

  sigmaSelfContraction sigma_self_data = reIm(sigma_self_data_Z, 0); //copy real part

  //Read pipi bubble as complex
  bubbleDataZ pipi_self_data_Z;
  getA1projectedSourcePiPiBubble(pipi_self_data_Z, args.data_dir, args.traj_start, args.traj_inc, args.traj_lessthan, args.tsep_pipi, args.Lt);

  bubbleData pipi_self_data = reIm(pipi_self_data_Z, 0); //real part

  //Get pipi->sigma data
  figureData pipitosigma_data;
  readPiPiToSigma(pipitosigma_data, args.data_dir, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
  
  //Reconstruct disconnected and connected part
  figureData pipitosigma_disconn_data_ReZZ;
  reconstructPiPiToSigmaDisconnected(pipitosigma_disconn_data_ReZZ, pipi_self_data_Z, sigma_self_data_Z); // Re ( pipi_bubble * sigma_bubble )
 
  figureData pipitosigma_disconn_data_ReZReZ;
  reconstructPiPiToSigmaDisconnected(pipitosigma_disconn_data_ReZReZ, pipi_self_data, sigma_self_data); // Re ( pipi_bubble ) * Re ( sigma_bubble )

  figureData pipitosigma_conn_data;
  reconstructPiPiToSigmaConnected(pipitosigma_conn_data, pipitosigma_data, pipitosigma_disconn_data_ReZZ, args.tstep_src);

  //(Very slightly) better statistics if we use the Re ( pipi_bubble ) * Re ( sigma_bubble ) for the disconnected part, taking advantage of the fact that the bubbles are real under the ensemble avg
  figureData &pipitosigma_disconn_data = pipitosigma_disconn_data_ReZReZ;

  //The code computes the disconnected component for all tsrc, but this option can be used to constrain the number of source timeslices to observe the effect
  if(cmdline.force_disconn_tstep_src){
    for(int tsrc=0;tsrc<args.Lt;tsrc++){
      if(tsrc % cmdline.disconn_tstep_src != 0)
	for(int t=0;t<args.Lt;t++) pipitosigma_disconn_data(tsrc, t).zero();
    }
  }

  //Source avg connected and disconnected parts and sum the contributions
  rawCorrelationFunction correlator_raw_conn = sourceAverage(pipitosigma_conn_data);
  rawCorrelationFunction correlator_raw_disconn = sourceAverage(pipitosigma_disconn_data);

  rawCorrelationFunction correlator_raw = correlator_raw_conn;
  for(int t=0;t<args.Lt;t++) correlator_raw.value(t) = correlator_raw.value(t) + correlator_raw_disconn.value(t);

  std::cout << "Raw data connected/disconnected parts:\n";
  for(int t=0;t<args.Lt;t++) std::cout << t << " " << correlator_raw_conn.value(t) << " " << correlator_raw_disconn.value(t) << std::endl;


  //Resample data
  doubleJackCorrelationFunction correlator_dj(args.Lt,
  					      [&](const int t){ 
  						return typename doubleJackCorrelationFunction::ElementType(t, doubleJackknifeDistributionD(correlator_raw.value(t).bin(args.bin_size)));
  					      }
  					      );

  const int nsample = correlator_dj.value(0).size();

  if(args.do_vacuum_subtraction){
    doubleJackCorrelationFunction vac_sub_dj = computePiPiToSigmaVacSub(sigma_self_data, pipi_self_data, args.bin_size);
    correlator_dj = correlator_dj - vac_sub_dj;
  }
  
  //Convert to single-jackknife
  jackknifeCorrelationFunction correlator_j(correlator_dj.size(),
  					    [&correlator_dj](const int i){
  					      return typename jackknifeCorrelationFunction::ElementType(double(correlator_dj.coord(i)), correlator_dj.value(i).toJackknife());
  					    });

  correlator_j = foldPiPiToSigma(correlator_j, args.Lt, args.tsep_pipi);
  correlator_dj = foldPiPiToSigma(correlator_dj, args.Lt, args.tsep_pipi);

  std::cout << "Resampled data:\n";
  for(int t=0;t<args.Lt;t++) std::cout << t << " " << correlator_j.value(t) << std::endl;

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
  PiPiToSigmaFitArgs fargs;
  args.transfer(fargs);
  cmdline.transfer(fargs);
  std::pair<jackknifeDistributionD,jackknifeDistributionD> Epipi_and_const = fit_pipi_to_sigma(correlator_j_inrange, correlator_dj_inrange, fargs);
  
  std::cout << "Done\n";
  return 0;
}
