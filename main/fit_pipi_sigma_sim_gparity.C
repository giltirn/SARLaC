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

using namespace CPSfit;

#include <fit_pipi_gparity/data_containers.h>
#include <fit_pipi_gparity/threemomentum.h>
#include <fit_pipi_gparity/mom_data_containers.h>
#include <fit_pipi_gparity/mom_project.h>
#include <fit_pipi_gparity/fitfunc.h>
#include <fit_pipi_gparity/raw_correlator.h>
#include <fit_pipi_gparity/read_data.h>
#include <fit_pipi_gparity/args.h>
#include <fit_pipi_gparity/cmdline.h>
#include <fit_pipi_gparity/raw_data.h>
#include <fit_pipi_gparity/resampled_data.h>

#include <fit_sigmasigma_gparity/read_data.h>
#include <fit_sigmasigma_gparity/resampled_correlator.h>

#include <fit_pipitosigma_gparity/read_data.h>
#include <fit_pipitosigma_gparity/resampled_correlator.h>

#include <fit_pipi_sigma_sim_gparity/fitfunc.h>
#include <fit_pipi_sigma_sim_gparity/fit.h>
#include <fit_pipi_sigma_sim_gparity/args.h>
#include <fit_pipi_sigma_sim_gparity/cmdline.h>
#include <fit_pipi_sigma_sim_gparity/read_data.h>
#include <fit_pipi_sigma_sim_gparity/resampled_correlator.h>


int main(const int argc, const char* argv[]){
  PiPiSigmaSimArgs args;
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

  PiPiSigmaSimCMDline cmdline(argc,argv,2);

  std::vector<threeMomentum> pion_mom = { {1,1,1}, {-1,-1,-1},
					  {-1,1,1}, {1,-1,-1},
					  {1,-1,1}, {-1,1,-1},
					  {1,1,-1}, {-1,-1,1} };

  rawCorrelationFunction pipi_raw, pipi_to_sigma_raw, sigma2pt_raw;
  bubbleDataAllMomenta pipi_self_data;
  sigmaSelfContraction sigma_self_data;
  
  if(cmdline.load_checkpoint){
    readCheckpoint(pipi_raw,pipi_to_sigma_raw,sigma2pt_raw,pipi_self_data,sigma_self_data,cmdline.load_checkpoint_file);
  }else{
    readData(pipi_raw,pipi_to_sigma_raw,sigma2pt_raw,pipi_self_data,sigma_self_data,pion_mom,args);
  }
 
  if(cmdline.save_checkpoint){
    writeCheckpoint(cmdline.save_checkpoint_file,pipi_raw,pipi_to_sigma_raw,sigma2pt_raw,pipi_self_data,sigma_self_data);
  }

  //Get double-jack data
  doubleJackCorrelationFunction pipi_to_sigma_dj = binDoubleJackResample(pipi_to_sigma_raw, args.bin_size);
  doubleJackCorrelationFunction sigma2pt_dj = binDoubleJackResample(sigma2pt_raw, args.bin_size);
  doubleJackCorrelationFunction pipi_dj = binDoubleJackResample(pipi_raw, args.bin_size);
  
  //Compute vacuum subtractions
  if(args.do_vacuum_subtraction){
    { //Pipi->sigma
      bubbleData pipi_self_proj = A1projectPiPiBubble(pipi_self_data, pion_mom, args.Lt, args.tsep_pipi);
      doubleJackCorrelationFunction vac_sub_dj = computePiPiToSigmaVacSub(sigma_self_data, pipi_self_proj, args.bin_size);
      pipi_to_sigma_dj = pipi_to_sigma_dj - vac_sub_dj;
    }
    { //sigma->sigma
      doubleJackCorrelationFunction vac_sub_dj = computeSigmaVacSub(sigma_self_data, args.bin_size);
      sigma2pt_dj = sigma2pt_dj - vac_sub_dj;
    }
    { //Pipi->pipi
      doubleJackCorrelationFunction vac_sub_dj = computePiPi2ptVacSub(pipi_self_data, args.bin_size, args.tsep_pipi, pion_mom);
      pipi_dj = pipi_dj - vac_sub_dj;
    }
  }
  

  //Fold data
  pipi_to_sigma_dj = foldPiPiToSigma(pipi_to_sigma_dj, args.Lt, args.tsep_pipi);
  sigma2pt_dj = foldSigma(sigma2pt_dj, args.Lt);
  pipi_dj = fold(pipi_dj, args.tsep_pipi);
  
  //Build the combined data set
  simFitCorrFuncDJ corr_comb_dj;
  doubleJackCorrelationFunction* dsets[3] = { &pipi_dj, &pipi_to_sigma_dj, &sigma2pt_dj };
  SimFitType dtype[3] = { PiPi2pt, PiPiToSigma, Sigma2pt };
  bool dincl[3] = { cmdline.include_pipi_2pt, cmdline.include_pipi_to_sigma, cmdline.include_sigma_2pt };

  for(int d=0;d<3;d++)
    if(dincl[d])
      for(int t=args.t_min; t<=args.t_max; t++) 
	corr_comb_dj.push_back( simFitCorrFuncDJ::ElementType(  SimFitCoord(dtype[d],t), dsets[d]->value(t) ) );
  
  simFitCorrFuncJ corr_comb_j(corr_comb_dj.size(), [&](const int i){ return simFitCorrFuncJ::ElementType(corr_comb_dj.coord(i), corr_comb_dj.value(i).toJackknife()); });

  //Run the fit
  SimFitArgs fargs;
  args.transfer(fargs);
  cmdline.transfer(fargs);

  fit(corr_comb_j, corr_comb_dj, fargs);

  std::cout << "Done\n";
  return 0;
}
