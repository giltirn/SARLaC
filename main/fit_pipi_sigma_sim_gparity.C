#include <pipi_common/pipi_common.h>

using namespace CPSfit;

#include <fit_pipi_gparity/args.h>
#include <fit_pipi_gparity/cmdline.h>

#include <fit_pipi_sigma_sim_gparity/fitfunc.h>
#include <fit_pipi_sigma_sim_gparity/fit.h>
#include <fit_pipi_sigma_sim_gparity/args.h>
#include <fit_pipi_sigma_sim_gparity/cmdline.h>
#include <fit_pipi_sigma_sim_gparity/read_data.h>


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

  rawCorrelationFunction pipi_raw, pipi_to_sigma_raw, sigma2pt_raw;
  bubbleDataAllMomenta pipi_self_data;
  sigmaSelfContraction sigma_self_data;
  
  if(cmdline.load_checkpoint){
    readCheckpoint(pipi_raw,pipi_to_sigma_raw,sigma2pt_raw,pipi_self_data,sigma_self_data,cmdline.load_checkpoint_file);
  }else{
    readData(pipi_raw,pipi_to_sigma_raw,sigma2pt_raw,pipi_self_data,sigma_self_data,
	     args.data_dir, 
	     args.pipi2pt_figure_file_fmt, args.sigma2pt_file_fmt, args.pipitosigma_file_fmt,
	     args.pipi_bubble_file_fmt, args.sigma_bubble_file_fmt,
	     args.tsep_pipi, args.tstep_pipi2pt, args.tstep_pipitosigma,
	     args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, !cmdline.use_pipitosigma_disconn_complex_prod);
  }

  if(cmdline.save_checkpoint){
    writeCheckpoint(cmdline.save_checkpoint_file,pipi_raw,pipi_to_sigma_raw,sigma2pt_raw,pipi_self_data,sigma_self_data);
  }

  //Get double-jack data
  auto pipi_to_sigma_j = binResample<jackknifeCorrelationFunction>(pipi_to_sigma_raw, args.bin_size);
  auto pipi_to_sigma_dj = binResample<doubleJackCorrelationFunction>(pipi_to_sigma_raw, args.bin_size);

  auto sigma2pt_j = binResample<jackknifeCorrelationFunction>(sigma2pt_raw, args.bin_size);
  auto sigma2pt_dj = binResample<doubleJackCorrelationFunction>(sigma2pt_raw, args.bin_size);

  auto pipi_j = binResample<jackknifeCorrelationFunction>(pipi_raw, args.bin_size);
  auto pipi_dj = binResample<doubleJackCorrelationFunction>(pipi_raw, args.bin_size);
  
  //Compute vacuum subtractions
  if(args.do_vacuum_subtraction){
    { //Pipi->sigma
      PiPiProjectorA1Basis111 proj_pipi;
      bubbleData pipi_self_proj = projectSourcePiPiBubble(pipi_self_data, proj_pipi);
      pipi_to_sigma_j = pipi_to_sigma_j - computePiPiToSigmaVacSub<jackknifeCorrelationFunction>(sigma_self_data, pipi_self_proj, args.bin_size);
      pipi_to_sigma_dj = pipi_to_sigma_dj - computePiPiToSigmaVacSub<doubleJackCorrelationFunction>(sigma_self_data, pipi_self_proj, args.bin_size);
    }
    { //sigma->sigma
      sigma2pt_j = sigma2pt_j - computeSigmaVacSub<jackknifeCorrelationFunction>(sigma_self_data, args.bin_size);
      sigma2pt_dj = sigma2pt_dj - computeSigmaVacSub<doubleJackCorrelationFunction>(sigma_self_data, args.bin_size);
    }
    { //Pipi->pipi
      pipi_j = pipi_j - computePiPi2ptVacSub<jackknifeCorrelationFunction>(pipi_self_data, args.bin_size, args.tsep_pipi);
      pipi_dj = pipi_dj - computePiPi2ptVacSub<doubleJackCorrelationFunction>(pipi_self_data, args.bin_size, args.tsep_pipi);
    }
  }
  
  //Fold data
  pipi_to_sigma_j = fold(pipi_to_sigma_j, args.tsep_pipi);
  pipi_to_sigma_dj = fold(pipi_to_sigma_dj, args.tsep_pipi);

  sigma2pt_j = fold(sigma2pt_j, 0);
  sigma2pt_dj = fold(sigma2pt_dj, 0);

  pipi_j = fold(pipi_j, 2*args.tsep_pipi);
  pipi_dj = fold(pipi_dj, 2*args.tsep_pipi);
  
  //Build the combined data set
  simFitCorrFuncJ corr_comb_j;
  simFitCorrFuncDJ corr_comb_dj;
  jackknifeCorrelationFunction const* dsets_j[3] = { &pipi_j, &pipi_to_sigma_j, &sigma2pt_j };
  doubleJackCorrelationFunction const* dsets_dj[3] = { &pipi_dj, &pipi_to_sigma_dj, &sigma2pt_dj };

  SimFitType dtype[3] = { SimFitType::PiPi2pt, SimFitType::PiPiToSigma, SimFitType::Sigma2pt };
  bool dincl[3] = { cmdline.include_pipi_2pt, cmdline.include_pipi_to_sigma, cmdline.include_sigma_2pt };

  for(int d=0;d<3;d++)
    if(dincl[d])
      for(int t=args.t_min; t<=args.t_max; t++){ 
	corr_comb_j.push_back(SimFitCoord(dtype[d],t), dsets_j[d]->value(t));
	corr_comb_dj.push_back(SimFitCoord(dtype[d],t), dsets_dj[d]->value(t));
      }

  std::cout << "Data in fit range:\n";
  for(int i=0;i<corr_comb_j.size();i++){
    std::cout << corr_comb_j.coord(i).type << " " << corr_comb_j.coord(i).t << " " << corr_comb_j.value(i) << std::endl;
  }

  //Run the fit
  SimFitArgs fargs;
  args.transfer(fargs);
  cmdline.transfer(fargs);

  fit(corr_comb_j, corr_comb_dj, fargs);

  std::cout << "Done\n";
  return 0;
}
