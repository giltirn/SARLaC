#include<pipi_common/pipi_common.h>

#include <fit_pipi_sigma_sim_gparity/read_data.h>
using namespace CPSfit;

#include<fit_pipi_gnd_exc_sim_gparity/args.h>
#include<fit_pipi_gnd_exc_sim_gparity/cmdline.h>
#include<fit_pipi_gnd_exc_sim_gparity/read_data.h>
#include<fit_pipi_gnd_exc_sim_gparity/resampled_data.h>
#include<fit_pipi_gnd_exc_sim_gparity/fit.h>
#include<fit_pipi_gnd_exc_sim_gparity/fitfunc.h>


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
    parse(args,arg_file);
    std::cout << "Read arguments: \n" << args << std::endl;
  }

  CMDline cmdline(argc,argv,2);

  std::map< std::pair<std::string,std::string>, SubFitFuncParameterMap > subfit_pmaps;
  ParamTagIdxMap param_map;
  Params guess;

  //Setup the fit parameter maps
  setupParameterMaps(subfit_pmaps, param_map, guess, args.fitfunc, args.Ascale);

  //Write guess template and exit if requested
  if(cmdline.save_guess_template){ saveGuess(guess, "guess_template.dat"); exit(0); }

  //Get raw data
  bubbleDataAllMomenta raw_bubble_gnd_gnd, raw_bubble_exc_exc, raw_bubble_gnd_exc;
  rawCorrelationFunction raw_data_gnd_gnd, raw_data_exc_exc, raw_data_gnd_exc;
  
  readRawDataOptions ropt;
  ropt.load_hdf5_data_checkpoint = cmdline.load_hdf5_data_checkpoint;
  ropt.load_hdf5_data_checkpoint_stub = cmdline.load_hdf5_data_checkpoint_stub;
  ropt.save_hdf5_data_checkpoint = cmdline.save_hdf5_data_checkpoint;
  ropt.save_hdf5_data_checkpoint_stub = cmdline.save_hdf5_data_checkpoint_stub;
  
  if(!cmdline.load_combined_data) getRawPiPiGndExcData(raw_bubble_gnd_gnd, raw_bubble_exc_exc, raw_bubble_gnd_exc,
					     raw_data_gnd_gnd, raw_data_exc_exc, raw_data_gnd_exc,
					     args.data_dir, args.figure_file_format, args.bubble_file_format,
					     args.Lt, args.tsep_pipi, args.tstep_pipi,
					     args.traj_start, args.traj_inc, args.traj_lessthan, ropt);
  

  //Get double-jack data
  jackknifeCorrelationFunction j_data_gnd_gnd, j_data_exc_exc, j_data_gnd_exc;
  doubleJackCorrelationFunction dj_data_gnd_gnd, dj_data_exc_exc, dj_data_gnd_exc;

  generateResampledDataOptions gopt; 
  gopt.load_combined_data = cmdline.load_combined_data;
  gopt.load_combined_data_file = cmdline.load_combined_data_file;
  gopt.save_combined_data = cmdline.save_combined_data;
  gopt.save_combined_data_file = cmdline.save_combined_data_file;

  generateResampledData(j_data_gnd_gnd, dj_data_gnd_gnd, 
			j_data_exc_exc, dj_data_exc_exc, 
			j_data_gnd_exc, dj_data_gnd_exc,
			raw_bubble_gnd_gnd, raw_bubble_exc_exc, raw_bubble_gnd_exc,
			raw_data_gnd_gnd, raw_data_exc_exc, raw_data_gnd_exc,
			args.tsep_pipi, args.bin_size, 
			args.do_vacuum_subtraction, gopt);

  int nsample = dj_data_gnd_gnd.value(0).size();

  //Put all the data in the same container with appropriate generalized coordinate
  std::cout << "Selecting data in fit range" << std::endl;
  correlationFunction<SimFitCoordGen,  jackknifeDistributionD> corr_comb_j;
  correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> corr_comb_dj;

  static const std::vector<std::string> op = { "gnd", "exc" };
  jackknifeCorrelationFunction const* dptrs_j[] = { &j_data_gnd_gnd, &j_data_gnd_exc, &j_data_exc_exc };
  doubleJackCorrelationFunction const* dptrs_dj[] = { &dj_data_gnd_gnd, &dj_data_gnd_exc, &dj_data_exc_exc };
  
  std::vector<std::string> vnm; //for printing

  int d=0;
  for(int i=0;i<2;i++){
    for(int j=i; j<2; j++){
      std::string nm = op[i] + "/" + op[j];
      for(int t=args.t_min;t<=args.t_max;t++){
	SimFitCoordGen coord(t, &subfit_pmaps[{op[i],op[j]}], 2*args.tsep_pipi);
	corr_comb_j.push_back(coord, dptrs_j[d]->value(t) );
	corr_comb_dj.push_back(coord, dptrs_dj[d]->value(t) );
	vnm.push_back(nm);
      }
      d++;
    }
  }

  std::cout << "Data in fit:" << std::endl;
  for(int i=0;i<corr_comb_j.size();i++){
    std::cout << vnm[i] << " " << corr_comb_j.coord(i).t << " " << corr_comb_j.value(i) << std::endl;
  }

  if(cmdline.load_guess) loadGuess(guess, cmdline.guess_file);

  std::cout << "Performing fit" << std::endl;

  jackknifeDistribution<Params> params(nsample, guess);
  jackknifeDistributionD chisq(nsample), chisq_per_dof(nsample);


  fitOptions opt;
  opt.load_frozen_fit_params = cmdline.load_frozen_fit_params;
  opt.load_frozen_fit_params_file = cmdline.load_frozen_fit_params_file;

  const bool correlated = true;
  fit(params, chisq, chisq_per_dof,
      corr_comb_j, corr_comb_dj, args.fitfunc, param_map,
      args.Lt, args.t_min, args.t_max, correlated, args.Ascale, args.Cscale, opt);

  std::cout << "Params:\n";
  {
    for(int p=0;p<params.sample(0).size();p++){
      std::string tag = params.sample(0).tag(p);
      jackknifeDistributionD tmp;
      standardIOhelper<jackknifeDistributionD, jackknifeDistribution<Params> >::extractStructEntry(tmp, params, p);
      std::cout << tag << " = " << tmp << std::endl;
    }
  }
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");

  std::cout << "Done\n";
  return 0;
}

