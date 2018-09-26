#include <utils.h>
#include <pipi_common/pipi_common.h>

using namespace CPSfit;

#include<fit_pipi_gnd_exc_sim_gparity/enums.h>
#include<fit_pipi_gnd_exc_sim_gparity/fit.h>

#include<fit_pipi_gnd_exc_sigma_sim_gparity/fitfunc.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/args.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/cmdline.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/raw_data.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/resampled_data.h>

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

  //Setup the fit func parameter maps
  std::map< std::pair<Operator,Operator>, SubFitFuncParameterMap > subfit_pmaps;
  ParamTagIdxMap param_map;
  Params guess;
  setupParameterMaps(subfit_pmaps, param_map, guess, args.fitfunc, args.Ascale);

  //Write guess template and exit if requested
  if(cmdline.save_guess_template){ saveGuess(guess, "guess_template.dat"); exit(0); }
  
  //Read data
  RawData raw_data;
  if(!cmdline.load_combined_data) raw_data.read(args.Lt, args.data_dir, args.traj_start, args.traj_inc, args.traj_lessthan,
						args.pipi_figure_file_format, args.pipi_bubble_file_format, args.tsep_pipi, args.tstep_pipi,
						args.pipi_to_sigma_file_format, args.tstep_pipi_to_sigma,
						args.sigma2pt_file_format, args.sigma_bubble_file_format);

  ResampledData data;
  if(cmdline.load_combined_data) data.loadCheckpoint(cmdline.load_combined_data_file);
  else data.generatedResampledData(raw_data, args.bin_size, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction);
  
  if(cmdline.save_combined_data) data.saveCheckpoint(cmdline.save_combined_data_file);

  int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc/args.bin_size;

  //Add resampled data to full data set with generalized coordinate set appropriately
  correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> corr_comb_dj;
  
  std::vector<std::string> vnm; //for printing
  static const std::vector<Operator> ops = {PiPiGnd, PiPiExc, Sigma};

  for(int i=0;i<3;i++){
    for(int j=i;j<3;j++){
      std::ostringstream nm; nm << ops[i] << " " << ops[j];
      for(int t=args.t_min;t<=args.t_max;t++){
	corr_comb_dj.push_back( 
			       SimFitCoordGen(t, &subfit_pmaps[{ops[i],ops[j]}] , foldOffsetMultiplier(ops[i],ops[j])*args.tsep_pipi),
			       data.correlator(ops[i],ops[j]).value(t) 
				);
	vnm.push_back(nm.str());
      }
    }
  }
  
  correlationFunction<SimFitCoordGen,  jackknifeDistributionD> corr_comb_j(corr_comb_dj.size(), 
		[&](const int i){
		     return correlationFunction<SimFitCoordGen,  jackknifeDistributionD>::ElementType( corr_comb_dj.coord(i), corr_comb_dj.value(i).toJackknife() );
		});

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

  fit(params, chisq, chisq_per_dof,
      corr_comb_j, corr_comb_dj, args.fitfunc, param_map,
      args.Lt, args.Ascale, args.Cscale, opt);

  std::cout << "Params:\n";
  {
    for(int p=0;p<params.sample(0).size();p++){
      std::string tag = params.sample(0).tag(p);
      jackknifeDistributionD tmp;
      standardIOhelper<jackknifeDistributionD, jackknifeDistribution<Params> >::extractStructEntry(tmp, params, p);
      std::cout << tag << " = " << tmp << std::endl;
    }
  }
  std::cout << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");

  std::cout << "Done\n";
  return 0;
}


