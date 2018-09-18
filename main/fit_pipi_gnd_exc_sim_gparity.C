#include<pipi_common/simfit_generic.h>

#include <fit_pipi_sigma_sim_gparity/read_data.h>
#include <fit_pipi_sigma_sim_gparity/resampled_correlator.h>
using namespace CPSfit;

#include<fit_pipi_gnd_exc_sim_gparity/args.h>
#include<fit_pipi_gnd_exc_sim_gparity/cmdline.h>
#include<fit_pipi_gnd_exc_sim_gparity/read_data.h>
#include<fit_pipi_gnd_exc_sim_gparity/resampled_data.h>
#include<fit_pipi_gnd_exc_sim_gparity/fit.h>



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
  doubleJackCorrelationFunction dj_data_gnd_gnd, dj_data_exc_exc, dj_data_gnd_exc;

  generateResampledDataOptions gopt; 
  gopt.load_combined_data = cmdline.load_combined_data;
  gopt.load_combined_data_file = cmdline.load_combined_data_file;
  gopt.save_combined_data = cmdline.save_combined_data;
  gopt.save_combined_data_file = cmdline.save_combined_data_file;

  generateResampledData(dj_data_gnd_gnd, dj_data_exc_exc, dj_data_gnd_exc,
			raw_bubble_gnd_gnd, raw_bubble_exc_exc, raw_bubble_gnd_exc,
			raw_data_gnd_gnd, raw_data_exc_exc, raw_data_gnd_exc,
			args.tsep_pipi, args.bin_size, 
			args.do_vacuum_subtraction, gopt);

  int nsample = dj_data_gnd_gnd.value(0).size();


  typedef std::unordered_map<std::string, std::string> SubFitFuncParameterMap;
  typedef std::unordered_map<std::string,size_t> ParamTagIdxMap;

  SubFitFuncParameterMap pipi2pt_gnd_gnd_oneexp_pmap, pipi2pt_exc_exc_oneexp_pmap, pipi2pt_gnd_exc_oneexp_pmap;
  ParamTagIdxMap param_map;

  typedef taggedValueContainer<double,std::string> Params;

  if(args.fitfunc == FitFuncType::FSimGenOneState){
    pipi2pt_gnd_gnd_oneexp_pmap = {  {"Asrc","Apipi_gnd"}, {"Asnk", "Apipi_gnd"}, {"E", "Epipi"}, {"Csys", "Cpipi_gnd_gnd"} };
    pipi2pt_exc_exc_oneexp_pmap = {  {"Asrc","Apipi_exc"}, {"Asnk", "Apipi_exc"}, {"E", "Epipi"}, {"Csys", "Cpipi_exc_exc"} };
    pipi2pt_gnd_exc_oneexp_pmap = {  {"Asrc","Apipi_gnd"}, {"Asnk", "Apipi_exc"}, {"E", "Epipi"}, {"Csys", "Cpipi_gnd_exc"} };
    
    param_map = { {"Apipi_gnd", 0}, {"Apipi_exc", 1},  {"Epipi",2}, {"Cpipi_gnd_gnd", 3}, {"Cpipi_exc_exc", 4},  {"Cpipi_gnd_exc", 5} };
  }else if(args.fitfunc == FitFuncType::FSimGenTwoState){
    pipi2pt_gnd_gnd_oneexp_pmap = {  {"Asrc0","Apipi_gnd_0"}, {"Asnk0", "Apipi_gnd_0"}, {"E0", "Epipi"},
				     {"Asrc1","Apipi_gnd_1"}, {"Asnk1", "Apipi_gnd_1"}, {"E1", "Eexc"},
				     {"Csys", "Cpipi_gnd_gnd"} };

    pipi2pt_exc_exc_oneexp_pmap = {  {"Asrc0","Apipi_exc_0"}, {"Asnk0", "Apipi_exc_0"}, {"E0", "Epipi"},
				     {"Asrc1","Apipi_exc_1"}, {"Asnk1", "Apipi_exc_1"}, {"E1", "Eexc"},
				     {"Csys", "Cpipi_exc_exc"} };

    pipi2pt_gnd_exc_oneexp_pmap = {  {"Asrc0","Apipi_gnd_0"}, {"Asnk0", "Apipi_exc_0"}, {"E0", "Epipi"},
				     {"Asrc1","Apipi_gnd_1"}, {"Asnk1", "Apipi_exc_1"}, {"E1", "Eexc"},
				     {"Csys", "Cpipi_gnd_exc"} };
    
    param_map = { {"Apipi_gnd_0", 0}, {"Apipi_exc_0", 1},  {"Epipi",2}, 
		  {"Apipi_gnd_1", 3}, {"Apipi_exc_1", 4},  {"Eexc",5},
		  {"Cpipi_gnd_gnd", 6}, {"Cpipi_exc_exc", 7},  {"Cpipi_gnd_exc", 8} };
  }else{
    assert(0);
  }

  std::unordered_map< SubFitFuncParameterMap const*, std::string > type_names = { {&pipi2pt_gnd_gnd_oneexp_pmap, "PiPi 2pt gnd/gnd"},
										  {&pipi2pt_exc_exc_oneexp_pmap, "PiPi 2pt exc/exc"},
										  {&pipi2pt_gnd_exc_oneexp_pmap, "PiPi 2pt gnd/exc"}  }; //just used for printing input data

  std::cout << "Selecting data in fit range" << std::endl;
  correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> corr_comb_dj;

  for(int t=args.t_min;t<=args.t_max;t++)
    corr_comb_dj.push_back( SimFitCoordGen(t, &pipi2pt_gnd_gnd_oneexp_pmap, 2*args.tsep_pipi), dj_data_gnd_gnd.value(t) );
  
  for(int t=args.t_min;t<=args.t_max;t++)
    corr_comb_dj.push_back( SimFitCoordGen(t, &pipi2pt_exc_exc_oneexp_pmap, 2*args.tsep_pipi), dj_data_exc_exc.value(t) );

  for(int t=args.t_min;t<=args.t_max;t++)
    corr_comb_dj.push_back( SimFitCoordGen(t, &pipi2pt_gnd_exc_oneexp_pmap, 2*args.tsep_pipi), dj_data_gnd_exc.value(t) );


  correlationFunction<SimFitCoordGen,  jackknifeDistributionD> corr_comb_j(corr_comb_dj.size(), 
    [&](const int i){
	return correlationFunction<SimFitCoordGen,  jackknifeDistributionD>::ElementType( corr_comb_dj.coord(i), corr_comb_dj.value(i).toJackknife() );
    });

  std::cout << "Data in fit:" << std::endl;
  for(int i=0;i<corr_comb_j.size();i++){
    std::cout << type_names[corr_comb_j.coord(i).param_map] << " " << corr_comb_j.coord(i).t << " " << corr_comb_j.value(i) << std::endl;
  }

  std::cout << "Setting up guess" << std::endl;

  Params guess(&param_map);
  if(args.fitfunc == FitFuncType::FSimGenOneState){
    guess("Apipi_gnd") = 0.2;
    guess("Apipi_exc") = 0.2;
    guess("Epipi") = 0.4;
    guess("Cpipi_gnd_gnd") = 0;
    guess("Cpipi_exc_exc") = 0;
  }else if(args.fitfunc == FitFuncType::FSimGenTwoState){
    guess("Apipi_gnd_0") = 0.2;
    guess("Apipi_exc_0") = 0.2;
    guess("Apipi_gnd_1") = 0.2;
    guess("Apipi_exc_1") = 0.2;
    guess("Epipi") = 0.4;
    guess("Eexc") = 0.7;
    guess("Cpipi_gnd_gnd") = 0;
    guess("Cpipi_exc_exc") = 0;
    guess("Cpipi_gnd_exc") = 0;
  }

  if(cmdline.load_guess){
    std::map<std::string,double> gmap;
    parse(gmap, cmdline.guess_file);
    for(auto it = param_map.begin(); it != param_map.end(); it++){
      const auto &tag = it->first;
      auto git = gmap.find(tag);
      if(git == gmap.end()) error_exit(std::cout << "In guess file " << cmdline.guess_file << " could not find tag " << tag << std::endl);
      guess(tag) = git->second;
    }
    std::cout << "Loaded guess: " << guess << std::endl;
  }

  std::cout << "Performing fit" << std::endl;

  jackknifeDistribution<Params> params(nsample, guess);
  jackknifeDistributionD chisq(nsample), chisq_per_dof(nsample);


  fitOptions opt;
  opt.load_frozen_fit_params = cmdline.load_frozen_fit_params;
  opt.load_frozen_fit_params_file = cmdline.load_frozen_fit_params_file;

  fit(params, chisq, chisq_per_dof,
      corr_comb_j, corr_comb_dj, args.fitfunc, param_map,
      args.Lt, args.Ascale, args.Cscale, opt);

  std::cout << "Params:\n" << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");

  std::cout << "Done\n";
  return 0;
}

