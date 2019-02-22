#include <utils.h>
#include <random.h>
#include <pipi_common/pipi_common.h>
#include <pipi_common/analyze_chisq.h>
using namespace CPSfit;

#include<fit_pipi_gnd_exc_sim_gparity/enums.h>
#include<fit_pipi_gnd_exc_sim_gparity/fit.h>

#include<fit_pipi_gnd_exc_sigma_sim_gparity/args.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/fitfunc.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/cmdline.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/raw_data.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/resampled_data.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/fit.h>

int main(const int argc, const char* argv[]){
  RNG.initialize(1234);

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

  CMDline cmdline(argc,argv,args.minimizer,2);

  //Setup the fit func parameter maps
  std::map< std::pair<Operator,Operator>, SubFitFuncParameterMap > subfit_pmaps;
  ParamTagIdxMap param_map;
  Params guess;
  setupParameterMaps(subfit_pmaps, param_map, guess, args.operators, args.fitfunc, args.Ascale, args.nstate);

  //Write guess template and exit if requested
  if(cmdline.save_guess_template){ saveGuess(guess, "guess_template.dat"); exit(0); }
  
  const std::vector<Operator> &ops = args.operators;

  //Read data
  RawData raw_data;
  if(!cmdline.load_combined_data){
    if(cmdline.load_raw_data){
      std::cout << "Reading raw data from " << cmdline.load_raw_data_file << std::endl;
      HDF5reader rd(cmdline.load_raw_data_file);  raw_data.read(rd, "raw_data");
      for(int i=0;i<ops.size();i++)
	for(int j=i;j<ops.size();j++)
	  assert(raw_data.haveData(ops[i],ops[j]));
    }else{
      raw_data.read(args.Lt, args.data_dir, args.traj_start, args.traj_inc, args.traj_lessthan,
		    args.pipi_figure_file_format, args.pipi_bubble_file_format, args.tsep_pipi, args.tstep_pipi,
		    args.pipi_to_sigma_file_format, args.tstep_pipi_to_sigma,
		    args.sigma2pt_file_format, args.sigma_bubble_file_format,
		    ops);
    }
    if(cmdline.save_raw_data){
      std::cout << "Saving raw data to " << cmdline.save_raw_data_file << std::endl;
      HDF5writer wr(cmdline.save_raw_data_file);  raw_data.write(wr, "raw_data");
    }

    if(cmdline.remove_samples_in_range)
      raw_data.removeSamplesInRange(cmdline.remove_samples_in_range_start, cmdline.remove_samples_in_range_lessthan);
  }

  ResampledData<jackknifeCorrelationFunction> data_j;
  ResampledData<doubleJackCorrelationFunction> data_dj;
  if(cmdline.load_combined_data){
    loadCheckpoint(data_j, data_dj, cmdline.load_combined_data_file);
    for(int i=0;i<ops.size();i++)
      for(int j=i;j<ops.size();j++)
	if(!data_j.haveData(ops[i],ops[j]))
	  error_exit(std::cout << "Loaded checkpoint does not contain data for (" << ops[i] << ", " << ops[j] << ")\n");
  }else{
    data_j.generatedResampledData(raw_data, args.bin_size, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction, args.timeslice_avg_vac_sub);
    data_dj.generatedResampledData(raw_data, args.bin_size, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction, args.timeslice_avg_vac_sub);
  }  

  if(cmdline.save_combined_data) saveCheckpoint(data_j, data_dj, cmdline.save_combined_data_file);

  const int nsample = data_j.getNsample();
  std::cout << "Number of binned samples is " << nsample << std::endl;

  //Add resampled data to full data set with generalized coordinate set appropriately
  correlationFunction<SimFitCoordGen,  jackknifeDistributionD> corr_comb_j;
  correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> corr_comb_dj;
  
  std::vector<std::string> vnm; //for printing

  std::map<std::unordered_map<std::string, std::string> const*, std::string> pmap_descr;
  for(int i=0;i<ops.size();i++){
    for(int j=i;j<ops.size();j++){
      std::ostringstream nm; nm << ops[i] << " " << ops[j];
      pmap_descr[&subfit_pmaps[{ops[i],ops[j]}]] = nm.str();

      for(int t=args.t_min;t<=args.t_max;t++){
	SimFitCoordGen coord(t, &subfit_pmaps[{ops[i],ops[j]}] , foldOffsetMultiplier(ops[i],ops[j])*args.tsep_pipi);
	corr_comb_j.push_back(coord, data_j.correlator(ops[i],ops[j]).value(t));
	corr_comb_dj.push_back(coord, data_dj.correlator(ops[i],ops[j]).value(t));
	vnm.push_back(nm.str());
      }
    }
  }

  std::cout << "Data in fit:" << std::endl;
  for(int i=0;i<corr_comb_j.size();i++){
    std::cout << vnm[i] << " " << corr_comb_j.coord(i).t << " " << corr_comb_j.value(i) << std::endl;
  }
  if(cmdline.write_fit_data){
    std::cout << "Writing fit data to data_in_fit.hdf5 (and key data_in_fit.key)" << std::endl;
    std::vector<jackknifeDistributionD> fd(corr_comb_j.size());
    std::ofstream of("data_in_fit.key");
    for(int i=0;i<corr_comb_j.size();i++){
      of << i << " " << vnm[i] << " " << (int)corr_comb_j.coord(i).t << std::endl;
      fd[i] = corr_comb_j.value(i);
    }
    writeParamsStandard(fd, "data_in_fit.hdf5");
  }
 
  std::cout << "Performing any data transformations required by the fit func" << std::endl;
  transformData(corr_comb_j, corr_comb_dj, args.t_min, args.t_max, args.fitfunc);
 
  std::cout << "Data post-transformation:" << std::endl;
  for(int i=0;i<corr_comb_j.size();i++){
    std::cout << pmap_descr[corr_comb_j.coord(i).param_map] << " " << corr_comb_j.coord(i).t << " " << corr_comb_j.value(i) << std::endl;
  }

  if(cmdline.load_guess) loadGuess(guess, cmdline.guess_file);

  std::cout << "Performing fit" << std::endl;

  jackknifeDistribution<Params> params(nsample, guess);
  jackknifeDistributionD chisq(nsample), chisq_per_dof(nsample);

  fitOptions opt;
  cmdline.exportOptions(opt);
  args.exportOptions(opt);

  fit(params, chisq, chisq_per_dof,
      corr_comb_j, corr_comb_dj, args.fitfunc, param_map,
      args.nstate, args.Lt, args.t_min, args.t_max, args.correlated, args.Ascale, args.Cscale, opt);

  std::cout << "Params:\n";
  {
    for(int p=0;p<params.sample(0).size();p++){
      std::string tag = params.sample(0).tag(p);
      jackknifeDistributionD tmp;
      standardIOhelper<jackknifeDistributionD, jackknifeDistribution<Params> >::extractStructEntry(tmp, params, p);
      std::cout << tag << " = " << tmp << std::endl;
    }
  }

  double dof = chisq.sample(0)/chisq_per_dof.sample(0);  
  jackknifeDistributionD pvalue(nsample, [&](const int s){ return chiSquareDistribution::pvalue(dof, chisq.sample(s)); });

  std::cout << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  std::cout << "P-value: " << pvalue << std::endl;

#ifdef HAVE_HDF5
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
  writeParamsStandard(pvalue, "pvalue.hdf5");
#endif

  if(args.correlated) analyzeChisq(corr_comb_j, params, args.fitfunc, args.nstate, args.Lt, args.t_min, args.t_max, args.Ascale, args.Cscale, pmap_descr); 

  std::cout << "Done\n";
  return 0;
}


