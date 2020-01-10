#include <utils.h>
#include <random.h>
#include <pipi_common/pipi_common.h>
#include <pipi_common/analyze_chisq.h>
using namespace CPSfit;

#include<fit_pipi_gnd_exc_sim_gparity/enums.h>
#include<fit_pipi_gnd_exc_sim_gparity/fit.h>

#include<fit_pipi_gnd_exc_sigma_sim_gparity/fitfunc.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/filters.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/raw_data.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/resampled_data.h>

#include<fit_pipi_gnd_exc_sigma_sim_gparity_bootstrap/args.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity_bootstrap/cmdline.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity_bootstrap/fit.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity_bootstrap/bootstrap_pvalue.h>

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

  assert(args.covariance_matrix == CovarianceMatrix::Regular || args.covariance_matrix == CovarianceMatrix::Frozen);

  const std::vector<Operator> &ops = args.operators;

  //Setup the fit func parameter maps
  std::map< std::pair<Operator,Operator>, SubFitFuncParameterMap > subfit_pmaps;
  ParamTagIdxMap param_map;
  Params guess;
  setupParameterMaps(subfit_pmaps, param_map, guess, ops, args.fitfunc, args.Ascale, args.nstate);

  pmapDescrType pmap_descr =  getPmapDescriptions(ops, subfit_pmaps); //description of the pmaps

  //Write guess template and exit if requested
  if(cmdline.save_guess_template){ saveGuess(guess, "guess_template.dat"); exit(0); }

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
		    ops, args.isospin);
    }
    if(cmdline.save_raw_data){
      std::cout << "Saving raw data to " << cmdline.save_raw_data_file << std::endl;
      HDF5writer wr(cmdline.save_raw_data_file);  raw_data.write(wr, "raw_data");
    }

    if(cmdline.remove_samples_in_range)
      raw_data.removeSamplesInRange(cmdline.remove_samples_in_range_start, cmdline.remove_samples_in_range_lessthan);
    if(cmdline.scramble_raw_data)
      raw_data.scrambleSamples();
  }

  ResampledData<bootstrapCorrelationFunctionD> data_b;
  ResampledData<bootJackknifeCorrelationFunctionD> data_bj;
  if(cmdline.load_combined_data){
    HDF5reader rd(cmdline.load_combined_data_file);
    read(rd, data_b, "data_b");
    read(rd, data_bj, "data_bj");
    for(int i=0;i<ops.size();i++)
      for(int j=i;j<ops.size();j++)
	if(!data_b.haveData(ops[i],ops[j]))
	  error_exit(std::cout << "Loaded checkpoint does not contain data for (" << ops[i] << ", " << ops[j] << ")\n");
  }else{
    int nsample = raw_data.correlator(ops[0],ops[0]).value(0).size();

    //Generate the resample table
    resampleTableOptions ropt;
    ropt.read_from_file = cmdline.load_boot_resample_table;
    ropt.read_file = cmdline.load_boot_resample_table_file;
    ropt.write_to_file = cmdline.save_boot_resample_table;
    ropt.write_file = cmdline.save_boot_resample_table_file;
    
    std::vector<std::vector<int> > rtable = generateResampleTable(nsample, args.nboot, args.resample_table_type, args.block_size, RNG, ropt);
    if(rtable.size() != args.nboot) error_exit(std::cout << "Expected resample table of size " << args.nboot << ", got " << rtable.size() << std::endl);

    if(args.isospin != 0 && args.do_vacuum_subtraction) error_exit(std::cout << "Vacuum subtraction is not appropriate to I!=0\n");
    bootstrapBlockResampler resampler(rtable);
    data_b.generatedResampledData(raw_data, resampler, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction, args.timeslice_avg_vac_sub);
    data_bj.generatedResampledData(raw_data, resampler, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction, args.timeslice_avg_vac_sub);
  }  

  if(cmdline.save_combined_data){
    HDF5writer wr(cmdline.save_combined_data_file);
    write(wr, data_b, "data_b");
    write(wr, data_bj, "data_bj");
  }

  //Load extra data filters
  Filters filters;
  if(cmdline.load_filters)
    parse(filters,cmdline.load_filters_file);

  //Find which data satisfy tmin, tmax, filters
  std::vector<DataDescr> keep = getFitDataElemIdx(data_b, ops, args.tsep_pipi, args.t_min, args.t_max, filters, cmdline.load_filters);

  //Add resampled data to full data set with generalized coordinate set appropriately
  correlationFunction<SimFitCoordGen,  bootstrapDistributionD> corr_comb_b;
  correlationFunction<SimFitCoordGen,  bootJackknifeDistributionD> corr_comb_bj;
  filterData(corr_comb_b, data_b, keep, subfit_pmaps, args.tsep_pipi);
  filterData(corr_comb_bj, data_bj, keep, subfit_pmaps, args.tsep_pipi);

  std::cout << "Data in fit:" << std::endl;
  for(int i=0;i<corr_comb_b.size();i++){
    std::cout << coordDescr(corr_comb_b.coord(i),pmap_descr) << " " << corr_comb_b.value(i) 
	      << " (err/mean = " << fabs(corr_comb_b.value(i).standardError()/corr_comb_b.value(i).mean()) << ")" <<   std::endl;
  }
  if(cmdline.write_fit_data){
    std::cout << "Writing fit data to data_in_fit.hdf5 (and key data_in_fit.key)" << std::endl;
    std::vector<bootstrapDistributionD> fd(corr_comb_b.size());
    std::ofstream of("data_in_fit.key");
    for(int i=0;i<corr_comb_b.size();i++){
      of << i << " " << coordDescr(corr_comb_b.coord(i),pmap_descr) << std::endl;
      fd[i] = corr_comb_b.value(i);
    }
    writeParamsStandard(fd, "data_in_fit.hdf5");
  }
 
  std::cout << "Performing any data transformations required by the fit func" << std::endl;
  transformData(corr_comb_b, args.t_min, args.t_max, args.fitfunc);
  transformData(corr_comb_bj, args.t_min, args.t_max, args.fitfunc);

  std::cout << "Data post-transformation:" << std::endl;
  for(int i=0;i<corr_comb_b.size();i++)
    std::cout << coordDescr(corr_comb_b.coord(i),pmap_descr) << " " << corr_comb_b.value(i) << std::endl;

  if(cmdline.load_guess) loadGuess(guess, cmdline.guess_file);

  std::cout << "Performing fit" << std::endl;
  bootstrapInitType binit(args.nboot);

  bootstrapDistribution<Params> params(guess, binit);
  bootstrapDistributionD chisq(binit), chisq_per_dof(binit);

  fitOptions opt;
  cmdline.exportOptions(opt);
  args.exportOptions(opt);

  fit(params, chisq, chisq_per_dof,
      corr_comb_b, corr_comb_bj, args.fitfunc, param_map,
      args.nstate, args.Lt, args.t_min, args.t_max, args.correlated, args.covariance_matrix, args.Ascale, args.Cscale, opt);

  std::cout << "Params:\n";
  {
    for(int p=0;p<params.sample(0).size();p++){
      std::string tag = params.sample(0).tag(p);
      bootstrapDistributionD tmp;
      standardIOhelper<bootstrapDistributionD, bootstrapDistribution<Params> >::extractStructEntry(tmp, params, p);
      std::cout << tag << " = " << tmp << std::endl;
    }
  }

  double dof = chisq.best()/chisq_per_dof.best();  
  bootstrapDistributionD pvalue_chisq(binit);
  for(int s=0;s<iterate<bootstrapDistributionD>::size(pvalue_chisq);s++) 
    iterate<bootstrapDistributionD>::at(s,pvalue_chisq) = chiSquareDistribution::pvalue(dof, iterate<bootstrapDistributionD>::at(s, chisq));

  std::cout << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  std::cout << "P-value(chi^2): " << pvalue_chisq << std::endl;

#ifdef HAVE_HDF5
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
  writeParamsStandard(pvalue_chisq, "pvalue_chisq.hdf5");
#endif

  bootstrapPvalue(chisq.best(), 
		  corr_comb_b, corr_comb_bj, params,
		  args.fitfunc, param_map,
		  args.nstate, args.Lt, args.t_min, args.t_max, 
		  args.correlated, args.covariance_matrix, args.Ascale, args.Cscale, opt);

  std::cout << "Done\n";
  return 0;
}


