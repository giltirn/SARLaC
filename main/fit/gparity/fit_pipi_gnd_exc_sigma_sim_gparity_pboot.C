#include <utils.h>
#include <random.h>
#include <pipi_common/pipi_common.h>
#include <fit/bootstrap_pvalue.h>
using namespace SARLaC;

#include<fit_pipi_gnd_exc_sim_gparity/enums.h>
#include<fit_pipi_gnd_exc_sim_gparity/fit_central.h>

#include<fit_pipi_gnd_exc_sigma_sim_gparity/fitfunc.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/filters.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/raw_data.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/resampled_data.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/fit.h>

#include<fit_pipi_gnd_exc_sigma_sim_gparity_pboot/args.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity_pboot/cmdline.h>

struct RawDataResampler{
  void operator()(RawData &raw, const std::vector<int> &map) const{
    raw.iterateOverRawDistributions([&](auto &dist){ bootstrapResampleRaw(dist,map); });
  }
};

struct PiPiFitter{
  const Args &args;
  const CMDline &cmdline;

  int nsample;
  Filters filters;
  const std::map< std::pair<Operator,Operator>, SubFitFuncParameterMap > &subfit_pmaps;
  const ParamTagIdxMap &param_map;

  Params guess;
  
  fitCentralOptions opt;

  std::unique_ptr<genericFitFuncBase> fitfunc;

  int nboot;
  bootstrapDistribution<Params> boot_params;

  PiPiFitter(const Args &args, const CMDline &cmdline,
	     const std::map< std::pair<Operator,Operator>, SubFitFuncParameterMap > &subfit_pmaps, const ParamTagIdxMap &param_map,
	     const Params &guess, const int nboot): args(args), cmdline(cmdline), subfit_pmaps(subfit_pmaps), 
						    param_map(param_map), guess(guess), boot_params(bootstrapDistribution<Params>::initType(nboot)){
    if(cmdline.load_filters)
      parse(filters,cmdline.load_filters_file);

    cmdline.exportOptions(opt);
    args.exportOptions(opt);

    fitfunc = getFitFunc(args.fitfunc, args.nstate, args.t_min, args.Lt, param_map.size(), args.Ascale, args.Cscale, guess);
  }

  correlationFunction<SimFitCoordGen,  jackknifeDistributionD> getData(const RawData &raw_data) const{
    ResampledData<jackknifeCorrelationFunctionD> data_j;
    basicBinResampler resampler(args.bin_size);
    data_j.generatedResampledData(raw_data, resampler, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction, args.timeslice_avg_vac_sub);
    
    //Find which data satisfy tmin, tmax, filters
    const std::vector<Operator> &ops = args.operators;
    std::vector<DataDescr> keep = getFitDataElemIdx(data_j, ops, args.tsep_pipi, args.t_min, args.t_max, filters, cmdline.load_filters);

    correlationFunction<SimFitCoordGen,  jackknifeDistributionD> corr_comb_j;
    filterData(corr_comb_j, data_j, keep, subfit_pmaps, args.tsep_pipi);

    transformData(corr_comb_j, args.t_min, args.t_max, args.fitfunc);

    return corr_comb_j;
  }
    
  double fit(Params &params, const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j) const{
    params = guess;
    double chisq, chisq_per_dof;

    fitCentral(params, chisq, chisq_per_dof,
	corr_comb_j, fitfunc, param_map,
	args.nstate, args.Lt, args.t_min, args.t_max, args.correlated, args.Ascale, args.Cscale, opt);
    
    std::cout << "Params = "<< params << std::endl;
    std::cout << "Chisq = " << chisq << std::endl;
    std::cout << "Chisq/dof = " << chisq_per_dof << std::endl; 
    return chisq;
  }

  //Called by the function that computes the bootstrap p-value
  double operator()(const RawData &raw_data, const correlationFunction<SimFitCoordGen, double> &corrections, const int b){
    correlationFunction<SimFitCoordGen,  jackknifeDistributionD> corr_comb_j = getData(raw_data);

    assert(corrections.size() == corr_comb_j.size());
    for(int i=0;i<corrections.size();i++){
      assert(corr_comb_j.coord(i) == corrections.coord(i));
      corr_comb_j.value(i) = corr_comb_j.value(i) + corrections.value(i);
    }
    Params params;
    double q2 = fit(params, corr_comb_j);

    boot_params.sample(b) = params;
    return q2;
  }
};



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

  const std::vector<Operator> &ops = args.operators;

  //Setup the fit func parameter maps
  std::map< std::pair<Operator,Operator>, SubFitFuncParameterMap > subfit_pmaps;
  ParamTagIdxMap param_map;
  Params guess;
  setupParameterMaps(subfit_pmaps, param_map, guess, args.operators, args.fitfunc, args.Ascale, args.nstate);

  pmapDescrType pmap_descr =  getPmapDescriptions(ops, subfit_pmaps); //description of the pmaps

  //Write guess template and exit if requested
  if(cmdline.save_guess_template){ saveGuess(guess, "guess_template.dat"); exit(0); }

  //Read data
  RawData raw_data;

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
		  ops, 0);
  }
  if(cmdline.save_raw_data){
    std::cout << "Saving raw data to " << cmdline.save_raw_data_file << std::endl;
    HDF5writer wr(cmdline.save_raw_data_file);  raw_data.write(wr, "raw_data");
  }
  


  if(cmdline.load_guess) loadGuess(guess, cmdline.guess_file);

  std::cout << "Performing fit" << std::endl;
  PiPiFitter fitter(args, cmdline, subfit_pmaps, param_map, guess, args.nboot);
  
  correlationFunction<SimFitCoordGen,  jackknifeDistributionD> orig_data = fitter.getData(raw_data);
  int nsample = orig_data.value(0).size();

  Params orig_params;
  double q2 = fitter.fit(orig_params, orig_data);
  
  fitter.guess = orig_params;

  parameterVector<double> orig_params_p(orig_params.size(), [&](const int i){ return orig_params(i); });
  
  std::cout << "Data mean vs fit:\n";
  int ndata = orig_data.size();
  correlationFunction<SimFitCoordGen,  double> orig_data_means(ndata), orig_data_fitvals(ndata);
  for(int i=0;i<ndata;i++){
    orig_data_means.coord(i) = orig_data_fitvals.coord(i) = orig_data.coord(i);
    orig_data_means.value(i) = orig_data.value(i).mean();
    orig_data_fitvals.value(i) = fitter.fitfunc->value(orig_data.coord(i), orig_params_p);
    
    std::cout << coordDescr(orig_data_means.coord(i), pmap_descr) << " " << orig_data_means.value(i) << " " << orig_data_fitvals.value(i) << " " << orig_data_means.value(i) - orig_data_fitvals.value(i) << std::endl;
  }
  
  std::cout << "Computing bootstrap p-value:" << std::endl;
 
  //Get the resample table
  resampleTableOptions ropt;
  ropt.read_from_file = cmdline.load_boot_resample_table;
  ropt.read_file = cmdline.load_boot_resample_table_file;
  ropt.write_to_file = cmdline.save_boot_resample_table;
  ropt.write_file = cmdline.save_boot_resample_table_file;

  std::vector<std::vector<int> > rtable = generateResampleTable(nsample, args.nboot, BootResampleTableType::NonOverlappingBlock, args.block_size, RNG, ropt);

  std::vector<double> q2_boot_dist;

  double p_boot = bootstrapPvalue(q2, raw_data, nsample, orig_data_means, orig_data_fitvals, RawDataResampler(), fitter, rtable, -1, &q2_boot_dist);

  std::cout << "Bootstrap p-value: " << p_boot << std::endl;

  {
    std::ofstream of("p_boot.dat");
    of << p_boot << std::endl;
  }

  rawDataDistributionD q2_boot_dist_r(args.nboot, [&](const int s){ return q2_boot_dist[s]; });
  std::cout << "Bootstrap distribution mean " << q2_boot_dist_r.mean() << " variance " << q2_boot_dist_r.variance() << std::endl;
  writeParamsStandard(q2_boot_dist_r, "boot_q2.hdf5");

  int nparam = orig_params.size();
  std::vector<bootstrapDistributionD> boot_params(nparam, bootstrapDistributionD(bootstrapDistributionD::initType(args.nboot)));
  for(int p=0;p<nparam;p++){
    for(int b=0;b<args.nboot;b++)
      boot_params[p].sample(b) = fitter.boot_params.sample(b)(p);
    boot_params[p].best() = orig_params(p); //take bootstrap central value from unresampled fit
  }

  writeParamsStandard(boot_params, "params_boot.hdf5");

  std::cout << "Bootstrap params: "<< std::endl;
  for(int p=0;p<nparam;p++)
    std::cout << orig_params.tag(p) << " " << boot_params[p] << std::endl;

  int dof = ndata - nparam;
  std::cout << "Dof: " << dof << std::endl;

  std::cout << "T^2(" << dof << ", " << nsample-1 << ") mean " << TsquareDistribution::mean(dof, nsample-1) << " variance " << TsquareDistribution::variance(dof, nsample-1) << std::endl;

  std::cout << "Done\n";
  return 0;
}


