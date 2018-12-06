#include <utils.h>
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

template<typename FitFunc>
void analyzeChisqFF(const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
		    const jackknifeDistribution<Params> &params, const FitFunc &fitfunc,
		    const std::map<std::unordered_map<std::string, std::string> const*, std::string> &pmap_descr){
  struct PP{
    typedef std::map<std::unordered_map<std::string, std::string> const*, std::string> const* PtrType;
    inline static PtrType & descr(){ static PtrType p; return p; }

    inline static void print(std::ostream &os, const SimFitCoordGen &c){ os << "(" << descr()->find(c.param_map)->second << "," << c.t << ")" ; }
  };
  PP::descr() = &pmap_descr;
  
  AnalyzeChisq<FitFunc,PP> chisq_analyze(corr_comb_j, fitfunc, params);
  chisq_analyze.printChisqContribs(Correlation);
  chisq_analyze.examineProjectedDeviationContribsEvalNorm(Correlation);
  chisq_analyze.printChisqContribs(Covariance);
  chisq_analyze.examineProjectedDeviationContribsEvalNorm(Covariance);
}

//nstate is for MultiState variants
void analyzeChisq(const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
		  const jackknifeDistribution<Params> &params, FitFuncType ffunc, const int nstate, const int Lt, double Ascale, double Cscale,
		  const std::map<std::unordered_map<std::string, std::string> const*, std::string> &pmap_descr){
  if(ffunc == FitFuncType::FSimGenOneState){
    typedef FitSimGenOneState FitFunc;
    FitFunc fitfunc(Lt, params.size(), Ascale, Cscale);
    return analyzeChisqFF<FitFunc>(corr_comb_j, params, fitfunc, pmap_descr);
  }else if(ffunc == FitFuncType::FSimGenTwoState){
    typedef FitSimGenTwoState FitFunc;
    FitFunc fitfunc(Lt, params.size(), Ascale, Cscale);
    return analyzeChisqFF<FitFunc>(corr_comb_j, params, fitfunc, pmap_descr);
  }else if(ffunc == FitFuncType::FSimGenThreeState){
    typedef FitSimGenThreeState FitFunc;
    FitFunc fitfunc(Lt, params.size(), Ascale, Cscale);
    return analyzeChisqFF<FitFunc>(corr_comb_j, params, fitfunc, pmap_descr);
  }else if(ffunc == FitFuncType::FSimGenThreeStateLogEdiff){
    typedef FitSimGenThreeStateLogEdiff FitFunc;
    FitFunc fitfunc(Lt, params.size(), Ascale, Cscale);
    return analyzeChisqFF<FitFunc>(corr_comb_j, params, fitfunc, pmap_descr);
  }else if(ffunc == FitFuncType::FSimGenMultiState){
    typedef FitSimGenMultiState FitFunc;
    FitFunc fitfunc(nstate, Lt, params.size(), Ascale, Cscale);
    return analyzeChisqFF<FitFunc>(corr_comb_j, params, fitfunc, pmap_descr);
  }else if(ffunc == FitFuncType::FSimGenMultiStateLogEdiff){
    typedef FitSimGenMultiStateLogEdiff FitFunc;
    FitFunc fitfunc(nstate, Lt, params.size(), Ascale, Cscale);
    return analyzeChisqFF<FitFunc>(corr_comb_j, params, fitfunc, pmap_descr);
  }else{
    assert(0);
  }
}


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
  setupParameterMaps(subfit_pmaps, param_map, guess, args.operators, args.fitfunc, args.Ascale, args.nstate);

  //Write guess template and exit if requested
  if(cmdline.save_guess_template){ saveGuess(guess, "guess_template.dat"); exit(0); }
  
  const std::vector<Operator> &ops = args.operators;

  //Read data
  RawData raw_data;
  if(!cmdline.load_combined_data) raw_data.read(args.Lt, args.data_dir, args.traj_start, args.traj_inc, args.traj_lessthan,
						args.pipi_figure_file_format, args.pipi_bubble_file_format, args.tsep_pipi, args.tstep_pipi,
						args.pipi_to_sigma_file_format, args.tstep_pipi_to_sigma,
						args.sigma2pt_file_format, args.sigma_bubble_file_format,
						ops);

  ResampledData<jackknifeCorrelationFunction> data_j;
  ResampledData<doubleJackCorrelationFunction> data_dj;
  if(cmdline.load_combined_data){
    loadCheckpoint(data_j, data_dj, cmdline.load_combined_data_file);
    for(int i=0;i<ops.size();i++)
      for(int j=i;j<ops.size();j++)
	assert(data_j.haveData(ops[i],ops[j]));
    
  }else{
    data_j.generatedResampledData(raw_data, args.bin_size, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction);
    data_dj.generatedResampledData(raw_data, args.bin_size, args.Lt, args.tsep_pipi, args.do_vacuum_subtraction);
  }  

  if(cmdline.save_combined_data) saveCheckpoint(data_j, data_dj, cmdline.save_combined_data_file);

  int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc/args.bin_size;

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
  
  if(cmdline.load_guess) loadGuess(guess, cmdline.guess_file);

  std::cout << "Performing fit" << std::endl;

  jackknifeDistribution<Params> params(nsample, guess);
  jackknifeDistributionD chisq(nsample), chisq_per_dof(nsample);


  fitOptions opt;
  cmdline.exportOptions(opt);

  fit(params, chisq, chisq_per_dof,
      corr_comb_j, corr_comb_dj, args.fitfunc, param_map,
      args.nstate, args.Lt, args.Ascale, args.Cscale, opt);

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

  analyzeChisq(corr_comb_j, params, args.fitfunc, args.nstate, args.Lt, args.Ascale, args.Cscale, pmap_descr); 

  std::cout << "Done\n";
  return 0;
}


