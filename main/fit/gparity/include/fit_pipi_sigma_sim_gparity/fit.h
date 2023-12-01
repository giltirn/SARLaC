#ifndef _PIPI_SIGMA_SIM_FIT_FIT_H_
#define _PIPI_SIGMA_SIM_FIT_FIT_H_

#include<config.h>
#include<utils/macros.h>
#include<pipi_common/analyze_chisq.h>

CPSFIT_START_NAMESPACE

struct SimFitArgs{
  bool correlated;
  int Lt;
  int tsep_pipi;
  double Ascale;
  double Cscale;

  bool load_guess;
  std::string guess_file;
  
  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  bool write_covariance_matrix;
  std::string write_covariance_matrix_file;

  bool ledoit_wolf;

  bool load_mlparams;
  std::string mlparams_file;

SimFitArgs(): correlated(true), Lt(64), tsep_pipi(4), Ascale(1e13), Cscale(1e13), load_guess(false), load_frozen_fit_params(false), write_covariance_matrix(false), ledoit_wolf(false), load_mlparams(false){}
};


template<typename importer>
struct doLedoitWolf{ inline static double doit(importer &im){ return 0.; }  };
template<typename FitPolicies>
struct doLedoitWolf<importCostFunctionParameters<correlatedFitPolicy,FitPolicies> >{ 
  inline static double doit(importCostFunctionParameters<correlatedFitPolicy,FitPolicies> &im){ return im.LedoitWolfShrinkage(); }
};


template<template<typename> class corrUncorrFitPolicy>
void fit_corr_uncorr(const simFitCorrFuncJ &data_j, const simFitCorrFuncDJ &data_dj, const SimFitArgs &args){
  typedef FitSim FitFunc;

  typedef typename FitFunc::Params Params;
  Params guess;
  if(args.load_guess){
    parse(guess, args.guess_file);
    std::cout << "Loaded guess: " << guess << std::endl;
  }
  const int nsample = data_dj.value(0).size();
  
  FitFunc fitfunc(args.Lt, args.tsep_pipi, args.Ascale, args.Cscale);

  typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, corrUncorrFitPolicy>::type FitPolicies;
  fitter<FitPolicies> fitter;

  if(args.load_mlparams){
    typename fitter<FitPolicies>::minimizerParamsType minparams;
    parse(minparams, args.mlparams_file);
    std::cout << "Loaded minimizer params: " << minparams << std::endl;
    fitter.setMinimizerParams(minparams);
  }

  fitter.importFitFunc(fitfunc);

  if(args.load_frozen_fit_params)
    readFrozenParams(fitter, args.load_frozen_fit_params_file, nsample);
  
  importCostFunctionParameters<corrUncorrFitPolicy,FitPolicies> prepare(fitter, data_dj);

  double ledoit_wolf_lambda = 0.;
  if(args.ledoit_wolf) ledoit_wolf_lambda = doLedoitWolf<decltype(prepare)>::doit(prepare);
    
  if(args.write_covariance_matrix) prepare.writeCovarianceMatrixHDF5(args.write_covariance_matrix_file);

  jackknifeDistribution<Params> params(nsample, guess);
  jackknifeDistributionD chisq;
  jackknifeDistributionD chisq_per_dof;
  fitter.fit(params, chisq, chisq_per_dof, data_j);

  double dof = chisq.sample(0)/chisq_per_dof.sample(0);
  
  jackknifeDistributionD pvalue(nsample, [&](const int s){ return chiSquareDistribution::pvalue(dof, chisq.sample(s)); });
  
  distributionPrint<jackknifeDistribution<Params> >::printer(new pipiParamsPrinter<FitFunc>);

  std::cout << "Params: " << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  std::cout << "Dof: " << dof << std::endl;
  std::cout << "P-value: " << pvalue << std::endl;

#ifdef HAVE_HDF5
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
  writeParamsStandard(pvalue, "pvalue.hdf5");
#endif

  struct PP{
    inline static void print(std::ostream &os, const SimFitCoord &c){ os << "(" << c.type << "," << c.t << ")" ; }
    inline static std::string typeInfo(const SimFitCoord &c){ std::ostringstream os; os << c.type; return os.str(); }   
  };

  AnalyzeChisq<FitFunc,PP> chisq_analyze(data_j, fitfunc, params, ledoit_wolf_lambda);

  if(args.ledoit_wolf){
    chisq_analyze.printChisqContribs(CorrelationLedoitWolf);
    chisq_analyze.examineProjectedDeviationContribsEvalNorm(CorrelationLedoitWolf);
  }
  chisq_analyze.printChisqContribs(Correlation);
  chisq_analyze.examineProjectedDeviationContribsEvalNorm(Correlation);

  chisq_analyze.printChisqContribs(Covariance);
  chisq_analyze.examineProjectedDeviationContribsEvalNorm(Covariance);
}

void fit(const simFitCorrFuncJ &data_j, const simFitCorrFuncDJ &data_dj, const SimFitArgs &args){
  if(args.correlated){
    fit_corr_uncorr<correlatedFitPolicy>(data_j,data_dj,args);
  }else{
    fit_corr_uncorr<uncorrelatedFitPolicy>(data_j,data_dj,args);
  }
}

CPSFIT_END_NAMESPACE

#endif
