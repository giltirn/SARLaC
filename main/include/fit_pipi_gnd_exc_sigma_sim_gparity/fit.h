#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_FIT_H
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_FIT_H


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
  }else if(ffunc == FitFuncType::FSimGenMultiStateCparam){
    typedef FitSimGenMultiStateCparam FitFunc;
    FitFunc fitfunc(nstate, Lt, params.size(), Ascale, Cscale);
    return analyzeChisqFF<FitFunc>(corr_comb_j, params, fitfunc, pmap_descr);
  }else{
    assert(0);
  }
}

#endif
