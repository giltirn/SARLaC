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
    inline static std::string typeInfo(const SimFitCoordGen &c){ return descr()->find(c.param_map)->second; }
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
		  const jackknifeDistribution<Params> &params, FitFuncType ffunc, const int nstate, const int Lt, 
		  const int t_min, const int t_max,
		  double Ascale, double Cscale,
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
  }else if(ffunc == FitFuncType::FSimGenMultiStateSub){
    typedef FitSimGenMultiStateSub FitFunc;
    FitFunc fitfunc(nstate, Lt, params.size(), Ascale, Cscale);
    return analyzeChisqFF<FitFunc>(corr_comb_j, params, fitfunc, pmap_descr);
  }else if(ffunc == FitFuncType::FSimGenMultiStateTminSub || ffunc == FitFuncType::FSimGenMultiStateTminSubForceZero){
    typedef FitSimGenMultiStateTminSub FitFunc;
    FitFunc fitfunc(nstate, Lt, params.size(), t_min, Ascale, Cscale);
    return analyzeChisqFF<FitFunc>(corr_comb_j, params, fitfunc, pmap_descr);
  }else{
    assert(0);
  }
}



//Implement fit-function specific data transformations
template<typename DistributionType>
void transformData(correlationFunction<SimFitCoordGen,  DistributionType> &corr_comb,
		   const int t_min, const int t_max,
		   FitFuncType ffunc){
  if(ffunc == FitFuncType::FSimGenMultiStateSub || ffunc == FitFuncType::FSimGenMultiStateTminSub || ffunc == FitFuncType::FSimGenMultiStateTminSubForceZero){
    typedef std::unordered_map<std::string, std::string> InnerParamMap;
    std::map< std::pair<InnerParamMap const*, int>, int > data_map;
    for(int d=0;d<corr_comb.size();d++){
      auto const &c = corr_comb.coord(d);
      data_map[ {c.param_map, (int)c.t} ] = d;
    }
    correlationFunction<SimFitCoordGen,  DistributionType> corr_comb_out;

    if(ffunc == FitFuncType::FSimGenMultiStateSub){
      for(int d=0;d<corr_comb.size();d++){
	auto const &c = corr_comb.coord(d);
	auto it = data_map.find( {c.param_map, (int)c.t+1 } );
	if(it != data_map.end()) //data point at t_max doesn't have anything we can subtract so drop it
	  corr_comb_out.push_back(c, DistributionType(corr_comb.value(it->second) - corr_comb.value(d)));	
      }
    }else if(ffunc == FitFuncType::FSimGenMultiStateTminSub || ffunc == FitFuncType::FSimGenMultiStateTminSubForceZero){
      for(int d=0;d<corr_comb.size();d++){
	auto const &c = corr_comb.coord(d);
	auto itmin = data_map.find( {c.param_map, t_min } );
	assert(itmin != data_map.end());
	int dtmin = itmin->second;
	if(dtmin != d){ //don't include the data point we are subtracting because it is exactly zero; the covariance matrix becomes singular 
	  DistributionType sub = corr_comb.value(dtmin);
	  corr_comb_out.push_back(c, DistributionType(corr_comb.value(d) - sub));
	}
	else if(ffunc == FitFuncType::FSimGenMultiStateTminSubForceZero){
	  //Here we force the fit to go through zero at t_min by adding a zero with a random distribution of some small width
	  if(!RNG.isInitialized()) RNG.initialize(1234);
	  double w = 1e-5;
	  int nsample = corr_comb.value(dtmin).size();
	  DistributionType zero = corr_comb.value(dtmin);
	  zeroit(zero);

	  for(int i=0;i<iterate<DistributionType>::size(zero);i++)
	    gaussianRandom(iterate<DistributionType>::at(i,zero), 0., w/sqrt(nsample-1)); 
	  
	  corr_comb_out.push_back(c, zero);
	}

      }
    }

    corr_comb = std::move(corr_comb_out);
  }
}

#endif

