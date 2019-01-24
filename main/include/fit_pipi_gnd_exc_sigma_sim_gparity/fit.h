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
void transformData(correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
		   correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> &corr_comb_dj,
		   const int t_min, const int t_max,
		   FitFuncType ffunc){
  if(ffunc == FitFuncType::FSimGenMultiStateSub || ffunc == FitFuncType::FSimGenMultiStateTminSub || ffunc == FitFuncType::FSimGenMultiStateTminSubForceZero){
    typedef std::unordered_map<std::string, std::string> InnerParamMap;
    std::map< std::pair<InnerParamMap const*, int>, int > data_map;
    for(int d=0;d<corr_comb_j.size();d++){
      auto const &c = corr_comb_j.coord(d);
      data_map[ {c.param_map, (int)c.t} ] = d;
    }
    correlationFunction<SimFitCoordGen,  jackknifeDistributionD> corr_comb_j_out;
    correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> corr_comb_dj_out;

    if(ffunc == FitFuncType::FSimGenMultiStateSub){
      for(int d=0;d<corr_comb_j.size();d++){
	auto const &c = corr_comb_j.coord(d);
	auto it = data_map.find( {c.param_map, (int)c.t+1 } );
	if(it != data_map.end()){ //data point at t_max doesn't have anything we can subtract so drop it
	  int dp = it->second;
	  corr_comb_j_out.push_back(c, jackknifeDistributionD(corr_comb_j.value(dp) - corr_comb_j.value(d)));
	  corr_comb_dj_out.push_back(c, doubleJackknifeDistributionD(corr_comb_dj.value(dp) - corr_comb_dj.value(d)));
	}
      }
    }else if(ffunc == FitFuncType::FSimGenMultiStateTminSub || ffunc == FitFuncType::FSimGenMultiStateTminSubForceZero){
      for(int d=0;d<corr_comb_j.size();d++){
	auto const &c = corr_comb_j.coord(d);
	auto itmin = data_map.find( {c.param_map, t_min } );
	assert(itmin != data_map.end());
	int dtmin = itmin->second;
	if(dtmin != d){ //don't include the data point we are subtracting because it is exactly zero; the covariance matrix becomes singular 
	  jackknifeDistributionD subj = corr_comb_j.value(dtmin);
	  doubleJackknifeDistributionD subdj = corr_comb_dj.value(dtmin);

	  corr_comb_j_out.push_back(c, jackknifeDistributionD(corr_comb_j.value(d) - subj));
	  corr_comb_dj_out.push_back(c, doubleJackknifeDistributionD(corr_comb_dj.value(d) - subdj));
	}
	else if(ffunc == FitFuncType::FSimGenMultiStateTminSubForceZero){
	  //Here we force the fit to go through zero at t_min by adding a zero with a random distribution of some small width
	  if(!RNG.isInitialized()) RNG.initialize(1234);
	  double w = 1e-5;
	  int nsample = corr_comb_j.value(dtmin).size();
	  jackknifeDistributionD zeroj(nsample);
	  doubleJackknifeDistributionD zerodj(nsample);
	  gaussianRandom(zeroj, 0., w/sqrt(nsample-1)); //jackknife error is sqrt(N-1)x std.dev
	  for(int i=0;i<nsample;i++) gaussianRandom(zerodj.sample(i), 0., w/sqrt(nsample-1)); 
	  
	  corr_comb_j_out.push_back(c, zeroj);
	  corr_comb_dj_out.push_back(c, zerodj);
	}

      }
    }

    corr_comb_j = std::move(corr_comb_j_out);
    corr_comb_dj = std::move(corr_comb_dj_out);    
  }
}

#endif



// else if(ffunc == FitFuncType::FSimGenMultiStateTminSub){
//       int tmin = getTmin(corr_comb_j);
//       for(int d=0;d<corr_comb_j.size();d++){
// 	auto const &c = corr_comb_j.coord(d);
// 	auto itmin = data_map.find( {c.param_map, tmin } );
// 	auto itminp1  = data_map.find( {c.param_map, tmin+1 } );
// 	assert(itmin != data_map.end() && itminp1 != data_map.end());

// 	int dtmin = itmin->second, dtminp1 = itminp1->second;

// 	jackknifeDistributionD subj = (corr_comb_j.value(dtmin) + corr_comb_j.value(dtminp1))/2.;
// 	doubleJackknifeDistributionD subdj = (corr_comb_dj.value(dtmin) + corr_comb_dj.value(dtminp1))/2.;

// 	corr_comb_j_out.push_back(c, jackknifeDistributionD(corr_comb_j.value(d) - subj));
// 	corr_comb_dj_out.push_back(c, doubleJackknifeDistributionD(corr_comb_dj.value(d) - subdj));
//       }
//     }
