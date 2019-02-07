#ifndef _FIT_SIMPLE_FIT_H_
#define _FIT_SIMPLE_FIT_H_

template<typename FitFunc, typename ArgsType>
void plotEffectiveMass(const ArgsType &args, const FitFunc &fitfunc, const jackknifeCorrelationFunctionD &data_j, 
		       const jackknifeDistribution<typename FitFunc::ParameterType> &params, const int params_mass_idx){
  double guess = params.best()(params_mass_idx);
  jackknifeCorrelationFunctionD effmass = effectiveMass2pt<jackknifeCorrelationFunctionD,FitFunc>(data_j,fitfunc,params.sample(0), params_mass_idx, args.Lt, guess);
  
  MatPlotLibScriptGenerate plotter;
  typedef MatPlotLibScriptGenerate::handleType Handle;
  typename MatPlotLibScriptGenerate::kwargsType plot_args;
    
  typedef DataSeriesAccessor<jackknifeCorrelationFunctionD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > Accessor;
  Accessor a(effmass);
  Handle ah = plotter.plotData(a);

  //   Plot the fitted mass as constant
  typename FitFunc::ParameterType mn = params.best();
  typename FitFunc::ParameterType err = params.standardError();
  const double m = mn(params_mass_idx);
  const double dm = err(params_mass_idx);
  
  std::vector<double> x = {double(args.t_min), double(args.t_max)};
  std::vector<double> upper = {m + dm, m + dm};
  std::vector<double> lower = {m - dm, m - dm};    
  BandVectorAccessor band(x,upper,lower);
  plot_args["alpha"] = 0.2;
  ah = plotter.errorBand(band, plot_args);
  
  
  plotter.createLegend();
  plotter.setXlabel("$t$");
  plotter.setYlabel("$m_{\\rm eff}(t)$");
  plotter.setXaxisBounds(-0.2,args.Lt+0.2);

  const double ymid = m;
  const double yw = 20 * dm;
  
  plotter.setYaxisBounds(ymid-yw, ymid+yw);

  std::cout << "Writing plot to 'effective_mass.py'\n";  
  plotter.write("effective_mass.py", "effective_mass.pdf");
}


template<typename TwoStateFitFunc, typename OneStateFitFunc, typename ArgsType>
void plotTwoStateEffectiveMass(const ArgsType &args, 
			       const jackknifeCorrelationFunctionD &data_j, 
			       const jackknifeDistribution<typename TwoStateFitFunc::ParameterType> &params, 
			       const OneStateFitFunc &fitfunc_1state, const TwoStateFitFunc &fitfunc_2state,
			       const int fitfunc_1state_params_mass_idx){

  //Generate the one-state effective mass from the data
  StandardFitParams base(1.0,1.0);

  int nsample = data_j.value(0).size();
  jackknifeCorrelationFunctionD effmass_data = effectiveMass2pt<jackknifeCorrelationFunctionD,OneStateFitFunc>(data_j,fitfunc_1state,base, fitfunc_1state_params_mass_idx, args.Lt);
  
  //Generate the fit curve from the 2 state fit
  jackknifeCorrelationFunctionD curve_2state(args.Lt);
  for(int t=0;t<args.Lt;t++){
    curve_2state.coord(t) = t;
    curve_2state.value(t) = jackknifeDistributionD(nsample, [&](const int s){ return fitfunc_2state.value(t, params.sample(s)); });
  }

  //Generate one-state effective mass from fit curve
  jackknifeCorrelationFunctionD effmass_curve = effectiveMass2pt<jackknifeCorrelationFunctionD,OneStateFitFunc>(curve_2state,fitfunc_1state,base, fitfunc_1state_params_mass_idx, args.Lt);

 
  //Plot everything
  MatPlotLibScriptGenerate plotter;
  typedef MatPlotLibScriptGenerate::handleType Handle;
  typename MatPlotLibScriptGenerate::kwargsType plot_args;
    
  typedef DataSeriesAccessor<jackknifeCorrelationFunctionD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > Accessor;
  Accessor a(effmass_data);
  Handle ah = plotter.plotData(a);
  
  //   Plot the fitted curve as error band 
  std::vector<double> x(args.t_max - args.t_min + 1);
  std::vector<double> upper(args.t_max - args.t_min + 1);
  std::vector<double> lower(args.t_max - args.t_min + 1);

  for(int t=args.t_min ; t <= args.t_max; t++){
    x[t-args.t_min] = t;

    double y = effmass_curve.value(t).best();
    double dy = effmass_curve.value(t).standardError();

    lower[t-args.t_min] = y-dy;
    upper[t-args.t_min] = y+dy;
  }
  BandVectorAccessor band(x,upper,lower);
  plot_args["alpha"] = 0.2;
  ah = plotter.errorBand(band, plot_args);
    
  plotter.createLegend();
  plotter.setXlabel("$t$");
  plotter.setYlabel("$m_{\\rm eff}(t)$");
  plotter.setXaxisBounds(-0.2,args.Lt+0.2);

  const double yw = 20 * (upper.back() - lower.back())/2.;
  const double ymid = (upper.back() + lower.back())/2;

  plotter.setYaxisBounds(ymid-yw, ymid+yw);

  std::cout << "Writing plot to 'effective_mass.py'\n";  
  plotter.write("effective_mass.py", "effective_mass.pdf");
}



template<typename FitFunc, typename ArgsType>
void plotRaw(const ArgsType &args, const FitFunc &fitfunc, const jackknifeCorrelationFunctionD &data_j, 
	     const jackknifeDistribution<typename FitFunc::ParameterType> &params){
  jackknifeCorrelationFunctionD fit_data(args.Lt,[&](const int t){ 
      jackknifeDistributionD val(params.size(), [&](const int s){ return fitfunc.value(t,params.sample(s)); });
      return jackknifeCorrelationFunctionD::ElementType(t,val);
    });

  MatPlotLibScriptGenerate plotter;
  typedef MatPlotLibScriptGenerate::handleType Handle;
  typename MatPlotLibScriptGenerate::kwargsType plot_args;
    
  typedef DataSeriesAccessor<jackknifeCorrelationFunctionD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > Accessor;
  Accessor a(data_j);
  Handle ah = plotter.plotData(a);
  
  Accessor band(fit_data);
  plot_args["alpha"] = 0.2;
  ah = plotter.errorBand(band, plot_args);
  
  plotter.createLegend();
  plotter.setXlabel("$t$");
  plotter.setYlabel("$y(t)$");
  plotter.setXaxisBounds(-0.2,args.Lt+0.2);

  std::cout << "Writing plot to 'plot.py'\n";  
  plotter.write("plot.py", "plot.pdf");
}

template<typename FitFunc, typename ArgsType>
struct FitFuncPolicy{
  static inline FitFunc* get(const ArgsType &args){ assert(0); }
  static inline void plot(const ArgsType &args, const FitFunc &fitfunc, const jackknifeCorrelationFunctionD &data_j, const jackknifeDistribution<typename FitFunc::ParameterType> &params){ assert(0); }
};

template<typename ArgsType>
struct FitFuncPolicy<FitCosh,ArgsType>{
  typedef FitCosh FitFunc;
  static inline FitCosh* get(const ArgsType &args){ return new FitCosh(args.Lt); }
  static inline void plot(const ArgsType &args, const FitFunc &fitfunc, const jackknifeCorrelationFunctionD &data_j, const jackknifeDistribution<typename FitFunc::ParameterType> &params){
    plotEffectiveMass(args,fitfunc,data_j,params,1);
  }
};

template<typename ArgsType>
struct FitFuncPolicy<FitSinh,ArgsType>{
  typedef FitSinh FitFunc;
  static inline FitSinh* get(const ArgsType &args){ return new FitSinh(args.Lt); }
  static inline void plot(const ArgsType &args, const FitFunc &fitfunc, const jackknifeCorrelationFunctionD &data_j, const jackknifeDistribution<typename FitFunc::ParameterType> &params){
    plotEffectiveMass(args,fitfunc,data_j,params,1);
  }
};

template<typename ArgsType>
struct FitFuncPolicy<FitExp,ArgsType>{
  typedef FitExp FitFunc;
  static inline FitExp* get(const ArgsType &args){ return new FitExp; }
  static inline void plot(const ArgsType &args, const FitFunc &fitfunc, const jackknifeCorrelationFunctionD &data_j, const jackknifeDistribution<typename FitFunc::ParameterType> &params){
    plotEffectiveMass(args,fitfunc,data_j,params,1);
  }
};

template<typename ArgsType>
struct FitFuncPolicy<FitConstant,ArgsType>{
  typedef FitConstant FitFunc;
  static inline FitConstant* get(const ArgsType &args){ return new FitConstant; }
  static inline void plot(const ArgsType &args, const FitFunc &fitfunc, const jackknifeCorrelationFunctionD &data_j, const jackknifeDistribution<typename FitFunc::ParameterType> &params){
    plotRaw(args,fitfunc,data_j,params);
  }
};

template<typename ArgsType>
struct FitFuncPolicy<FitTwoStateCosh,ArgsType>{
  typedef FitTwoStateCosh FitFunc;
  static inline FitTwoStateCosh* get(const ArgsType &args){ return new FitTwoStateCosh(args.Lt); }
  static inline void plot(const ArgsType &args, const FitFunc &fitfunc, const jackknifeCorrelationFunctionD &data_j, const jackknifeDistribution<typename FitFunc::ParameterType> &params){
    //Plot the one-state effective mass
    //jackknifeDistribution<StandardFitParams> p1exp(params.size(),[&](const int s){ return StandardFitParams(params.sample(s)(0), params.sample(s)(1)); });
    //FitCosh fcosh(args.Lt);
    //plotEffectiveMass(args,fcosh,data_j,p1exp,1);

    FitCosh fcosh(args.Lt);
    plotTwoStateEffectiveMass<FitTwoStateCosh, FitCosh, ArgsType>(args, data_j, params, fcosh, fitfunc, 1); 
  }
};


template<typename FitFunc, template<typename> class CostFunctionPolicy, typename ArgsType, typename CMDlineType>
void fitSpecFFcorr(const jackknifeCorrelationFunctionD &data_j, const doubleJackknifeCorrelationFunctionD &data_dj, const ArgsType &args, const CMDlineType &cmdline);

template<typename FitFunc, typename ArgsType, typename CMDlineType>
inline void fitSpecFF(const jackknifeCorrelationFunctionD &data_j, const doubleJackknifeCorrelationFunctionD &data_dj, const ArgsType &args, const CMDlineType &cmdline){
  switch(args.covariance_strategy){
  case CovarianceStrategy::Correlated:
    return fitSpecFFcorr<FitFunc,correlatedFitPolicy,ArgsType,CMDlineType>(data_j,data_dj,args,cmdline);
  case CovarianceStrategy::Uncorrelated:
    return fitSpecFFcorr<FitFunc,uncorrelatedFitPolicy,ArgsType,CMDlineType>(data_j,data_dj,args,cmdline);
  case CovarianceStrategy::FrozenCorrelated:
    return fitSpecFFcorr<FitFunc,frozenCorrelatedFitPolicy,ArgsType,CMDlineType>(data_j,data_dj,args,cmdline);
  default:
    error_exit(std::cout << "fitSpecFF unknown CovarianceStrategy " << args.covariance_strategy << std::endl);
  }
};

template<typename ArgsType, typename CMDlineType>
inline void fitResampled(const jackknifeCorrelationFunctionD &data_j, const doubleJackknifeCorrelationFunctionD &data_dj, const ArgsType &args, const CMDlineType &cmdline){
  switch(args.fitfunc){
  case FitFuncType::FCosh:
    return fitSpecFF<FitCosh,ArgsType,CMDlineType>(data_j, data_dj,args,cmdline);
  case FitFuncType::FSinh:
    return fitSpecFF<FitSinh,ArgsType,CMDlineType>(data_j, data_dj,args,cmdline);
  case FitFuncType::FExp:
    return fitSpecFF<FitExp,ArgsType,CMDlineType>(data_j, data_dj,args,cmdline);
  case FitFuncType::FConstant:
    return fitSpecFF<FitConstant,ArgsType,CMDlineType>(data_j, data_dj,args,cmdline); 
  case FitFuncType::FTwoStateCosh:
    return fitSpecFF<FitTwoStateCosh,ArgsType,CMDlineType>(data_j, data_dj,args,cmdline); 
  default:
    error_exit(std::cout << "fit: Invalid fitfunc " << args.fitfunc << std::endl);
  };
}
template<typename ArgsType, typename CMDlineType>
inline void fit(const jackknifeCorrelationFunctionD &data_j, const doubleJackknifeCorrelationFunctionD &data_dj, const ArgsType &args, const CMDlineType &cmdline){
  if(cmdline.save_combined_data){
#ifdef HAVE_HDF5
    std::cout << "Writing resampled data to " << cmdline.save_combined_data_file << std::endl;
    HDF5writer writer(cmdline.save_combined_data_file);
    write(writer, data_j, "data_j");
    write(writer, data_dj, "data_dj");
#else
    error_exit("fitSpecFFcorr: Saving amplitude data requires HDF5\n");
#endif
  }
  fitResampled<ArgsType,CMDlineType>(data_j, data_dj, args, cmdline);
}
  
#endif
