#ifndef _FIT_SIMPLE_FIT_H_
#define _FIT_SIMPLE_FIT_H_

template<typename FitFunc, typename ArgsType>
void plotEffectiveMass(const ArgsType &args, const FitFunc &fitfunc, const jackknifeCorrelationFunctionD &data_j, 
		       const jackknifeDistribution<typename FitFunc::ParameterType> &params, const int params_mass_idx){
  jackknifeCorrelationFunctionD effmass = effectiveMass2pt<jackknifeCorrelationFunctionD,FitFunc>(data_j,fitfunc,params.sample(0), params_mass_idx, args.Lt);
  
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
    jackknifeDistribution<StandardFitParams> p1exp(params.size(),[&](const int s){ return StandardFitParams(params.sample(s)(0), params.sample(s)(1)); });
    FitCosh fcosh(args.Lt);
    plotEffectiveMass(args,fcosh,data_j,p1exp,1);
  }
};


template<typename FitFunc, template<typename> class CostFunctionPolicy, typename ArgsType, typename CMDlineType>
void fitSpecFFcorr(const jackknifeCorrelationFunctionD &data_j, const doubleJackknifeCorrelationFunctionD &data_dj, const ArgsType &args, const CMDlineType &cmdline);

template<typename FitFunc, typename ArgsType, typename CMDlineType>
inline void fitSpecFF(const jackknifeCorrelationFunctionD &data_j, const doubleJackknifeCorrelationFunctionD &data_dj, const ArgsType &args, const CMDlineType &cmdline){
  return args.correlated ? 
    fitSpecFFcorr<FitFunc,correlatedFitPolicy,ArgsType,CMDlineType>(data_j,data_dj,args,cmdline) : 
    fitSpecFFcorr<FitFunc,uncorrelatedFitPolicy,ArgsType,CMDlineType>(data_j, data_dj,args,cmdline);
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
