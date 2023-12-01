#ifndef _FIT_SIMPLE_FIT_PLOT_H_
#define _FIT_SIMPLE_FIT_PLOT_H_

template<typename FitFunc, template<typename,template<typename> class> class DistributionType, template<typename> class V>
void plotEffectiveMass(const FitFunc &fitfunc, 
		       const correlationFunction<double, DistributionType<double,V> > &data_j, 
		       const DistributionType<typename FitFunc::ParameterType, V> &params, const int params_mass_idx,
		       const int Lt, const int t_min, const int t_max){
  double guess = params.best()(params_mass_idx);
  typedef DistributionType<double,V> DistributionD;
  typedef correlationFunction<double, DistributionD> DistributionCorrelationFunctionD;

  DistributionCorrelationFunctionD effmass = effectiveMass2pt<DistributionCorrelationFunctionD,FitFunc>(data_j,fitfunc,params.sample(0), params_mass_idx, Lt, guess);
  
  {
    std::ofstream os("effective_mass.key"); 
    std::vector<DistributionD> em(effmass.size()); 
    for(int t=0;t<effmass.size();t++){
      em[t] = effmass.value(t);
      os << t << " " << effmass.coord(t) << std::endl;
    }
    writeParamsStandard(em, "effective_mass.hdf5");
  }

  MatPlotLibScriptGenerate plotter;
  typedef MatPlotLibScriptGenerate::handleType Handle;
  typename MatPlotLibScriptGenerate::kwargsType plot_args;
    
  typedef DataSeriesAccessor<DistributionCorrelationFunctionD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<DistributionD> > Accessor;
  Accessor a(effmass);
  Handle ah = plotter.plotData(a);

  //   Plot the fitted mass as constant
  typename FitFunc::ParameterType mn = params.best();
  typename FitFunc::ParameterType err = params.standardError();
  const double m = mn(params_mass_idx);
  const double dm = err(params_mass_idx);
  
  std::vector<double> x = {double(t_min), double(t_max)};
  std::vector<double> upper = {m + dm, m + dm};
  std::vector<double> lower = {m - dm, m - dm};    
  BandVectorAccessor band(x,upper,lower);
  plot_args["alpha"] = 0.2;
  ah = plotter.errorBand(band, plot_args);
  
  
  plotter.createLegend();
  plotter.setXlabel("$t$");
  plotter.setYlabel("$m_{\\rm eff}(t)$");
  plotter.setXaxisBounds(-0.2,Lt+0.2);

  const double ymid = m;
  const double yw = 20 * dm;
  
  plotter.setYaxisBounds(ymid-yw, ymid+yw);

  std::cout << "Writing plot to 'effective_mass.py'\n";  
  plotter.write("effective_mass.py", "effective_mass.pdf");
}


template<typename TwoStateFitFunc, typename OneStateFitFunc, template<typename,template<typename> class> class DistributionType, template<typename> class V>
void plotTwoStateEffectiveMass(const correlationFunction<double, DistributionType<double,V> > &data_j, 
			       const DistributionType<typename TwoStateFitFunc::ParameterType, V> &params, 
			       const OneStateFitFunc &fitfunc_1state, const TwoStateFitFunc &fitfunc_2state,
			       const int fitfunc_1state_params_mass_idx,
			       const int Lt, const int t_min, const int t_max){
  typedef DistributionType<double,V> DistributionD;
  typedef DistributionType<typename TwoStateFitFunc::ParameterType, V> DistributionParams;
  typedef correlationFunction<double, DistributionD> DistributionCorrelationFunctionD;

  //Generate the one-state effective mass from the data
  StandardFitParams base(1.0,1.0);

  DistributionCorrelationFunctionD effmass_data = effectiveMass2pt<DistributionCorrelationFunctionD,OneStateFitFunc>(data_j,fitfunc_1state,base, fitfunc_1state_params_mass_idx, Lt);
  
  DistributionD zero(data_j.value(0)); zeroit(zero);
  int niter = iterate<DistributionD>::size(zero);

  //Generate the fit curve from the 2 state fit
  DistributionCorrelationFunctionD curve_2state(Lt);
  for(int t=0;t<Lt;t++){
    curve_2state.coord(t) = t;
    curve_2state.value(t) = zero;
    for(int s=0;s<niter;s++) iterate<DistributionD>::at(s, curve_2state.value(t)) = fitfunc_2state.value(t, iterate<DistributionParams>::at(s, params));
  }

  //Generate one-state effective mass from fit curve
  DistributionCorrelationFunctionD effmass_curve = effectiveMass2pt<DistributionCorrelationFunctionD,OneStateFitFunc>(curve_2state,fitfunc_1state,base, fitfunc_1state_params_mass_idx, Lt);

 
  //Plot everything
  MatPlotLibScriptGenerate plotter;
  typedef MatPlotLibScriptGenerate::handleType Handle;
  typename MatPlotLibScriptGenerate::kwargsType plot_args;
    
  typedef DataSeriesAccessor<DistributionCorrelationFunctionD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<DistributionD> > Accessor;
  Accessor a(effmass_data);
  Handle ah = plotter.plotData(a);
  
  //   Plot the fitted curve as error band 
  std::vector<double> x(t_max - t_min + 1);
  std::vector<double> upper(t_max - t_min + 1);
  std::vector<double> lower(t_max - t_min + 1);

  for(int t=t_min ; t <= t_max; t++){
    x[t-t_min] = t;

    double y = effmass_curve.value(t).best();
    double dy = effmass_curve.value(t).standardError();

    lower[t-t_min] = y-dy;
    upper[t-t_min] = y+dy;
  }
  BandVectorAccessor band(x,upper,lower);
  plot_args["alpha"] = 0.2;
  ah = plotter.errorBand(band, plot_args);
    
  plotter.createLegend();
  plotter.setXlabel("$t$");
  plotter.setYlabel("$m_{\\rm eff}(t)$");
  plotter.setXaxisBounds(-0.2,Lt+0.2);

  const double yw = 20 * (upper.back() - lower.back())/2.;
  const double ymid = (upper.back() + lower.back())/2;

  plotter.setYaxisBounds(ymid-yw, ymid+yw);

  std::cout << "Writing plot to 'effective_mass.py'\n";  
  plotter.write("effective_mass.py", "effective_mass.pdf");
}



template<typename FitFunc, template<typename,template<typename> class> class DistributionType, template<typename> class V>
void plotRaw(const FitFunc &fitfunc, const correlationFunction<double, DistributionType<double,V> > &data_j, 
	     const DistributionType<typename FitFunc::ParameterType,V> &params, const int Lt, const int t_min, const int t_max){

  typedef DistributionType<double,V> DistributionD;
  typedef DistributionType<typename FitFunc::ParameterType, V> DistributionParams;
  typedef correlationFunction<double, DistributionD> DistributionCorrelationFunctionD;


  DistributionD zero(data_j.value(0)); zeroit(zero);
  int niter = iterate<DistributionD>::size(zero);

  DistributionCorrelationFunctionD fit_data(t_max-t_min+1,[&](const int tt){ 
      int t = tt + t_min;
      DistributionD val(zero);    
      for(int s=0;s<niter;s++) iterate<DistributionD>::at(s, val) = fitfunc.value(t, iterate<DistributionParams>::at(s, params));
      return typename DistributionCorrelationFunctionD::ElementType(t,val);
    });

  MatPlotLibScriptGenerate plotter;
  typedef MatPlotLibScriptGenerate::handleType Handle;
  typename MatPlotLibScriptGenerate::kwargsType plot_args;
    
  typedef DataSeriesAccessor<DistributionCorrelationFunctionD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<DistributionD> > Accessor;
  Accessor a(data_j);
  Handle ah = plotter.plotData(a);
  
  Accessor band(fit_data);
  plot_args["alpha"] = 0.2;
  ah = plotter.errorBand(band, plot_args);
  
  plotter.createLegend();
  plotter.setXlabel("$t$");
  plotter.setYlabel("$y(t)$");
  plotter.setXaxisBounds(-0.2,Lt+0.2);

  std::cout << "Writing plot to 'plot.py'\n";  
  plotter.write("plot.py", "plot.pdf");
}

#endif
