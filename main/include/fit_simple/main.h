#ifndef _FIT_SIMPLE_MAIN_H_
#define _FIT_SIMPLE_MAIN_H_

template<typename FitFunc, template<typename> class CostFunctionPolicy>
void fitSpecFFcorr(const jackknifeCorrelationFunctionD &data_j, const doubleJackknifeCorrelationFunctionD &data_dj, const Args &args, const CMDline &cmdline){
  typedef typename composeFitPolicy<FitFunc, standardFitFuncPolicy, CostFunctionPolicy>::type FitPolicies;
  typedef fitter<FitPolicies> Fitter;
  typedef typename FitPolicies::FitParameterDistribution FitParameterDistribution;
  
  const int nt_fit = args.t_max - args.t_min + 1;

  doubleJackknifeCorrelationFunctionD data_dj_inrange(nt_fit, 
						      [&](const int i){ return typename doubleJackknifeCorrelationFunctionD::ElementType(args.t_min + i, data_dj.value(args.t_min + i)); }
						      );
  
  jackknifeCorrelationFunctionD data_j_inrange(nt_fit, 
					       [&](const int i){ return typename jackknifeCorrelationFunctionD::ElementType(args.t_min + i, data_j.value(args.t_min + i)); }
					       );

  FitFunc* fitfunc = getFitFunc<FitFunc>(args);
  Fitter fit;
  fit.importFitFunc(*fitfunc);

  importCostFunctionParameters<CostFunctionPolicy,FitPolicies> prepare(fit,data_dj_inrange);

  typename FitFunc::ParameterType guess;
  if(cmdline.load_guess){
    parse(guess,cmdline.guess_file);
  }else{
    for(int i=0;i<guess.size();i++)
      guess(i) = 1;
  }

  const int nsample = data_j.value(0).size();
  FitParameterDistribution params(nsample, guess);
  jackknifeDistributionD chisq(nsample);
  jackknifeDistributionD chisq_per_dof(nsample);

  fit.fit(params, chisq, chisq_per_dof, data_j_inrange);

  std::cout << "Params: " << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  
#ifdef HAVE_HDF5
  writeParamsStandard(params, "params.hdf5");
#endif

  //Get all data
  int params_mass_idx = getMassParamIdx(args.fitfunc);
  
  jackknifeCorrelationFunctionD effmass = effectiveMass2pt<jackknifeCorrelationFunctionD,FitFunc>(data_j,*fitfunc,params.sample(0), params_mass_idx, args.Lt);
  
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
  
  delete fitfunc;
}

#endif
