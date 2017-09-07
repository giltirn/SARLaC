#ifndef _FIT_SIMPLE_MAIN_H_
#define _FIT_SIMPLE_MAIN_H_

template<typename FitFunc, template<typename> class CostFunctionPolicy>
void fitSpecFFcorr(const rawDataCorrelationFunctionD &data, const Args &args, const CMDline &cmdline){
  typedef typename composeFitPolicy<double,FitFunc, standardFitFuncPolicy, CostFunctionPolicy>::type FitPolicies;
  typedef fitter<FitPolicies> Fitter;
  typedef typename FitPolicies::jackknifeFitParameters jackknifeFitParameters;
  
  const int nt_fit = args.t_max - args.t_min + 1;
  const int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc;
  doubleJackknifeCorrelationFunctionD data_dj(nt_fit, [&](const int i){
      typename doubleJackknifeCorrelationFunctionD::ElementType out(args.t_min + i,  doubleJackknifeDistributionD(nsample));
      out.second.resample( data.value(args.t_min + i) );
      return out;
    }
    );

  FitFunc* fitfunc = getFitFunc<FitFunc>(args);
  Fitter fit;
  fit.importFitFunc(*fitfunc);
  
  prepareFitter<FitPolicies,CostFunctionPolicy> prepare(fit, data_dj);

  jackknifeCorrelationFunctionD data_j(nt_fit, [&](const int i){ return typename jackknifeCorrelationFunctionD::ElementType( data_dj.coord(i), data_dj.value(i).toJackknife()); });

  typename FitFunc::ParameterType guess;
  if(cmdline.load_guess){
    parse(guess,cmdline.guess_file);
  }else{
    for(int i=0;i<guess.size();i++)
      guess(i) = 1;
  }

  jackknifeFitParameters params(nsample, guess);
  jackknifeDistributionD chisq(nsample);
  jackknifeDistributionD chisq_per_dof(nsample);

  fit.fit(params, chisq, chisq_per_dof, data_j);

  std::cout << "Params: " << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  
#ifdef HAVE_HDF5
  writeParamsStandard(params, "params.hdf5");
#endif

  //Get all data
  int params_mass_idx = getMassParamIdx(args.fitfunc);
  
  data_j = jackknifeCorrelationFunctionD(args.Lt,
					 [&](const int t){
					   jackknifeDistributionD jack(nsample); jack.resample(data.value(t)); return typename jackknifeCorrelationFunctionD::ElementType(t,std::move(jack));
					 });
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
  plotter.write("effective_mass.py");
  
  delete fitfunc;
}

#endif
