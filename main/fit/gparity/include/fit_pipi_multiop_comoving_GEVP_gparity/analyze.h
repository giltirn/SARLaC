#ifndef _FIT_PIPI_COMOVING_GEVP_GPARITY_ANALYZE_H
#define _FIT_PIPI_COMOVING_GEVP_GPARITY_ANALYZE_H

#include<type_traits>
#include<fit.h>
#include<distribution.h>

template<typename GEVPsolverType>
void checkResiduals(const GEVPsolverType &gevp, const int t_max, const int nop){
  typedef typename std::decay<decltype(gevp.residuals(0,0)->operator[](0) )>::type DistributionType;
  for(int t0=0;t0<=t_max-1;t0++){
    for(int t=t0;t<=t_max;t++){
      std::vector<DistributionType> const *resid = gevp.residuals(t0,t);
      typedef iterate<DistributionType> it_t;
      if(resid != NULL){
	for(int i=0;i<nop;i++){
	  bool fail = false;
	  const DistributionType &d = (*resid)[i];
	  for(int s=0;s<it_t::size(d);s++) if(it_t::at(s,d) > 1e-3){ fail = true; break; }
	  if(fail)
	    std::cout << "WARNING: Large residuals found for t0="<<t0 <<" t="<<t << " state=" << i << " : " << d << std::endl;
	}
      }
    }
  }
}

template<typename GEVPsolverType, typename DistributionType>
void computeEnergies(const GEVPsolverType &gevp, const correlationFunction<double, NumericSquareMatrix<DistributionType> > &C, 
		     const int t_max, const int nop){
  std::cout << "Energies:" << std::endl;
  for(int t0=0;t0<=t_max-1;t0++){
    for(int t=t0;t<=t_max;t++){
      auto E = gevp.effectiveEnergy(t0,t);
      if(E.size() > 0){
  	std::cout << t0 << " " << t;
  	for(int n=0;n<nop;n++)
  	  std::cout << " " << E[n];
  	std::cout << std::endl;
	
  	std::ostringstream pth; pth << "effective_energies/" << t0 << "/" << t;
  	createDirectory(pth.str());
  	pth << "/E.hdf5";
  	writeParamsStandard(E,pth.str());
      }
    }
  }
}

template<typename GEVPsolverType, typename DistributionType>
void computeAmplitudes(const GEVPsolverType &gevp, const correlationFunction<double, NumericSquareMatrix<DistributionType> > &C,
		       const int t_max, const int nop, const double Ascale){
  std::cout << "Amplitudes (outer index is operator, inner is state):\n";
  for(int t0=0; t0<=t_max-1; t0++){
    for(int t=t0+1; t<=t_max; t++){
      std::vector<std::vector<DistributionType> > Coeffs_all = gevp. effectiveAmplitude(t0,t,C);     
      std::cout << t0 << " " << t << std::endl;
      for(int op=0;op<nop;op++){
	for(int state=0;state<nop;state++){
	  if(Coeffs_all.size() != 0){
	    Coeffs_all[op][state] = Coeffs_all[op][state]/sqrt(Ascale);
	    std::cout << Coeffs_all[op][state] << " ";
	  }else{
	    std::cout << "- ";
	  }
	}
	std::cout << std::endl;
      }
      if(Coeffs_all.size() != 0){
	std::ostringstream pth; pth << "effective_amplitudes/" << t0 << "/" << t;
	createDirectory(pth.str());
	pth << "/A.hdf5";
	writeParamsStandard(Coeffs_all,pth.str());
      }
    }
  }
}

template<typename GEVPsolverType>
void fitConstantFullRangeFrozenCov(const GEVPsolverType &gevp,
				   const int nop, const int fit_tmin, const int fit_tmax, const int fit_t0min, const int fit_t0max){
  typedef typename std::decay<decltype(gevp.effectiveEnergy(0,0)[0] )>::type DistributionType;
  typedef iterate<DistributionType> it_t;

  std::cout << "Fitting all data within supplied fit ranges that exists and has t0>=t/2" << std::endl;
  std::vector< correlationFunction<double,  DistributionType> > fitdata(nop);
    
  for(int t0=fit_t0min; t0<=fit_t0max; t0++){
    for(int t=fit_tmin; t<=fit_tmax; t++){
      if(t0 >= t/2){
	auto E = gevp.effectiveEnergy(t0,t);
	if(E.size() > 0){	
	  for(int n=0;n<nop;n++){
	    if(!std::isnan(E[n].mean()) && !std::isnan(E[n].standardError())){
	      std::cout << "State " << n << " including " << t0 << " " << t << " with value " << E[n] << std::endl;
	      fitdata[n].push_back(t0, E[n]); //note the coordinate is irrelevant because we fit to a constant
	    }
	  }
	}
      }
    }
  }
  
  auto dist_init = fitdata[0].value(0).getInitializer();
    
  typedef FitFuncLinearMultiDim<double,double,0> FitFunc;    
  for(int n=0;n<nop;n++){
    std::cout << "Doing frozen correlated fit for operator " << n << std::endl;

    if(fitdata[n].size() > 0){
      FitFunc fitfunc;

      typedef typename composeFitPolicy<FitFunc, standardFitFuncPolicy, frozenCorrelatedFitPolicy, MarquardtLevenbergMinimizerPolicy, DistributionType>::type FitPolicies;
      typename fitter<FitPolicies>::minimizerParamsType mlparams;
      mlparams.verbose = true;
      mlparams.lambda_max = 1e10;
      mlparams.lambda_factor = 3;
      fitter<FitPolicies> fit(mlparams);
      fit.importFitFunc(fitfunc);

      const correlationFunction<double,  DistributionType> &d = fitdata[n];
      importCostFunctionParameters<frozenCorrelatedFitPolicy, FitPolicies> import(fit, d);
	
      parameterVector<double> guess(1, fitdata[n].value(0).mean());

      typedef typename DistributionType::template rebase<parameterVector<double> > ParameterDistributionType;

      ParameterDistributionType params(dist_init,guess);
	
      DistributionType chisq(dist_init), chisq_per_dof(dist_init);
	
      fit.fit(params, chisq, chisq_per_dof, fitdata[n]);

      double dof = chisq.sample(0)/chisq_per_dof.sample(0);
      DistributionType pvalue(dist_init);
      for(int s=0;s<it_t::size(pvalue);s++) it_t::at(s, pvalue)=  chiSquareDistribution::pvalue(dof, it_t::at(s, chisq));

      writeParamsStandard(params, stringize("fit_params_all_state%d.hdf5", n));
      writeParamsStandard(chisq, stringize("fit_chisq_all_state%d.hdf5", n));
      writeParamsStandard(chisq_per_dof, stringize("fit_chisq_per_dof_state%d.hdf5", n));
      writeParamsStandard(pvalue, stringize("fit_pvalue_state%d.hdf5", n));

      std::cout << "Params: " << params << std::endl;
      std::cout << "Chisq: " << chisq << std::endl;
      std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
      std::cout << "Dof: " << dof << std::endl;
      std::cout << "P-value: " << pvalue << std::endl;
    }
  }
}




template<typename GEVPsolverType>
void plotFixedtmt0(const GEVPsolverType &gevp, const int nop, const int t_max, const int Cmin=1, const int Cmax=5){
  typedef typename std::decay<decltype(gevp.effectiveEnergy(0,0)[0] )>::type DistributionType;

  //If t-t0=C  where C is a constant, and we restrict  t0 >= t/2   then any t0 >= 2C is valid
  for(int C=Cmin; C<=Cmax; C++){
    std::cout << "Examining energies with t-t0=" << C << " and t0 >= t/2   i.e. t0 >= 2C = " << 2*C << "\n";

    std::vector< correlationFunction<double,  DistributionType> > state_t0dep(nop);

    for(int t0=2*C; t0 < t_max; t0++){
      int t = C + t0;
      
      auto E = gevp.effectiveEnergy(t0,t);
      if(E.size() > 0){	
	std::cout << t0 << " " << t;
	for(int n=0;n<nop;n++){
	  std::cout << " " << E[n];	
	  if(!std::isnan(E[n].mean()) && !std::isnan(E[n].standardError())){
	    state_t0dep[n].push_back(t0, E[n]);
	  }
	}
	std::cout << std::endl;
      }
    }

    MatPlotLibScriptGenerate plot;
    typedef DataSeriesAccessor< correlationFunction<double,  DistributionType>, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<DistributionType> > accessor;
    for(int n=0;n<nop;n++){
      if(state_t0dep[n].size() > 0){
	plot.plotData(accessor(state_t0dep[n]));
      }
    }
    std::string stub = stringize("plot_fixedtmt0%d", C);
    plot.write(stub + ".py", stub + ".pdf");    
  }
}


template<typename GEVPsolverType>
void fitFixedtmt0frozenCov(const GEVPsolverType &gevp, const int nop, const int fit_t0min, const int fit_t0max, const int Cmin=1, const int Cmax=5){
  typedef typename std::decay<decltype(gevp.effectiveEnergy(0,0)[0] )>::type DistributionType;
  typedef iterate<DistributionType> it_t;

  //If t-t0=C  where C is a constant, and we restrict  t0 >= t/2   then any t0 >= 2C is valid
  for(int C=Cmin; C<=Cmax; C++){
    std::cout << "Fitting energies with t-t0=" << C << " and t0 >= t/2   i.e. t0 >= 2C = " << 2*C << "\n";

    std::vector< correlationFunction<double,  DistributionType> > data_inrange_n(nop);

    for(int t0=fit_t0min; t0 <= fit_t0max; t0++){
      int t = C + t0;
      auto E = gevp.effectiveEnergy(t0,t);
      if(E.size() > 0){	
	for(int n=0;n<nop;n++){
	  if(!std::isnan(E[n].mean()) && !std::isnan(E[n].standardError())){
	    data_inrange_n[n].push_back(t0, E[n]);
	  }
	}
      }
    }

    auto dist_init = data_inrange_n[0].value(0).getInitializer();


    //Do a frozen correlated fit to each operator
    typedef FitFuncLinearMultiDim<double,double,0> FitFunc;    
    for(int n=0;n<nop;n++){
      std::cout << "Doing frozen correlated fit for operator " << n << " for C= "<< C << std::endl;

      const correlationFunction<double,  DistributionType> &data_inrange = data_inrange_n[n];
  
      if(data_inrange.size() > 0){
	FitFunc fitfunc;

	typedef typename composeFitPolicy<FitFunc, standardFitFuncPolicy, frozenCorrelatedFitPolicy, MarquardtLevenbergMinimizerPolicy, DistributionType>::type FitPolicies;

	typename fitter<FitPolicies>::minimizerParamsType mlparams;
	mlparams.verbose = false;
	fitter<FitPolicies> fit(mlparams);
	fit.importFitFunc(fitfunc);

	importCostFunctionParameters<frozenCorrelatedFitPolicy, FitPolicies> import(fit, data_inrange);

	parameterVector<double> guess(1, 0.3);

	typedef typename DistributionType::template rebase<parameterVector<double> > ParameterDistributionType;

	ParameterDistributionType params(dist_init,guess);
	
	DistributionType chisq(dist_init), chisq_per_dof(dist_init);
	
	fit.fit(params, chisq, chisq_per_dof, data_inrange);

	double dof = chisq.sample(0)/chisq_per_dof.sample(0);
	DistributionType pvalue(dist_init);
 	for(int s=0;s<it_t::size(pvalue);s++) it_t::at(s, pvalue) = chiSquareDistribution::pvalue(dof, it_t::at(s,chisq));

	writeParamsStandard(params, stringize("fit_params_fixedtmt0%d_state%d.hdf5", C,n));
	writeParamsStandard(chisq, stringize("fit_chisq_fixedtmt0%d_state%d.hdf5", C,n));
	writeParamsStandard(chisq_per_dof, stringize("fit_chisq_per_dof_fixedtmt0%d_state%d.hdf5", C,n));
	writeParamsStandard(pvalue, stringize("fit_pvalue_fixedtmt0%d_state%d.hdf5", C,n));

	std::cout << "Params: " << params << std::endl;
	std::cout << "Chisq: " << chisq << std::endl;
	std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
	std::cout << "Dof: " << dof << std::endl;
	std::cout << "P-value: " << pvalue << std::endl;
      }
    }
  }
}




template<typename GEVPsolverType, typename GEVPsolverTypeDJ>
void fitConstantFullRange(const GEVPsolverType &gevp, const GEVPsolverTypeDJ &gevp_dj,
			  const int nop, const int fit_tmin, const int fit_tmax, const int fit_t0min, const int fit_t0max){
  typedef typename std::decay<decltype(gevp.effectiveEnergy(0,0)[0] )>::type DistributionType;
  typedef typename std::decay<decltype(gevp_dj.effectiveEnergy(0,0)[0] )>::type CovDistributionType;
  typedef iterate<DistributionType> it_t;
  

  std::cout << "Fitting all data within supplied fit ranges that exists and has t0>=t/2" << std::endl;
  std::vector< correlationFunction<double,  DistributionType> > fitdata_j(nop);
  std::vector< correlationFunction<double,  CovDistributionType> > fitdata_dj(nop);
    
  for(int t0=fit_t0min; t0<=fit_t0max; t0++){
    for(int t=fit_tmin; t<=fit_tmax; t++){
      if(t0 >= t/2){
	auto E_j = gevp.effectiveEnergy(t0,t);
	auto E_dj = gevp_dj.effectiveEnergy(t0,t);

	if(E_j.size() > 0 && E_dj.size() > 0){	
	  for(int n=0;n<nop;n++){
	    if(!std::isnan(E_j[n].mean()) && !std::isnan(E_j[n].standardError())){
	      std::cout << "State " << n << " including " << t0 << " " << t << " with value " << E_j[n] << std::endl;
	      fitdata_j[n].push_back(t0, E_j[n]); //note the coordinate is irrelevant because we fit to a constant
	      fitdata_dj[n].push_back(t0, E_dj[n]);
	    }
	  }
	}
      }
    }
  }
    
  auto dist_init = fitdata_j[0].value(0).getInitializer();


  typedef FitFuncLinearMultiDim<double,double,0> FitFunc;    
  for(int n=0;n<nop;n++){
    std::cout << "Doing frozen correlated fit for operator " << n << std::endl;

    if(fitdata_j[n].size() > 0){
      FitFunc fitfunc;

      typedef typename composeFitPolicy<FitFunc, standardFitFuncPolicy, frozenCorrelatedFitPolicy, MarquardtLevenbergMinimizerPolicy, DistributionType>::type FitPolicies;
      typename fitter<FitPolicies>::minimizerParamsType mlparams;
      mlparams.verbose = true;
      mlparams.lambda_max = 1e10;
      mlparams.lambda_factor = 3;
      fitter<FitPolicies> fit(mlparams);
      fit.importFitFunc(fitfunc);

      importCostFunctionParameters<correlatedFitPolicy, FitPolicies> import(fit, fitdata_dj[n]);
	
      parameterVector<double> guess(1, fitdata_j[n].value(0).mean());

      typedef typename DistributionType::template rebase<parameterVector<double> > ParameterDistributionType;

      ParameterDistributionType params(dist_init,guess);
	
      DistributionType chisq(dist_init), chisq_per_dof(dist_init);
	
      fit.fit(params, chisq, chisq_per_dof, fitdata_j[n]);

      double dof = chisq.sample(0)/chisq_per_dof.sample(0);
      DistributionType pvalue(dist_init);
      for(int s=0;s<it_t::size(pvalue);s++) it_t::at(s, pvalue) = chiSquareDistribution::pvalue(dof, it_t::at(s,chisq));

      writeParamsStandard(params, stringize("fit_params_all_state%d.hdf5", n));
      writeParamsStandard(chisq, stringize("fit_chisq_all_state%d.hdf5", n));
      writeParamsStandard(chisq_per_dof, stringize("fit_chisq_per_dof_state%d.hdf5", n));
      writeParamsStandard(pvalue, stringize("fit_pvalue_state%d.hdf5", n));

      std::cout << "Params: " << params << std::endl;
      std::cout << "Chisq: " << chisq << std::endl;
      std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
      std::cout << "Dof: " << dof << std::endl;
      std::cout << "P-value: " << pvalue << std::endl;
    }
  }
}




//energies, chisq, dof returned for each state {0..nop-1}
template<typename GEVPsolverType, typename GEVPsolverTypeDJ, typename DistributionType>
void fitEigenvaluesNExpFixedt0(std::vector<DistributionType> &energies,
			       std::vector<DistributionType> &chisq, 
			       std::vector<int> &dof,
			       const GEVPsolverType &gevp, const GEVPsolverTypeDJ &gevp_dj,
			       const int nop, const int fit_tmin, const int fit_tmax, const int t0, 
			       const MarquardtLevenbergParameters<double> &mlp){
  typedef typename GEVPsolverTypeDJ::DistributionType CovDistributionType;
  typedef iterate<DistributionType> it_t;
  energies.resize(nop);
  chisq.resize(nop);
  dof.resize(nop);

  bool verbose = mlp.verbose;

  if(verbose) std::cout << "Fitting eigenvalues" << std::endl;
  std::vector< correlationFunction<double,  DistributionType> > fitdata_j(nop);
  std::vector< correlationFunction<double,  CovDistributionType> > fitdata_dj(nop);
    
  bool error = false;
  for(int t=fit_tmin; t<=fit_tmax; t++){
    std::vector<DistributionType> const* evals_j = gevp.evals(t0, t);
    std::vector<CovDistributionType> const* evals_dj = gevp_dj.evals(t0, t);
    if(evals_j == NULL){
      std::cout << "ERROR: (base) eigenvalue vector not available for t0=" << t0 << " t=" << t << std::endl; error = true;
    }
    if(evals_dj == NULL){
      std::cout << "ERROR: (cov) eigenvalue vector not available for t0=" << t0 << " t=" << t << std::endl; error = true;
    }
    if(!error){
      for(int n=0;n<nop;n++){
	fitdata_j[n].push_back(t-t0, (*evals_j)[n]);
	fitdata_dj[n].push_back(t-t0, (*evals_dj)[n]);
      }
    }
  }
  if(error) exit(-1);

  if(verbose){
    std::cout << fitdata_j[0].size() << "data in fit per state. t0=" << t0 << ":" << std::endl;
    for(int n=0;n<nop;n++){
      std::cout << "State " << n << ": " << std::endl;
      for(int tt=0;tt<fitdata_j[n].size();tt++)
	std::cout << "t-t0=" << fitdata_j[n].coord(tt) << " value=" << fitdata_j[n].value(tt) << std::endl;
    }
  }

  auto dist_init = fitdata_j[0].value(0).getInitializer();
  DistributionType one(dist_init,1.);
  typedef FitNStateExp FitFunc;
  
  parameterVector<double> guess({1,0.5});
  FitFunc fitfunc(1); //1 state
  genericFitFuncWrapper<FitFunc> fwrap(fitfunc, guess);
  simpleFitWrapper<DistributionType> fitter(fwrap, MinimizerType::MarquardtLevenberg, generalContainer(mlp)  );
  fitter.freeze(0, one); //freeze A[0] = 1
  
  for(int n=0;n<nop;n++){
    chisq[n] = one;
    
    if(verbose) std::cout << "Doing frozen correlated fit for operator " << n << std::endl;

    if(fitdata_j[n].size() > 0){
      fitter.generateCovarianceMatrix(fitdata_dj[n]);
	
      typedef typename DistributionType::template rebase<parameterVector<double> > ParameterDistributionType;

      ParameterDistributionType params(dist_init,guess);
	
      DistributionType chisq_per_dof(dist_init);
	
      fitter.fit(params, chisq[n], chisq_per_dof, dof[n], fitdata_j[n]);

      energies[n] = distributionStructPeek(params, 1); //extract energy

      if(verbose){
	std::cout << "State: " << n << std::endl;
	std::cout << "Energy: " << energies[n] << std::endl;
	std::cout << "Chisq: " << chisq[n] << std::endl;
	std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
	std::cout << "Dof: " << dof[n] << std::endl;
      }
    }
  }
}
template<typename GEVPsolverType, typename GEVPsolverTypeDJ, typename DistributionType>
void fitEigenvaluesNExpFixedt0(std::vector<DistributionType> &energies,
			       std::vector<DistributionType> &chisq, 
			       std::vector<int> &dof,
			       const GEVPsolverType &gevp, const GEVPsolverTypeDJ &gevp_dj,
			       const int nop, const int fit_tmin, const int fit_tmax, const int t0, bool verbose = true){
  MarquardtLevenbergParameters<double> mlp;
  mlp.verbose = verbose;
  fitEigenvaluesNExpFixedt0(energies, chisq, dof, gevp, gevp_dj, nop, fit_tmin, fit_tmax, t0, mlp);
}




template<typename GEVPsolverType, typename DistributionType>
void analyze_GEVP(const GEVPsolverType &gevp,
		  const correlationFunction<double, NumericSquareMatrix<DistributionType> > &C,
		  const int t_max, const int fit_tmin, const int fit_tmax, const int fit_t0min, const int fit_t0max, const double Ascale){

  const int nop = C.value(0).size();

  checkResiduals(gevp, t_max, nop);
  computeEnergies(gevp, C, t_max, nop);
  computeAmplitudes(gevp, C, t_max, nop, Ascale);

  fitConstantFullRangeFrozenCov(gevp,nop,fit_tmin,fit_tmax,fit_t0min,fit_t0max);

  plotFixedtmt0(gevp,nop,t_max,1,5);

  fitFixedtmt0frozenCov(gevp,nop,fit_t0min,fit_t0max,1,5);		       
}




template<typename GEVPsolverType, typename GEVPsolverTypeDJ, typename DistributionType>
void analyze_GEVP(const GEVPsolverType &gevp, const GEVPsolverTypeDJ &gevp_dj,
		  const correlationFunction<double, NumericSquareMatrix<DistributionType> > &C,
		  const int t_max, const int fit_tmin, const int fit_tmax, const int fit_t0min, const int fit_t0max, const double Ascale){

  const int nop = C.value(0).size();

  checkResiduals(gevp, t_max, nop);
  computeEnergies(gevp, C, t_max, nop);
  computeAmplitudes(gevp, C, t_max, nop, Ascale);

  fitConstantFullRange(gevp,gevp_dj,nop,fit_tmin,fit_tmax,fit_t0min,fit_t0max);

  plotFixedtmt0(gevp,nop,t_max,1,5);
}


#endif
