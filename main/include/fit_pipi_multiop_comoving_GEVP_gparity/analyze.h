#ifndef _FIT_PIPI_COMOVING_GEVP_GPARITY_ANALYZE_H
#define _FIT_PIPI_COMOVING_GEVP_GPARITY_ANALYZE_H

#include<fit/GEVP.h>

void analze_GEVP(const correlationFunction<double, NumericSquareMatrix<jackknifeDistributionD> > &C,
		 const int t_max, const int fit_tmin, const int fit_tmax, const double Ascale,
		 bool verbose_solver){

  const int nop = C.value(0).size();
  const int nsample = C.value(0)(0,0).size();

  //Perform GEVP
  std::cout << "Solving GEVP" << std::endl;
  GEVPsolver<jackknifeDistributionD> gevp(verbose_solver);
  gevp.solve(C, t_max);

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
    
  std::cout << "Amplitudes (outer index is operator, inner is state):\n";
  for(int t0=0; t0<=t_max-1; t0++){
    for(int t=t0+1; t<=t_max; t++){
      std::vector<std::vector<jackknifeDistributionD> > Coeffs_all = gevp. effectiveAmplitude(t0,t,C);     
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

  //If t-t0=C  where C is a constant, and we restrict  t0 >= t/2   then any t0 >= 2C is valid
  for(int C=1; C<5; C++){
    std::cout << "Examining energies with t-t0=" << C << " and t0 >= t/2   i.e. t0 >= 2C = " << 2*C << "\n";

    std::vector< correlationFunction<double,  jackknifeDistributionD> > state_t0dep(nop);

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
    typedef DataSeriesAccessor< correlationFunction<double,  jackknifeDistributionD>, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > accessor;
    for(int n=0;n<nop;n++){
      if(state_t0dep[n].size() > 0){
	plot.plotData(accessor(state_t0dep[n]));
      }
    }
    std::string stub = stringize("plot_fixedtmt0%d", C);
    plot.write(stub + ".py", stub + ".pdf");    

    //Do a frozen correlated fit to each operator
    typedef FitFuncLinearMultiDim<double,double,0> FitFunc;    
    for(int n=0;n<nop;n++){
      std::cout << "Doing frozen correlated fit for operator " << n << " for C= "<< C << std::endl;

      correlationFunction<double,  jackknifeDistributionD> data_inrange;
      for(int i=0;i<state_t0dep[n].size();i++){
	int t = (int)state_t0dep[n].coord(i);
	if(t >= fit_tmin && t <= fit_tmax){
	  data_inrange.push_back(t, state_t0dep[n].value(i));
	}
      }
      if(data_inrange.size() > 0){
	FitFunc fitfunc;

	typedef typename composeFitPolicy<FitFunc, standardFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
	typename fitter<FitPolicies>::minimizerParamsType mlparams;
	mlparams.verbose = false;
	fitter<FitPolicies> fit(mlparams);
	fit.importFitFunc(fitfunc);

	int N = data_inrange.size();
	NumericSquareMatrix<jackknifeDistributionD> cov(N);
	for(int i=0;i<N;i++)
	  for(int j=i;j<N;j++){
	    double cv = jackknifeDistributionD::covariance(data_inrange.value(i), data_inrange.value(j));
	    cov(i,j) = cov(j,i) = jackknifeDistributionD(nsample,cv);
	  }

	NumericVector<jackknifeDistributionD> sigma(N);
	for(int i=0;i<N;i++)
	  sigma(i) = sqrt(cov(i,i));

	NumericSquareMatrix<jackknifeDistributionD> corr(N);
	for(int i=0;i<N;i++)
	  for(int j=i;j<N;j++)
	    corr(i,j) = corr(j,i) = cov(i,j)/sigma(i)/sigma(j);
	
	NumericSquareMatrix<jackknifeDistributionD> inv_corr(corr);
	svd_inverse(inv_corr, corr);
	
	fit.importCostFunctionParameters(inv_corr,sigma);
	
	parameterVector<double> guess(1,0.3);

	jackknifeDistribution<parameterVector<double> > params(nsample,guess);
	
	jackknifeDistributionD chisq(nsample);
	jackknifeDistributionD chisq_per_dof(nsample);
	
	fit.fit(params, chisq, chisq_per_dof, data_inrange);

	double dof = chisq.sample(0)/chisq_per_dof.sample(0);
	jackknifeDistributionD pvalue(nsample, [&](const int s){ return chiSquareDistribution::pvalue(dof, chisq.sample(s)); });

	writeParamsStandard(params, stringize("fit_params_fixedtmt0%d", C));
	writeParamsStandard(chisq, stringize("fit_chisq_fixedtmt0%d", C));
	writeParamsStandard(chisq_per_dof, stringize("fit_chisq_per_dof_fixedtmt0%d", C));
	writeParamsStandard(pvalue, stringize("fit_pvalue_fixedtmt0%d", C));

	std::cout << "Params: " << params << std::endl;
	std::cout << "Chisq: " << chisq << std::endl;
	std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
	std::cout << "Dof: " << dof << std::endl;
	std::cout << "P-value: " << pvalue << std::endl;
      }
    }
  }

}
#endif
