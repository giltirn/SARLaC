#ifndef FIT_SIMPLE_BOOTSTRAP_BOOTSTRAP_PVALUE_H_
#define FIT_SIMPLE_BOOTSTRAP_BOOTSTRAP_PVALUE_H_

#include<fit/bootstrap_pvalue.h>

void bootstrapPvalue(const bootstrapCorrelationFunctionD &data_b,
		     const bootJackknifeCorrelationFunctionD &data_bj,
		     const bootstrapDistribution<parameterVectorD> &params,
		     const bootstrapDistributionD &chisq,
		     const Args &args, const CMDline &cmdline){

  std::cout << "Computing bootstrap p-value" << std::endl;
      
  //Recenter the data
  std::unique_ptr< FitFuncManagerBase > ffman = getFitFuncManager(args.fitfunc, args.Lt, args.t_min, args.t_max);
  const genericFitFuncBase &fitfunc = *ffman->getFitFunc();
  
  correlationFunctionD corrections(args.Lt, 
				   [&](const int t){ 
				     double fit_val = fitfunc.value( genericFitFuncBase::getWrappedCoord((double)t), params.best() ); 
				     double data_val = data_b.value(t).best();
				     return correlationFunctionD::ElementType(t, fit_val - data_val);
				   });
  bootstrapCorrelationFunctionD data_b_recentered(data_b);
  recenter(data_b_recentered, corrections);
  
  //Do the fit to the recentered data
  bootstrapDistribution<parameterVectorD> params_rc;
  bootstrapDistributionD chisq_rc;
  int dof;
  fit(params_rc, chisq_rc, dof, data_b_recentered, data_bj, args, cmdline, false);

  //Note we shouldn't use the "best" value here because it is not obtained using a bootstrap resampling; instead use the mean
  std::cout << "Bootstrap distribution of q^2 has mean " << chisq_rc.mean() << " and std.dev " << chisq_rc.standardDeviation() << std::endl;

  std::vector<double> q2_dist(args.nboot); for(int i=0;i<args.nboot;i++) q2_dist[i] = chisq_rc.sample(i);
  std::sort(q2_dist.begin(), q2_dist.end(), [&](const double a, const double b){ return a<b; });

  double p_boot = computePvalue(chisq.best(), q2_dist);
  
  {
    std::ofstream of("p_boot.dat");
    of << p_boot << std::endl;
  }

  //Save sorted bootstrap q2 distribution as rawDataDistribution
  rawDataDistributionD q2_dist_r(args.nboot, [&](const int i){ return q2_dist[i]; });
  writeParamsStandard(q2_dist_r, "boot_q2.hdf5");
  
  std::cout << "Bootstrap p-value " << p_boot << std::endl;
}

#endif
