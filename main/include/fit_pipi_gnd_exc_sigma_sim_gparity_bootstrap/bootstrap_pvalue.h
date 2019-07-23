#ifndef _PIPI_GND_EXC_SIM_FIT_BOOTSTRAP_PVALUE_H
#define _PIPI_GND_EXC_SIM_FIT_BOOTSTRAP_PVALUE_H

#include<config.h>
#include<utils/macros.h>
#include<fit/bootstrap_pvalue.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity_bootstrap/fit.h>

CPSFIT_START_NAMESPACE

void bootstrapPvalue(const double q2,
		     const correlationFunction<SimFitCoordGen,  bootstrapDistributionD> &corr_comb_b,
		     const correlationFunction<SimFitCoordGen,  bootJackknifeDistributionD> &corr_comb_bj,
		     const bootstrapDistribution<Params> &base_params,
		     FitFuncType ffunc, const std::unordered_map<std::string,size_t> &param_map,
		     const int nstate, const int Lt, const int t_min, const int t_max,
		     const bool correlated, const CovarianceMatrix covariance_matrix,
		     const double Ascale, const double Cscale, const fitOptions &fopt){
  
  std::unique_ptr<genericFitFuncBase> fitfunc = getFitFunc(ffunc, nstate, t_min, Lt, param_map.size(), Ascale, Cscale, base_params.best());

  //Recenter
  correlationFunction<SimFitCoordGen,  bootstrapDistributionD> corr_comb_b_recentered(corr_comb_b);

  for(int i=0;i<corr_comb_b.size();i++){
    auto c = genericFitFuncBase::getWrappedCoord(corr_comb_b.coord(i));
    auto p = genericFitFuncBase::getWrappedParams(base_params.best());
    double fval = fitfunc->value(c, p);
    double correction = fval - corr_comb_b.value(i).best();
    
    corr_comb_b_recentered.value(i) = corr_comb_b_recentered.value(i) + correction;
  }
  
  //Fit recentered data
  bootstrapDistribution<Params> rparams(base_params);
  bootstrapDistributionD chisq(rparams.getInitializer()), chisq_per_dof(rparams.getInitializer());

  fit(rparams, chisq, chisq_per_dof,
      corr_comb_b_recentered, corr_comb_bj, ffunc, param_map,
      nstate, Lt, t_min, t_max, correlated, covariance_matrix, Ascale, Cscale, fopt);


  //Compute p-value
  int nboot = rparams.size();

  //Note we shouldn't use the "best" value here because it is not obtained using a bootstrap resampling; instead use the mean
  std::cout << "Bootstrap distribution of q^2 has mean " << chisq.mean() << " and std.dev " << chisq.standardDeviation() << std::endl;

  std::vector<double> q2_dist(nboot); for(int i=0;i<nboot;i++) q2_dist[i] = chisq.sample(i);
  std::sort(q2_dist.begin(), q2_dist.end(), [&](const double a, const double b){ return a<b; });

  double p_boot = computePvalue(q2, q2_dist);
  
  {
    std::ofstream of("p_boot.dat");
    of << p_boot << std::endl;
  }

  //Save sorted bootstrap q2 distribution as rawDataDistribution
  rawDataDistributionD q2_dist_r(nboot, [&](const int i){ return q2_dist[i]; });
  writeParamsStandard(q2_dist_r, "boot_q2.hdf5");
  
  std::cout << "Bootstrap p-value " << p_boot << std::endl;
}


CPSFIT_END_NAMESPACE

#endif
