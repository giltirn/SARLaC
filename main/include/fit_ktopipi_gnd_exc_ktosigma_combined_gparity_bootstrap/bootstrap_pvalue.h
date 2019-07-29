#ifndef _FIT_KTOSIGMA_KTOPIPI_GPARITY_BOOTSTRAP_BOOTSTRAP_PVALUE_H_
#define _FIT_KTOSIGMA_KTOPIPI_GPARITY_BOOTSTRAP_BOOTSTRAP_PVALUE_H_

#include<config.h>
#include<utils/macros.h>

#include<fit/bootstrap_pvalue.h>

CPSFIT_START_NAMESPACE

void bootstrapPvalue(const std::vector<double> &q2, //one for each q!
		     const ResampledData<bootstrapDistributionD> &data_b,
		     const ResampledData<bootJackknifeDistributionD> &data_bj,
		     const std::vector< bootstrapDistribution<taggedValueContainer<double,std::string> > > &base_params,
		     simultaneousFitBase<bootstrapDistribution> const* fitter,
		     const std::vector<PiPiOperator> &operators,
		     const int Lt, const int tmin_k_op, const int tmin_op_snk, 
		     bool correlated, CovarianceMatrix covariance_matrix){
  
  std::cout << "Computing bootstrap p-values" << std::endl;

  typedef taggedValueContainer<double,std::string> Params;
  typedef iterate<bootstrapDistributionD> iter;

  ResampledDataContainers<bootstrapDistribution> fit_data(data_b, data_bj);

  SimFitDataContainers<bootstrapDistribution> simfit_data;
  fitter->generateSimFitData(simfit_data, fit_data, operators, Lt, tmin_k_op, tmin_op_snk, covariance_matrix);

  typedef correlationFunction<SimFitCoordGen, bootstrapDistributionD> SimFitCorrFuncBoot;
  std::vector<SimFitCorrFuncBoot> & simfit_data_b = simfit_data.A0_sim_b;
  
  int nq = simfit_data_b.size();

  auto const &pmap_descr = fitter->getParameterMapDescr();

  //Recenter
  std::cout << "Recentering" << std::endl;
  for(int q=0;q<nq;q++){
    for(int i=0;i<simfit_data_b[q].size();i++){
      double fitval = fitter->evaluateFitFunc(simfit_data_b[q].coord(i), base_params[q].best());
      double cenval = simfit_data_b[q].value(i).best();
      double correction = fitval - cenval;
      
      std::cout << "Q" << q+1 << " " << printCoord(simfit_data_b[q].coord(i), pmap_descr) << " "
		<< " fit value " << fitval << " central value " << cenval << " correction " << correction << std::endl;

      simfit_data_b[q].value(i) = simfit_data_b[q].value(i) + correction;
    }
  }

  //Fit
  std::cout << "Re-fitting" << std::endl;

  std::vector<bootstrapDistribution<Params> > params;
  std::vector<bootstrapDistributionD> q2_boot;
  
  fitter->fit(params, q2_boot, simfit_data, operators, Lt, tmin_k_op, tmin_op_snk, correlated, covariance_matrix, false);

  //Compute p-value
  std::cout << "Computing p-values" << std::endl;
  int nboot = q2_boot[0].size();

  std::vector<rawDataDistributionD> q2_dist_r(nq, rawDataDistributionD(nboot)); //for output write as rawDataDistribution
  std::vector<double> pvalue(nq);

  for(int q=0;q<nq;q++){
    //Note we shouldn't use the "best" value here because it is not obtained using a bootstrap resampling; instead use the mean
    std::cout << "Q" << q+1 << " bootstrap distribution of q^2 has mean " << q2_boot[q].mean() << " and std.dev " << q2_boot[q].standardDeviation() << std::endl;

    std::vector<double> q2_dist(nboot); for(int i=0;i<nboot;i++) q2_dist[i] = q2_boot[q].sample(i);
    std::sort(q2_dist.begin(), q2_dist.end(), [&](const double a, const double b){ return a<b; });

    pvalue[q] = computePvalue(q2[q], q2_dist);

    for(int s=0;s<nboot;s++) q2_dist_r[q].sample(s) = q2_dist[s];
    
    std::cout << "Q" << q+1 << " q^2 " << q2[q] << ", bootstrap p-value " << pvalue[q] << std::endl;
  }
  
  {
    std::ofstream of("p_boot.dat");
    for(int q=0;q<nq;q++)
      of << pvalue[q] << std::endl;
  }

  //Save sorted bootstrap q2 distribution as rawDataDistribution
  writeParamsStandard(q2_dist_r, "boot_q2.hdf5");
}
 





CPSFIT_END_NAMESPACE

#endif
