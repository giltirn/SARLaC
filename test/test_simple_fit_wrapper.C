#include<random.h>
#include<fit.h>
#include<common.h>

using namespace CPSfit;

int main(const int argc, const char** argv){
  RNG.initialize(1234);
  
  int nsample =100;
  int Lt = 10;

  std::cout << "Generating raw data" << std::endl;
  correlationFunction<double, rawDataDistributionD> raw(Lt);
  
  double A = 1.0;
  double m = 0.3;

  for(int t=0;t<Lt;t++){
    rawDataDistributionD mult(nsample);
    gaussianRandom(mult, 1.0, 0.5);

    raw.coord(t) = t;
    raw.value(t) = A*exp(-m*t)*mult;
  }
  
  std::cout << "Resampling data" << std::endl;
  correlationFunction<double, jackknifeDistributionD> data_j(Lt);
  correlationFunction<double, doubleJackknifeDistributionD> data_dj(Lt);
  
  for(int t=0;t<Lt;t++){
    data_j.coord(t) = data_dj.coord(t) = raw.coord(t);
    data_j.value(t).resample(raw.value(t));
    data_dj.value(t).resample(raw.value(t));
  }
  
  typedef FitExp FitFunc;
  FitFunc fitfunc;

  simpleFitFuncWrapper<FitFunc> fwrap(fitfunc); //make the fit function conform the the standard interface
  
  MarquardtLevenbergParameters<double> mlp;
  mlp.verbose = true;
  mlp.delta_cost_min = 1e-10;

  std::cout << "Performing fit with ML" << std::endl;
  simpleFitWrapper fitter(fwrap, MinimizerType::MarquardtLevenberg, mlp);

  //For testing add a frozen fit parameter and a gaussian prior
  fitter.freeze({0}, { jackknifeDistributionD(nsample, 0.7) });
  fitter.addPrior(0.2, 0.001, 1);

  //Generate the covariance matrix
  fitter.generateCovarianceMatrix(data_dj);

  //Fit with standard ML
  StandardFitParams guess(0.9,0.5);
	  
  jackknifeDistribution<StandardFitParams> params(nsample,guess);
  jackknifeDistributionD chisq(nsample), chisq_per_dof(nsample);
  int dof;

  fitter.fit(params, chisq, chisq_per_dof, dof, data_j);

  //Fit with GSL
  GSLtrsMinimizerParams gslp;
  gslp.verbose = true;
  gslp.stopping_conditions[0] = {GSLtrsStoppingCondition::StopCostDelta, 1e-10};
  fitter.setMinimizer(MinimizerType::GSLtrs, gslp);
  
  jackknifeDistribution<StandardFitParams> params_gsl(nsample,guess);
  jackknifeDistributionD chisq_gsl(nsample), chisq_per_dof_gsl(nsample);
  int dof_gsl;

  fitter.fit(params_gsl, chisq_gsl, chisq_per_dof_gsl, dof_gsl, data_j);
  
  std::cout << "ML Params = " << params << std::endl;
  std::cout << "ML Chisq = " << chisq << std::endl;
  std::cout << "ML Chisq/dof = " << chisq_per_dof << std::endl;

  std::cout << "GSL Params = " << params_gsl << std::endl;
  std::cout << "GSL Chisq = " << chisq_gsl << std::endl;
  std::cout << "GSL Chisq/dof = " << chisq_per_dof_gsl << std::endl;

  

  std::vector< std::pair<jackknifeDistributionD,jackknifeDistributionD> > comp = {  { distributionStructPeek(params, &StandardFitParams::A),
										      distributionStructPeek(params_gsl, &StandardFitParams::A) },
										    { distributionStructPeek(params, &StandardFitParams::m),
										      distributionStructPeek(params_gsl, &StandardFitParams::m) },
										    { chisq, chisq_gsl },
										    { chisq_per_dof, chisq_per_dof_gsl } };
  std::vector<std::string> descr = { "A", "m", "chisq", "chisq/dof" };

  std::cout << "Relative differences:\n";
  double tol = 1e-4;
  for(int i=0;i<comp.size();i++){
    jackknifeDistributionD reldiff = 2*(comp[i].first - comp[i].second)/(comp[i].first + comp[i].second);
    std::cout << descr[i] << " " << reldiff << std::endl;
    for(int s=0;s<nsample;s++) if(fabs(reldiff.sample(s)) > tol) error_exit(std::cout << descr[i] << " reldiff sample " << s << " value " << fabs(reldiff.sample(s)) << " is greater than tolerance " << tol << std::endl);
    
  }

  std::cout << "Done\n";
  return 0;
}


