#include<random.h>
#include<fit.h>
#include<common.h>

using namespace CPSfit;

struct results{
  jackknifeDistribution<StandardFitParams> params;
  jackknifeDistributionD chisq, chisq_per_dof;
  int dof;
  
  results(const int nsample, const StandardFitParams &guess): params(nsample,guess), chisq(nsample), chisq_per_dof(nsample){}

  static void compare(const results &a, const results &b, const std::string pre, const double tol = 1e-4){
    static std::vector<std::string> descr = { "A", "m", "chisq", "chisq/dof" };
    int nsample = a.chisq.size();

    std::vector< std::pair<jackknifeDistributionD,jackknifeDistributionD> > comp = {  { distributionStructPeek(a.params, &StandardFitParams::A),
											distributionStructPeek(b.params, &StandardFitParams::A) },
										      { distributionStructPeek(a.params, &StandardFitParams::m),
											distributionStructPeek(b.params, &StandardFitParams::m) },
										      { a.chisq, b.chisq },
										      { a.chisq_per_dof, b.chisq_per_dof } };
    
    std::cout << pre << "Relative differences:\n";
    for(int i=0;i<comp.size();i++){
      jackknifeDistributionD reldiff = 2*(comp[i].first - comp[i].second)/(comp[i].first + comp[i].second);
      std::cout << descr[i] << " " << reldiff << std::endl;
      for(int s=0;s<nsample;s++) if(fabs(reldiff.sample(s)) > tol) error_exit(std::cout << descr[i] << " reldiff sample " << s << " value " << fabs(reldiff.sample(s)) << " is greater than tolerance " << tol << std::endl);
      
    }
  }
  
  void print(const std::string pre) const{
    std::cout << pre << "Params = " << params << std::endl;
    std::cout << pre << "Chisq = " << chisq << std::endl;
    std::cout << pre << "Chisq/dof = " << chisq_per_dof << std::endl;
  }
};





int main(const int argc, const char** argv){
  int arg=1;
  while(arg < argc){
    std::string sarg(argv[arg]);
    if(sarg == "-nthread"){
      omp_set_num_threads(strToAny<int>(argv[arg+1]));
      arg+=2;
    }else{
      error_exit(std::cout << "Unrecognized argument: " << sarg << std::endl);
    }
  }

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
  simpleFitWrapper fitter(fwrap, MinimizerType::MarquardtLevenberg);

  //For testing add a frozen fit parameter and a gaussian prior
  fitter.freeze({0}, { jackknifeDistributionD(nsample, 0.7) });
  fitter.addPrior(0.2, 0.001, 1);
  
  //Generate the covariance matrix
  fitter.generateCovarianceMatrix(data_dj);
  
  StandardFitParams guess(0.9,0.5);
    
  //Fit with standard ML
  results results_std(nsample, guess);
  {
    MarquardtLevenbergParameters<double> mlp;
    mlp.verbose = true;
    mlp.delta_cost_min = 1e-10;
    fitter.setMinimizer(MinimizerType::MarquardtLevenberg, mlp);
  
    std::cout << "Performing fit with ML" << std::endl;
    results &r = results_std;
    fitter.fit(r.params, r.chisq, r.chisq_per_dof, r.dof, data_j);
  }
  
  //Fit with GSL trs
  results results_gsl_trs(nsample, guess);
  {
    GSLtrsMinimizerParams gslp;
    gslp.verbose = true;
    gslp.stopping_conditions[0] = {GSLtrsStoppingCondition::StopCostDelta, 1e-10};
    fitter.setMinimizer(MinimizerType::GSLtrs, gslp);
  
    std::cout << "Performing fit with GSL ML implementation" << std::endl;
    results &r = results_gsl_trs;
    fitter.fit(r.params, r.chisq, r.chisq_per_dof, r.dof, data_j);
  }  

  //Fit with GSL multimin fdf minimizer
  results results_gsl_multimin_fdf(nsample, guess);
  {
    GSLmultidimMinimizerParams gslp;
    gslp.verbose = true;
    gslp.algorithm = GSLmultiminAlgorithm::ConjugateGradientFR;
    //gslp.stopping_conditions = { {GSLmultiminStoppingCondition::StopGradient, 1e-3} };
    gslp.stopping_conditions = { {GSLmultiminStoppingCondition::StopRelativeGradients, 1e-3} };
    gslp.line_search_tol = 0.; //0.01;

    fitter.setMinimizer(MinimizerType::GSLmultimin, gslp);

    std::cout << "Performing fit with GSL multimin fdf implementation" << std::endl;
    results &r = results_gsl_multimin_fdf;
    fitter.fit(r.params, r.chisq, r.chisq_per_dof, r.dof, data_j);
  }

  //Fit with GSL multimin f minimizer
  results results_gsl_multimin_f(nsample, guess);
  {
    GSLmultidimMinimizerParams gslp;
    gslp.verbose = true;
    gslp.algorithm = GSLmultiminAlgorithm::NMsimplex2;
    gslp.stopping_conditions = { {GSLmultiminStoppingCondition::StopGradient, 1e-3} };
    gslp.step_size = { 0.01 };

    fitter.setMinimizer(MinimizerType::GSLmultimin, gslp);

    std::cout << "Performing fit with GSL multimin f implementation" << std::endl;
    results &r = results_gsl_multimin_f;
    fitter.fit(r.params, r.chisq, r.chisq_per_dof, r.dof, data_j);
  }


  //Fit with Minuit2
#ifdef HAVE_MINUIT2
  results results_minuit2(nsample, guess);
  {
    Minuit2minimizerParams minp;
    minp.verbose = true;
    minp.tolerance = 1e-10;

    fitter.setMinimizer(MinimizerType::Minuit2, minp);

    std::cout << "Performing fit with Minuit2 implementation" << std::endl;
    results &r = results_minuit2;
    fitter.fit(r.params, r.chisq, r.chisq_per_dof, r.dof, data_j);
  }
#endif

  results_std.print("ML ");
  results_gsl_trs.print("GSL TRS ");
  results_gsl_multimin_fdf.print("GSL multimin fdf ");
  results_gsl_multimin_f.print("GSL multimin f ");
#ifdef HAVE_MINUIT2
  results_minuit2.print("Minuit2 ");
#endif

  results::compare(results_gsl_trs, results_std, "GSL TRS ");
  results::compare(results_gsl_multimin_fdf, results_std, "GSL multimin fdf ");
  results::compare(results_gsl_multimin_f, results_std, "GSL multimin f ");
#ifdef HAVE_MINUIT2
  results::compare(results_minuit2, results_std, "Minuit2 ");
#endif

  std::cout << "Done\n";
  return 0;
}


